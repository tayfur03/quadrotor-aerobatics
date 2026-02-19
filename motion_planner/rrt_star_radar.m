function [path, info] = rrt_star_radar(start, goal, terrain, threat, params)
%RRT_STAR_RADAR Radar-aware RRT* path planner (optimized).
%
% Major performance upgrades:
% 1) O(N log N) nearest/range queries via KD-tree (createns/knnsearch/rangesearch)
% 2) Dynamic rewiring radius r_n = min(r_max, gamma*(log(n)/n)^(1/d))
%    or k-nearest rewiring k_n = ceil(k_const*log(n))
% 3) Vectorized edge collision + edge-cost evaluation on sampled edge matrices
% 4) Preallocated tree buffers + child adjacency list for O(descendants)
% 5) Optional parfor acceleration for neighbor parent/rewire evaluations

 tic;
 if nargin < 5
     params = struct();
 end

 % Parameters
 max_iter = get_param(params, 'max_iter', 3000);
 base_step_size = get_param(params, 'base_step_size', get_param(params, 'step_size', 50));
 max_step_size = get_param(params, 'max_step_size', base_step_size * 3);
 goal_bias = get_param(params, 'goal_bias', 0.1);
 goal_bias_after_goal = get_param(params, 'goal_bias_after_goal', max(0.02, 0.35 * goal_bias));

 rewire_radius_max = get_param(params, 'rewire_radius', base_step_size * 3);
 rewire_radius_min = get_param(params, 'rewire_radius_min', max(base_step_size * 1.2, 1.0));
 rewire_gamma_dyn = get_param(params, 'rewire_gamma_dyn', rewire_radius_max * 2.0);

 spatial_dim = 3;
 rewire_mode_raw = get_param(params, 'rewire_mode', 'dynamic');
 if isstring(rewire_mode_raw)
     rewire_mode_raw = char(rewire_mode_raw);
 end
 if ~ischar(rewire_mode_raw)
     rewire_mode_raw = 'dynamic';
 end
 rewire_mode = lower(rewire_mode_raw);
 if ~ismember(rewire_mode, {'dynamic', 'fixed', 'knearest'})
     warning('RRT*: Unknown rewire_mode="%s". Falling back to "dynamic".', rewire_mode);
     rewire_mode = 'dynamic';
 end
 rewire_k_const = get_param(params, 'rewire_k_const', 2.0 * exp(1) * (1 + 1/spatial_dim));
 rewire_k_min = max(spatial_dim + 1, round(get_param(params, 'rewire_k_min', 16)));
 rewire_k_max = max(rewire_k_min, round(get_param(params, 'rewire_k_max', 256)));

 alpha = get_param(params, 'alpha', 1.2);
 beta = get_param(params, 'beta', 100);
 gamma_pref = get_param(params, 'gamma', 0.3);

 min_clearance = get_param(params, 'min_clearance', 20);
 preferred_agl = get_param(params, 'preferred_agl', max(min_clearance + 20, 60));
 max_climb = get_param(params, 'max_climb', 30);
 max_descent = get_param(params, 'max_descent', 30);
 edge_sample_spacing = get_param(params, 'edge_sample_spacing', 10);

 shadow_bias = get_param(params, 'shadow_bias', 0.7);
 animate = get_param(params, 'animate', false);
 plot_interval = max(1, round(get_param(params, 'plot_interval', 200)));
 informed_sampling = get_param(params, 'informed_sampling', true);
 informed_max_attempts = max(10, round(get_param(params, 'informed_max_attempts', 100)));

 radar_hard_constraint = get_param(params, 'radar_hard_constraint', false);
 radar_visibility_threshold = get_param(params, 'radar_visibility_threshold', 0.5);

 kd_rebuild_interval = max(1, round(get_param(params, 'kd_rebuild_interval', 64)));

 use_parallel_rewire = get_param(params, 'use_parallel_rewire', false);
 parallel_neighbor_threshold = max(2, round(get_param(params, 'parallel_neighbor_threshold', 96)));
 if use_parallel_rewire && ~can_use_parfor()
     warning('RRT*: parfor unavailable. Falling back to serial rewiring loops.');
     use_parallel_rewire = false;
 end

 enable_shortcut_prune = get_param(params, 'enable_shortcut_prune', false);
 shortcut_trials = max(0, round(get_param(params, 'shortcut_trials', 64)));
 shortcut_min_index_gap = max(2, round(get_param(params, 'shortcut_min_index_gap', 2)));

 % Height source integration:
 % - default: terrain_map.get_height
 % - optional: params.terrain_mesh.get_height
 using_mesh_height = isfield(params, 'terrain_mesh') && ~isempty(params.terrain_mesh);
 if using_mesh_height
     height_source = params.terrain_mesh;
 else
     height_source = terrain;
 end

 mesh_const = [];
 if use_parallel_rewire && using_mesh_height
     try
         mesh_const = parallel.pool.Constant(height_source);
     catch ME
         warning('RRT*: Could not wrap terrain_mesh in parallel.pool.Constant (%s). Serial rewiring will be used.', ME.message);
         use_parallel_rewire = false;
     end
 end

 threat_available = ~isempty(threat) && isprop(threat, 'computed') && threat.computed;
 threat_const = [];
 use_threat_const = false;
 if use_parallel_rewire && threat_available
     try
         threat_const = parallel.pool.Constant(threat);
         use_threat_const = true;
     catch ME
         warning('RRT*: Could not wrap threat_map in parallel.pool.Constant (%s). Serial rewiring will be used.', ME.message);
         use_parallel_rewire = false;
     end
 end

 cleanup_mesh_const = []; %#ok<NASGU>
 cleanup_threat_const = []; %#ok<NASGU>
 if ~isempty(mesh_const)
     cleanup_mesh_const = onCleanup(@() delete(mesh_const));
 end
 if ~isempty(threat_const)
     cleanup_threat_const = onCleanup(@() delete(threat_const));
 end

 % Bounds
 if isfield(params, 'bounds')
     bounds = params.bounds;
 else
     tb = terrain.bounds;
     if isprop(terrain, 'Z')
         terrain_min = min(terrain.Z(:));
         terrain_max = max(terrain.Z(:));
     else
         terrain_min = 0;
         terrain_max = 500;
     end
     max_flight_alt = get_param(params, 'max_flight_alt', 500);
     min_alt = terrain_min + min_clearance;
     max_alt = terrain_max + max_flight_alt;
     bounds = [tb(1), tb(2), tb(3), tb(4), -max_alt, -min_alt];
 end

 start = start(:);
 goal = goal(:);

 % Optional live plotting handles
 anim_axes = [];
 if animate
     if isfield(params, 'anim_axes') && ~isempty(params.anim_axes) && isgraphics(params.anim_axes, 'axes')
         anim_axes = params.anim_axes;
     elseif isfield(params, 'anim_fig') && ~isempty(params.anim_fig) && isgraphics(params.anim_fig, 'figure')
         anim_axes = get(params.anim_fig, 'CurrentAxes');
     end

     if isempty(anim_axes) || ~isgraphics(anim_axes, 'axes')
         figure('Name', 'RRT* Live', 'Position', [120, 120, 1000, 760]);
         anim_axes = axes('Parent', gcf);
         hold(anim_axes, 'on');
         grid(anim_axes, 'on');
         xlabel(anim_axes, 'North [m]');
         ylabel(anim_axes, 'East [m]');
         zlabel(anim_axes, 'Altitude [m]');
         view(anim_axes, 35, 30);
     end
     hold(anim_axes, 'on');
     h_tree = plot3(anim_axes, nan, nan, nan, '-', 'Color', [0.15 0.55 1.0], ...
         'LineWidth', 0.8, 'DisplayName', 'RRT* Tree');
     h_best = plot3(anim_axes, nan, nan, nan, 'm-', 'LineWidth', 2.4, 'DisplayName', 'Best Path');
 else
     h_tree = [];
     h_best = [];
 end

 % Tree preallocation
 max_nodes = max_iter + 2;
 nodes = nan(3, max_nodes);
 parents = zeros(1, max_nodes);
 costs = inf(1, max_nodes);

 % Child adjacency list for O(descendants) propagation
 children = repmat({zeros(1, 0)}, 1, max_nodes);

 nodes(:, 1) = start;
 parents(1) = 0;
 costs(1) = 0;
 node_count = 1;

 goal_reached = false;
 goal_node_idx = 0;
 goal_tolerance = base_step_size;

 c_min = norm(goal - start);
 x_center = (start + goal) / 2;
 C = rotation_matrix_3d(start, goal);

 % Radar cache for vectorized dynamic RCS modifier
 radar_cache = struct('count', 0, 'positions_ned', zeros(3, 0));
 if threat_available && isprop(threat, 'radars') && ~isempty(threat.radars)
     n_rad = length(threat.radars);
     radar_positions = zeros(3, n_rad);
     for r = 1:n_rad
         rp = threat.radars{r}.position(:);
         radar_positions(:, r) = [rp(1); rp(2); -rp(3)];
     end
     radar_cache.count = n_rad;
     radar_cache.positions_ned = radar_positions;
 end

 % KD-tree state
 kd = struct();
 kd.enabled = can_use_kdtree();
 kd.tree = [];
 kd.count = 0;
 kd.pending = inf;  % force first build
 kd.rebuild_interval = kd_rebuild_interval;
 if ~kd.enabled
     warning('RRT*: KD-tree functions unavailable. Falling back to vectorized full scans.');
 end

 iter = 0;
 for iter = 1:max_iter
     goal_updated_this_iter = false;

     c_best = inf;
     if goal_reached
         c_best = costs(goal_node_idx);
     end

     goal_bias_curr = goal_bias;
     if goal_reached
         goal_bias_curr = goal_bias_after_goal;
     end

     if rand < goal_bias_curr
         q_rand = goal;
     else
         q_rand = sample_informed(bounds, height_source, threat, min_clearance, ...
             start, goal, c_best, c_min, x_center, C, shadow_bias, ...
             radar_hard_constraint, radar_visibility_threshold, using_mesh_height, ...
             informed_sampling, informed_max_attempts);
     end

     % O(N log N): nearest query via KD-tree
     kd = rebuild_kd_if_needed(nodes, node_count, kd, false);
     [q_near_idx, ~] = query_nearest_kd(nodes, node_count, q_rand, kd);
     q_near = nodes(:, q_near_idx);

     q_new = steer_adaptive(q_near, q_rand, base_step_size, max_step_size, ...
         max_climb, max_descent, height_source, using_mesh_height);

     if ~is_valid_point(q_new, height_source, min_clearance, bounds, using_mesh_height)
         continue;
     end

     [edge_valid, edge_cost] = evaluate_edge_worker(q_near, q_new, ...
         height_source, mesh_const, using_mesh_height, ...
         threat, threat_const, use_threat_const, threat_available, ...
         edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
         min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
     if ~edge_valid || isinf(edge_cost)
         continue;
     end

     new_cost = costs(q_near_idx) + edge_cost;

     % Rewiring neighborhood selection (radius or k-nearest)
     rewire_radius = NaN;
     rewire_k = NaN;
     switch rewire_mode
         case 'fixed'
             rewire_radius = rewire_radius_max;
             near_indices = query_near_kd(nodes, node_count, q_new, rewire_radius, kd);
         case 'knearest'
             rewire_k = compute_dynamic_rewire_k(node_count, spatial_dim, ...
                 rewire_k_const, rewire_k_min, rewire_k_max);
             near_indices = query_knear_kd(nodes, node_count, q_new, rewire_k, kd);
         otherwise  % 'dynamic'
             rewire_radius = compute_dynamic_rewire_radius(node_count, spatial_dim, ...
                 rewire_radius_max, rewire_radius_min, rewire_gamma_dyn);
             near_indices = query_near_kd(nodes, node_count, q_new, rewire_radius, kd);
     end

     % Keep nearest parent candidate available in all modes.
     near_indices = unique([near_indices, q_near_idx], 'stable');
     if isempty(near_indices)
         near_indices = q_near_idx;
     end

     best_parent = q_near_idx;
     best_cost = new_cost;

     cand_indices = near_indices(near_indices ~= q_near_idx);
     if ~isempty(cand_indices)
         cand_nodes = nodes(:, cand_indices);
         cand_base_costs = costs(cand_indices);
         cand_total_costs = inf(1, numel(cand_indices));

         if use_parallel_rewire && numel(cand_indices) >= parallel_neighbor_threshold
             parfor i_c = 1:numel(cand_indices)
                 [ok_c, c_edge_cost] = evaluate_edge_worker(cand_nodes(:, i_c), q_new, ...
                     height_source, mesh_const, using_mesh_height, ...
                     threat, threat_const, use_threat_const, threat_available, ...
                     edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
                     min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
                 if ok_c
                     cand_total_costs(i_c) = cand_base_costs(i_c) + c_edge_cost;
                 end
             end
         else
             for i_c = 1:numel(cand_indices)
                 [ok_c, c_edge_cost] = evaluate_edge_worker(cand_nodes(:, i_c), q_new, ...
                     height_source, mesh_const, using_mesh_height, ...
                     threat, threat_const, use_threat_const, threat_available, ...
                     edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
                     min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
                 if ok_c
                     cand_total_costs(i_c) = cand_base_costs(i_c) + c_edge_cost;
                 end
             end
         end

         [cand_best_cost, cand_best_loc] = min(cand_total_costs);
         if cand_best_cost < best_cost
             best_cost = cand_best_cost;
             best_parent = cand_indices(cand_best_loc);
         end
     end

    if node_count + 1 > max_nodes
        warning('RRT*: Reached max preallocated nodes (%d). Stopping early.', max_nodes);
        break;
    end

    node_count = node_count + 1;
    new_idx = node_count;
    nodes(:, new_idx) = q_new;
    parents(new_idx) = best_parent;
    costs(new_idx) = best_cost;
    children{new_idx} = zeros(1, 0);
    children{best_parent}(end+1) = new_idx;
    kd.pending = kd.pending + 1;

    rew_indices = near_indices(near_indices ~= best_parent & near_indices ~= new_idx);
    if ~isempty(rew_indices)
        rew_nodes = nodes(:, rew_indices);
        rew_new_costs = inf(1, numel(rew_indices));
        base_new_cost = costs(new_idx);

        if use_parallel_rewire && numel(rew_indices) >= parallel_neighbor_threshold
            parfor i_r = 1:numel(rew_indices)
                [ok_r, r_edge_cost] = evaluate_edge_worker(q_new, rew_nodes(:, i_r), ...
                    height_source, mesh_const, using_mesh_height, ...
                    threat, threat_const, use_threat_const, threat_available, ...
                    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
                    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
                if ok_r
                    rew_new_costs(i_r) = base_new_cost + r_edge_cost;
                end
            end
        else
            for i_r = 1:numel(rew_indices)
                [ok_r, r_edge_cost] = evaluate_edge_worker(q_new, rew_nodes(:, i_r), ...
                    height_source, mesh_const, using_mesh_height, ...
                    threat, threat_const, use_threat_const, threat_available, ...
                    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
                    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
                if ok_r
                    rew_new_costs(i_r) = base_new_cost + r_edge_cost;
                end
            end
        end

        for i_r = 1:numel(rew_indices)
            idx = rew_indices(i_r);
            candidate_cost = rew_new_costs(i_r);

            if candidate_cost + 1e-9 < costs(idx)
                old_parent = parents(idx);
                old_cost = costs(idx);

                parents(idx) = new_idx;
                children = remove_child(children, old_parent, idx);
                children{new_idx}(end+1) = idx;

                delta = candidate_cost - old_cost;
                costs(idx) = candidate_cost;
                costs = propagate_descendant_delta(children, costs, idx, delta);
            end
        end
    end

    dist_to_goal = norm(q_new - goal);
    if dist_to_goal < goal_tolerance
        [goal_ok, goal_edge_cost] = evaluate_edge_worker(q_new, goal, ...
            height_source, mesh_const, using_mesh_height, ...
            threat, threat_const, use_threat_const, threat_available, ...
            edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
            min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);

        if goal_ok && ~isinf(goal_edge_cost)
            total_goal_cost = costs(new_idx) + goal_edge_cost;

            if ~goal_reached
                if node_count + 1 > max_nodes
                    warning('RRT*: Could not add goal node (capacity reached).');
                else
                    node_count = node_count + 1;
                    goal_node_idx = node_count;
                    nodes(:, goal_node_idx) = goal;
                    parents(goal_node_idx) = new_idx;
                    costs(goal_node_idx) = total_goal_cost;
                    children{goal_node_idx} = zeros(1, 0);
                    children{new_idx}(end+1) = goal_node_idx;
                    kd.pending = kd.pending + 1;
                    goal_reached = true;
                    fprintf('Goal reached at iter %d! Cost: %.2f\n', iter, total_goal_cost);
                    goal_updated_this_iter = true;
                end
            elseif total_goal_cost < costs(goal_node_idx)
                old_parent = parents(goal_node_idx);
                old_cost = costs(goal_node_idx);
                parents(goal_node_idx) = new_idx;
                children = remove_child(children, old_parent, goal_node_idx);
                children{new_idx}(end+1) = goal_node_idx;
                delta = total_goal_cost - old_cost;
                costs(goal_node_idx) = total_goal_cost;
                costs = propagate_descendant_delta(children, costs, goal_node_idx, delta);
                fprintf('Goal improved at iter %d! Cost: %.2f\n', iter, total_goal_cost);
                goal_updated_this_iter = true;
            end
        end
    end

    if mod(iter, 500) == 0
        best_disp = NaN;
        if goal_reached
            best_disp = costs(goal_node_idx);
        end
        if strcmp(rewire_mode, 'knearest')
            fprintf('RRT* iter %d, tree: %d, best cost: %.2f, k_rewire: %d\n', ...
                iter, node_count, best_disp, rewire_k);
        else
            fprintf('RRT* iter %d, tree: %d, best cost: %.2f, r_rewire: %.2f\n', ...
                iter, node_count, best_disp, rewire_radius);
        end
    end

    if animate && isgraphics(anim_axes, 'axes') && (mod(iter, plot_interval) == 0 || goal_updated_this_iter)
        [x_tree, y_tree, z_tree] = build_tree_lines(nodes, parents, node_count);
        set(h_tree, 'XData', x_tree, 'YData', y_tree, 'ZData', z_tree);
        if goal_reached && goal_node_idx > 0
            path_live = extract_path(nodes, parents, goal_node_idx);
            set(h_best, 'XData', path_live(1, :), 'YData', path_live(2, :), 'ZData', -path_live(3, :));
        end
        best_disp = NaN;
        if goal_reached
            best_disp = costs(goal_node_idx);
        end
        title(anim_axes, sprintf('RRT* Live | iter=%d | tree=%d | best cost=%.2f', iter, node_count, best_disp));
        drawnow limitrate nocallbacks;
    end
 end

 if animate && isgraphics(anim_axes, 'axes')
     [x_tree, y_tree, z_tree] = build_tree_lines(nodes, parents, node_count);
     set(h_tree, 'XData', x_tree, 'YData', y_tree, 'ZData', z_tree);
     drawnow;
 end

 if goal_reached
     path = extract_path(nodes, parents, goal_node_idx);
     final_cost = costs(goal_node_idx);
 else
     [closest_idx, ~] = query_nearest_kd(nodes, node_count, goal, kd);
     path = extract_path(nodes, parents, closest_idx);
     final_cost = costs(closest_idx);
     warning('Goal not reached.');
 end

 shortcut_info = struct('enabled', enable_shortcut_prune, ...
     'attempts', 0, 'accepted', 0, 'node_reduction', 0);
 if enable_shortcut_prune && size(path, 2) >= 3
     [path, shortcut_info] = shortcut_prune_path(path, shortcut_trials, shortcut_min_index_gap, ...
         height_source, mesh_const, using_mesh_height, ...
         threat, threat_const, use_threat_const, threat_available, ...
         edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
         min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);

     [path_ok, smoothed_cost] = evaluate_path_cost(path, ...
         height_source, mesh_const, using_mesh_height, ...
         threat, threat_const, use_threat_const, threat_available, ...
         edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
         min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
     if path_ok
         final_cost = smoothed_cost;
     end
 end

 planning_time = toc;
 path_length = compute_path_length(path);
 path_risk = compute_path_risk(path, threat);

 info.success = goal_reached;
 info.iterations = iter;
 info.tree_size = node_count;
 info.path_cost = final_cost;
 info.path_length = path_length;
 info.path_risk = path_risk;
 info.planning_time = planning_time;
 info.shortcut = shortcut_info;

 fprintf('RRT* Finished. Time: %.2fs. Length: %.1fm. Risk: %.2f.\n', planning_time, path_length, path_risk);

end

function q = sample_informed(bounds, height_source, threat, min_clearance, ...
    start, goal, c_best, c_min, x_center, C, shadow_bias, ...
    radar_hard_constraint, radar_visibility_threshold, using_mesh_height, ...
    informed_sampling, max_attempts)
% Informed sampling with shadow-bias / hard-visibility rejection.

for k = 1:max_attempts
    if ~informed_sampling || ~isfinite(c_best) || (c_best <= c_min + 1e-9)
        N = bounds(1) + rand * (bounds(2) - bounds(1));
        E = bounds(3) + rand * (bounds(4) - bounds(3));
        D = bounds(5) + rand * (bounds(6) - bounds(5));
        cand = [N; E; D];
    else
        major = c_best / 2;
        minor_term = max(c_best^2 - c_min^2, 0);
        minor = sqrt(minor_term) / 2;
        L = diag([major, minor, minor]);
        x_ball = sample_unit_ball();
        cand = C * L * x_ball + x_center;
    end

    if cand(1) < bounds(1) || cand(1) > bounds(2) || ...
       cand(2) < bounds(3) || cand(2) > bounds(4) || ...
       cand(3) < bounds(5) || cand(3) > bounds(6)
        continue;
    end

    alt = -cand(3);
    terrain_h = query_height_source(cand(1), cand(2), height_source, [], using_mesh_height);
    terrain_h = terrain_h(1);
    if isnan(terrain_h) || alt <= terrain_h + min_clearance
        continue;
    end

    if ~isempty(threat) && isprop(threat, 'computed') && threat.computed
        point_risk = 0;
        try
            risk_val = threat.get_risk(cand(1), cand(2), alt);
            point_risk = risk_val(1);
        catch
            point_risk = 0;
        end
        visible_to_any = point_risk >= radar_visibility_threshold;

        if radar_hard_constraint && visible_to_any
            continue;
        end
        if ~radar_hard_constraint && rand < shadow_bias && visible_to_any
            continue;
        end
    end

    q = cand;
    return;
end

N = bounds(1) + rand * (bounds(2) - bounds(1));
E = bounds(3) + rand * (bounds(4) - bounds(3));
h = query_height_source(N, E, height_source, [], using_mesh_height);
h = h(1);
if isnan(h)
    h = 0;
end
q = [N; E; -(h + min_clearance + 20)];
end

function [x, y, z] = build_tree_lines(nodes, parents, node_count)
if node_count < 2
    x = nan;
    y = nan;
    z = nan;
    return;
end

child_idx = 2:node_count;
par_idx = parents(child_idx);
valid = par_idx > 0 & par_idx <= node_count;
child_idx = child_idx(valid);
par_idx = par_idx(valid);

if isempty(child_idx)
    x = nan;
    y = nan;
    z = nan;
    return;
end

n_edges = numel(child_idx);
x = nan(1, 3 * n_edges);
y = nan(1, 3 * n_edges);
z = nan(1, 3 * n_edges);

x(1:3:end) = nodes(1, par_idx);
x(2:3:end) = nodes(1, child_idx);
y(1:3:end) = nodes(2, par_idx);
y(2:3:end) = nodes(2, child_idx);
z(1:3:end) = -nodes(3, par_idx);
z(2:3:end) = -nodes(3, child_idx);
end

function x = sample_unit_ball()
while true
    x = 2 * rand(3, 1) - 1;
    if norm(x) <= 1
        return;
    end
end
end

function R = rotation_matrix_3d(a, b)
u = (b - a) / norm(b - a);
id1 = [1; 0; 0];
v = cross(id1, u);
s = norm(v);
c = dot(id1, u);

if s < 1e-6
    R = eye(3);
else
    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + vx + vx^2 * ((1 - c) / s^2);
end
end

function q_new = steer_adaptive(q_from, q_to, base_step, max_step, ...
    max_climb, max_descent, height_source, using_mesh_height)
% Adaptive steering based on local terrain slope.

direction = q_to - q_from;
dist = norm(direction);
if dist < 1e-9
    q_new = q_to;
    return;
end
direction = direction / dist;

probe_dist = min(10, base_step / 5);
h0 = query_height_source(q_from(1), q_from(2), height_source, [], using_mesh_height);
h_n = query_height_source(q_from(1) + probe_dist, q_from(2), height_source, [], using_mesh_height);
h_e = query_height_source(q_from(1), q_from(2) + probe_dist, height_source, [], using_mesh_height);
h0 = h0(1);
h_n = h_n(1);
h_e = h_e(1);

if any(isnan([h0, h_n, h_e]))
    slope_angle = 0;
else
    slope_n = (h_n - h0) / probe_dist;
    slope_e = (h_e - h0) / probe_dist;
    slope_angle = atan(hypot(slope_n, slope_e));
end

blend = 1 / (1 + exp(-10 * (slope_angle - deg2rad(20))));
step_size = max_step * (1 - blend) + base_step * blend;
step = min(dist, step_size) * direction;

horizontal_dist = norm(step(1:2));
vertical_dist = step(3);
if horizontal_dist > 1e-9
    max_climb_dist = horizontal_dist * tand(max_climb);
    max_descent_dist = horizontal_dist * tand(max_descent);
    if -vertical_dist > max_climb_dist
        step(3) = -max_climb_dist;
    elseif vertical_dist > max_descent_dist
        step(3) = max_descent_dist;
    end
end

q_new = q_from + step;
end

function valid = is_valid_point(q, height_source, min_clearance, bounds, using_mesh_height)
if q(1) < bounds(1) || q(1) > bounds(2) || ...
   q(2) < bounds(3) || q(2) > bounds(4) || ...
   q(3) < bounds(5) || q(3) > bounds(6)
    valid = false;
    return;
end

alt = -q(3);
terrain_h = query_height_source(q(1), q(2), height_source, [], using_mesh_height);
terrain_h = terrain_h(1);
if isnan(terrain_h) || alt < terrain_h + min_clearance
    valid = false;
    return;
end
valid = true;
end

function [valid, cost] = evaluate_edge_worker(q1, q2, ...
    height_source, mesh_const, using_mesh_height, ...
    threat, threat_const, use_threat_const, threat_available, ...
    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache)
% Vectorized edge validation and edge-cost evaluation.

[samples, dist] = sample_edge_points(q1, q2, edge_sample_spacing, 3);
if dist < 1e-9
    valid = true;
    cost = 0;
    return;
end

N = samples(1, :);
E = samples(2, :);
alt = -samples(3, :);

terrain_h = query_height_source(N, E, height_source, mesh_const, using_mesh_height);
if any(isnan(terrain_h))
    valid = false;
    cost = inf;
    return;
end

if any(alt < (terrain_h + min_clearance))
    valid = false;
    cost = inf;
    return;
end

cost_dist = alpha * dist;
cost_radar = 0;

if threat_available
    if use_threat_const
        threat_obj = threat_const.Value;
    else
        threat_obj = threat;
    end

    try
        P_base = threat_obj.get_risk(N, E, alt);
        P_base = reshape(P_base, 1, []);
    catch
        P_base = zeros(1, numel(N));
    end

    if radar_hard_constraint && any(P_base >= radar_visibility_threshold)
        valid = false;
        cost = inf;
        return;
    end

    if ~radar_hard_constraint && beta > 0
        rcs_mod = compute_dynamic_rcs_mod(samples, q2 - q1, radar_cache.positions_ned);
        P_effective = min(P_base .* rcs_mod, 0.99);
        radar_penalty = -log1p(-P_effective);
        cost_radar = beta * sum(radar_penalty) * (dist / max(1, numel(P_effective)));
    end
end

cost_alt = 0;
if gamma_pref > 0
    agl = alt - terrain_h;
    agl_scale = max(preferred_agl, 1);
    err_norm = (agl - preferred_agl) / agl_scale;
    cost_alt = gamma_pref * dist * mean(err_norm.^2);
end

cost = cost_dist + cost_radar + cost_alt;
valid = true;
end

function [samples, dist] = sample_edge_points(q1, q2, spacing, min_samples)
dist = norm(q2 - q1);
if dist < 1e-9
    samples = q1;
    return;
end

n_samples = max(min_samples, ceil(dist / max(spacing, 1e-6)) + 1);
t = linspace(0, 1, n_samples);
samples = q1 + (q2 - q1) * t;
end

function rcs_mod = compute_dynamic_rcs_mod(samples, path_vec, radar_positions_ned)
% Vectorized nearest-radar aspect-based RCS scaling.

n_samples = size(samples, 2);
n_radars = size(radar_positions_ned, 2);
if n_samples == 0 || n_radars == 0
    rcs_mod = ones(1, n_samples);
    return;
end

path_norm = norm(path_vec);
if path_norm < 1e-9
    path_unit = [1; 0; 0];
else
    path_unit = path_vec / path_norm;
end

sample_exp = permute(samples, [3, 2, 1]);
radar_exp = permute(radar_positions_ned.', [1, 3, 2]);
deltas = radar_exp - sample_exp;
d2 = sum(deltas.^2, 3);
[~, nearest_idx] = min(d2, [], 1);

rcs_mod = ones(1, n_samples);
for r = 1:n_radars
    idx_mask = find(nearest_idx == r);
    if isempty(idx_mask)
        continue;
    end

    v_los = radar_positions_ned(:, r) - samples(:, idx_mask);
    los_norm = sqrt(sum(v_los.^2, 1));
    valid_local = find(los_norm > 1e-9);
    if isempty(valid_local)
        continue;
    end

    v_los_valid = v_los(:, valid_local) ./ los_norm(valid_local);
    v_path = repmat(path_unit, 1, size(v_los_valid, 2));
    sin_theta = vecnorm(cross(v_path, v_los_valid, 1), 2, 1);
    rcs_mod(idx_mask(valid_local)) = 0.1 + 0.9 * (sin_theta.^2);
end
end

function h = query_height_source(N, E, height_source, mesh_const, using_mesh_height)
if nargin < 4
    mesh_const = [];
end
if nargin < 5
    using_mesh_height = false;
end

if using_mesh_height && ~isempty(mesh_const)
    h = mesh_const.Value.get_height(N, E);
else
    h = height_source.get_height(N, E);
end
h = reshape(h, 1, []);
end

function r = compute_dynamic_rewire_radius(n_nodes, dim, r_max, r_min, gamma_dyn)
% r_n = min(r_max, max(r_min, gamma*(log(n)/n)^(1/d)))

n_eff = max(2, n_nodes);
r = gamma_dyn * (log(n_eff) / n_eff)^(1 / dim);
r = min(r_max, max(r_min, r));
end

function k = compute_dynamic_rewire_k(n_nodes, dim, k_const, k_min, k_max)
% k_n = clamp(ceil(k_const * log(n)), [k_min, k_max])
% dim is kept for interface symmetry and future dimensional tuning.

n_eff = max(2, n_nodes);
k = ceil(k_const * log(n_eff));
k = max(k_min, k);
k = min([k, k_max, n_eff]);

% Ensure at least enough neighbors for geometric branching.
k = max(k, dim + 1);
k = min(k, n_eff);
end

function kd = rebuild_kd_if_needed(nodes, node_count, kd, force_rebuild)
% Periodic KD-tree rebuild to amortize insertion cost.
if ~kd.enabled
    return;
end

if force_rebuild || isempty(kd.tree) || kd.pending >= kd.rebuild_interval
    data = nodes(:, 1:node_count).';
    kd.tree = createns(data, 'NSMethod', 'kdtree');
    kd.count = node_count;
    kd.pending = 0;
end
end

function [idx, dist] = query_nearest_kd(nodes, node_count, q, kd)
% Nearest query from KD-tree + exact check of unreindexed tail.
if ~kd.enabled
    diffs = nodes(:, 1:node_count) - q;
    d2 = sum(diffs.^2, 1);
    [d2_min, idx] = min(d2);
    dist = sqrt(d2_min);
    return;
end

idx = 1;
dist = inf;
if ~isempty(kd.tree) && kd.count > 0
    [idx_kd, dist_kd] = knnsearch(kd.tree, q.');
    idx = idx_kd;
    dist = dist_kd;
end

if kd.count < node_count
    tail_idx = (kd.count + 1):node_count;
    tail_diffs = nodes(:, tail_idx) - q;
    tail_d2 = sum(tail_diffs.^2, 1);
    [tail_min_d2, tail_loc] = min(tail_d2);
    tail_dist = sqrt(tail_min_d2);
    if tail_dist < dist
        idx = tail_idx(tail_loc);
        dist = tail_dist;
    end
end
end

function indices = query_near_kd(nodes, node_count, q, radius, kd)
% Radius query from KD-tree + exact check of unreindexed tail.
if ~kd.enabled
    diffs = nodes(:, 1:node_count) - q;
    d = vecnorm(diffs, 2, 1);
    indices = find(d <= radius);
    return;
end

indices = [];
if ~isempty(kd.tree) && kd.count > 0
    in_kd = rangesearch(kd.tree, q.', radius);
    indices = in_kd{1};
end

if kd.count < node_count
    tail_idx = (kd.count + 1):node_count;
    tail_diffs = nodes(:, tail_idx) - q;
    tail_d = vecnorm(tail_diffs, 2, 1);
    tail_in = tail_idx(tail_d <= radius);
    indices = [indices, tail_in]; %#ok<AGROW>
end

if isempty(indices)
    return;
end
indices = unique(indices, 'stable');
end

function indices = query_knear_kd(nodes, node_count, q, k_neighbors, kd)
% K-nearest query from KD-tree + exact merge with unreindexed tail.

k_neighbors = max(1, min(round(k_neighbors), node_count));

if ~kd.enabled
    diffs = nodes(:, 1:node_count) - q;
    d2 = sum(diffs.^2, 1);
    [~, ord] = sort(d2, 'ascend');
    indices = ord(1:k_neighbors);
    return;
end

idx_all = [];
dist_all = [];

if ~isempty(kd.tree) && kd.count > 0
    k_kd = min(k_neighbors, kd.count);
    [idx_kd, dist_kd] = knnsearch(kd.tree, q.', 'K', k_kd);
    idx_kd = reshape(idx_kd, 1, []);
    dist_kd = reshape(dist_kd, 1, []);
    idx_all = [idx_all, idx_kd]; %#ok<AGROW>
    dist_all = [dist_all, dist_kd]; %#ok<AGROW>
end

if kd.count < node_count
    tail_idx = (kd.count + 1):node_count;
    tail_d = vecnorm(nodes(:, tail_idx) - q, 2, 1);
    idx_all = [idx_all, tail_idx]; %#ok<AGROW>
    dist_all = [dist_all, tail_d]; %#ok<AGROW>
end

if isempty(idx_all)
    indices = [];
    return;
end

[~, ord] = sort(dist_all, 'ascend');
idx_sorted = idx_all(ord);
idx_unique = unique(idx_sorted, 'stable');
indices = idx_unique(1:min(k_neighbors, numel(idx_unique)));
end

function children = remove_child(children, parent_idx, child_idx)
if parent_idx <= 0 || parent_idx > numel(children)
    return;
end
if isempty(children{parent_idx})
    return;
end
children{parent_idx} = children{parent_idx}(children{parent_idx} ~= child_idx);
end

function costs = propagate_descendant_delta(children, costs, root_idx, delta)
% Update all descendants with constant subtree delta in O(descendants).
if abs(delta) < 1e-12
    return;
end

queue = children{root_idx};
head = 1;
while head <= numel(queue)
    cidx = queue(head);
    head = head + 1;
    costs(cidx) = costs(cidx) + delta;
    kids = children{cidx};
    if ~isempty(kids)
        queue = [queue, kids]; %#ok<AGROW>
    end
end
end

function path = extract_path(nodes, parents, goal_idx)
if goal_idx <= 0
    path = nodes(:, 1);
    return;
end

idx_chain = zeros(1, numel(parents));
k = 0;
idx = goal_idx;
while idx > 0
    k = k + 1;
    idx_chain(k) = idx;
    idx = parents(idx);
end
idx_chain = idx_chain(k:-1:1);
path = nodes(:, idx_chain);
end

function [path_out, info] = shortcut_prune_path(path_in, n_trials, min_gap, ...
    height_source, mesh_const, using_mesh_height, ...
    threat, threat_const, use_threat_const, threat_available, ...
    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache)
% Optional post-process shortcutting after tree convergence.

path_out = path_in;
info = struct('enabled', true, 'attempts', 0, 'accepted', 0, 'node_reduction', 0);

if n_trials <= 0 || size(path_out, 2) < (min_gap + 1)
    return;
end

for trial = 1:n_trials
    n_pts = size(path_out, 2);
    if n_pts < (min_gap + 1)
        break;
    end

    i = randi([1, n_pts - min_gap]);
    j = randi([i + min_gap, n_pts]);
    info.attempts = info.attempts + 1;

    [ok_short, short_cost] = evaluate_edge_worker(path_out(:, i), path_out(:, j), ...
        height_source, mesh_const, using_mesh_height, ...
        threat, threat_const, use_threat_const, threat_available, ...
        edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
        min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
    if ~ok_short || isinf(short_cost)
        continue;
    end

    [ok_old, old_cost] = evaluate_path_window_cost(path_out, i, j, ...
        height_source, mesh_const, using_mesh_height, ...
        threat, threat_const, use_threat_const, threat_available, ...
        edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
        min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
    if ~ok_old
        continue;
    end

    if short_cost + 1e-9 < old_cost
        path_out = [path_out(:, 1:i), path_out(:, j:end)];
        info.accepted = info.accepted + 1;
    end
end

info.node_reduction = size(path_in, 2) - size(path_out, 2);
end

function [ok, total_cost] = evaluate_path_window_cost(path, idx_start, idx_end, ...
    height_source, mesh_const, using_mesh_height, ...
    threat, threat_const, use_threat_const, threat_available, ...
    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache)
% Cost of path(idx_start:idx_end) with the planner's edge objective.

ok = true;
total_cost = 0;
for k = idx_start:(idx_end - 1)
    [edge_ok, edge_cost] = evaluate_edge_worker(path(:, k), path(:, k + 1), ...
        height_source, mesh_const, using_mesh_height, ...
        threat, threat_const, use_threat_const, threat_available, ...
        edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
        min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
    if ~edge_ok || isinf(edge_cost)
        ok = false;
        total_cost = inf;
        return;
    end
    total_cost = total_cost + edge_cost;
end
end

function [ok, total_cost] = evaluate_path_cost(path, ...
    height_source, mesh_const, using_mesh_height, ...
    threat, threat_const, use_threat_const, threat_available, ...
    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache)

if size(path, 2) < 2
    ok = true;
    total_cost = 0;
    return;
end

[ok, total_cost] = evaluate_path_window_cost(path, 1, size(path, 2), ...
    height_source, mesh_const, using_mesh_height, ...
    threat, threat_const, use_threat_const, threat_available, ...
    edge_sample_spacing, alpha, beta, gamma_pref, preferred_agl, ...
    min_clearance, radar_hard_constraint, radar_visibility_threshold, radar_cache);
end

function len = compute_path_length(path)
if size(path, 2) < 2
    len = 0;
    return;
end
d = diff(path, 1, 2);
len = sum(vecnorm(d, 2, 1));
end

function risk = compute_path_risk(path, threat)
risk = 0;
if isempty(threat) || ~isprop(threat, 'computed') || ~threat.computed
    return;
end
try
    p = threat.get_risk(path(1, :), path(2, :), -path(3, :));
    risk = sum(p);
catch
    risk = 0;
end
end

function val = get_param(params, name, default)
if isfield(params, name)
    val = params.(name);
else
    val = default;
end
end

function result = ternary(condition, true_val, false_val)
if condition
    result = true_val;
else
    result = false_val;
end
end

function tf = can_use_parfor()
tf = license('test', 'Distrib_Computing_Toolbox') && ...
     (exist('parfor', 'builtin') > 0 || exist('parfor', 'file') > 0);
end

function tf = can_use_kdtree()
tf = (exist('createns', 'builtin') > 0 || exist('createns', 'file') > 0) && ...
     (exist('knnsearch', 'builtin') > 0 || exist('knnsearch', 'file') > 0) && ...
     (exist('rangesearch', 'builtin') > 0 || exist('rangesearch', 'file') > 0);
end
