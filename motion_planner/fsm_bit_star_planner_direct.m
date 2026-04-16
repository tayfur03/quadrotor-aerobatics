function [path, cost, info] = fsm_bit_star_planner_direct(x_start, x_goal, env, params)
%FSM_BIT_STAR_PLANNER_DIRECT Direct BIT* variant using FSM heuristics.
%
% Inputs
%   x_start   [3x1] start point [N;E;Alt]
%   x_goal    [3x1] goal point  [N;E;Alt]
%   env       struct with direct, corridor-free planning data:
%     .terrain or .terrain_map      terrain_map object or terrain struct
%     .risk_grid or .radar_los_map  3D threat / LOS volume
%     .N_vec, .E_vec, .alt_vec      world axes for the 3D volume
%     .visibility_threshold         threshold for binary visible/hidden classification
%     .bounds_world (optional)      [Nmin Nmax Emin Emax Altmin Altmax]
%   params    optional BIT*/heuristic settings
%
% Outputs
%   path  [3xM] planned path [N;E;Alt] (empty if failure)
%   cost  scalar path cost
%   info  timing and convergence metrics

if nargin < 4
    params = struct();
end

validateattributes(x_start, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_start', 1);
validateattributes(x_goal, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_goal', 2);
x_start = x_start(:);
x_goal = x_goal(:);

[terrain_input, terrain_query, radar_los_map, risk_query, vis_threshold, bounds, point_clearance, env_meta] = ...
    parse_direct_env(env, params, x_start, x_goal);

meta = struct('fallback', true, 'fsm_time_s', 0, 'warmstart_cost', Inf, ...
    'warmstart_path', [], 'grid_x', [], 'grid_y', [], 'grid_z', [], ...
    'dx', [], 'dy', [], 'dz', []);
T_fwd = [];
T_bwd = [];
try
    fsm_params = params;
    if ~isfield(fsm_params, 'terrain_margin')
        fsm_params.terrain_margin = point_clearance;
    end
    if ~isfield(fsm_params, 'fsm_bounds_world')
        fsm_params.fsm_bounds_world = bounds;
    end
    [T_fwd, T_bwd, ~, meta] = compute_fsm_heuristic( ...
        terrain_input, radar_los_map, world_to_ned(x_start), world_to_ned(x_goal), fsm_params);
catch ME
    warning('fsm_bit_star_planner_direct:FSMError', ...
        'FSM preprocessing failed (%s). Falling back to Euclidean heuristic.', ME.message);
    meta.fallback = true;
end

max_batches = round(get_opt(params, 'max_batches', 80));
batch_size = round(get_opt(params, 'batch_size', 240));
max_nodes = round(get_opt(params, 'max_nodes', 20000));
eps_global = get_opt(params, 'eps_global', 0.12);
stealth_weight = get_opt(params, 'stealth_weight', 2.0);
edge_check_samples = round(get_opt(params, 'edge_check_samples', 7));
post_solution_batches = round(get_opt(params, 'post_solution_batches', 12));
max_time = get_opt(params, 'max_time', inf);

max_batches = max(1, max_batches);
batch_size = max(16, batch_size);
max_nodes = max(batch_size + 2, max_nodes);
eps_global = min(max(eps_global, 0.0), 0.5);
edge_check_samples = max(3, edge_check_samples);
post_solution_batches = max(0, post_solution_batches);
strict_zero_risk_tol = max(0, get_opt(params, 'strict_zero_risk_tol', 1e-6));

default_goal_radius = max(5.0, 0.02 * norm(bounds([2 4 6]) - bounds([1 3 5])));
goal_radius = get_opt(params, 'goal_radius', default_goal_radius);
min_r = get_opt(params, 'min_connection_radius', max(5.0, 0.012 * norm(bounds([2 4 6]) - bounds([1 3 5]))));
max_r = get_opt(params, 'max_connection_radius', max(40, 0.08 * norm(bounds([2 4 6]) - bounds([1 3 5]))));
eta_radius = get_opt(params, 'eta_radius', 1.8);

if isinf(stealth_weight)
    risk_limit = strict_zero_risk_tol;
else
    risk_limit = vis_threshold;
end

nodes = nan(3, max_nodes);
parents = zeros(1, max_nodes, 'uint32');
g_cost = inf(1, max_nodes);

nodes(:, 1) = x_start;
parents(1) = uint32(0);
g_cost(1) = 0.0;
node_count = 1;

warmstart_bound_cost = Inf;
warmstart_world = [];
if ~meta.fallback && ~isempty(meta.warmstart_path)
    warmstart_world = [meta.warmstart_path(1, :); meta.warmstart_path(2, :); -meta.warmstart_path(3, :)];
    warmstart_bound_cost = evaluate_world_path_cost(warmstart_world, risk_query, risk_limit, stealth_weight, point_clearance, terrain_query);
end

best_goal_cost = warmstart_bound_cost;
goal_parent = uint32(0);
first_solution_batch = NaN;
first_solution_time = NaN;
if isfinite(warmstart_bound_cost)
    first_solution_batch = 0;
    first_solution_time = 0;
end
samples_generated = 0;
samples_accepted = 0;
volume_draws = 0;
volume_accept_fsm = 0;
volume_accept_eucl = 0;

t_start = tic;
batches_run = 0;
stop_reason = 'max_batches_reached';

for batch_idx = 1:max_batches
    if toc(t_start) > max_time
        stop_reason = 'max_time';
        break;
    end

    batches_run = batch_idx;
    S_all = sample_uniform_bounds(bounds, batch_size);
    samples_generated = samples_generated + batch_size;

    in_bounds = points_in_bounds(S_all, bounds);
    risk_all = risk_query(S_all);
    safe_mask = in_bounds & isfinite(risk_all) & (risk_all <= risk_limit);
    safe_mask = safe_mask & points_clear_of_terrain(S_all, terrain_query, point_clearance);

    if isfinite(best_goal_cost)
        if ~meta.fallback
            tf_all = trilinear_interp_batch(T_fwd, S_all, meta);
            tb_all = trilinear_interp_batch(T_bwd, S_all, meta);
            informed_mask = isfinite(tf_all) & isfinite(tb_all) & ((tf_all + tb_all) < best_goal_cost);
        else
            informed_mask = (vecnorm(S_all - x_start, 2, 1) + vecnorm(S_all - x_goal, 2, 1)) < best_goal_cost;
        end
        eucl_mask = (vecnorm(S_all - x_start, 2, 1) + vecnorm(S_all - x_goal, 2, 1)) < best_goal_cost;
    else
        informed_mask = true(1, batch_size);
        eucl_mask = true(1, batch_size);
    end

    volume_draws = volume_draws + batch_size;
    volume_accept_fsm = volume_accept_fsm + sum(safe_mask & informed_mask);
    volume_accept_eucl = volume_accept_eucl + sum(safe_mask & eucl_mask);

    n_global = round(eps_global * batch_size);
    global_mask = false(1, batch_size);
    if n_global > 0
        global_mask(randperm(batch_size, n_global)) = true;
    end

    sample_keep = safe_mask & (global_mask | informed_mask);
    S = S_all(:, sample_keep);
    if isempty(S)
        continue;
    end
    samples_accepted = samples_accepted + size(S, 2);

    tree_idx = 1:node_count;
    if isfinite(best_goal_cost)
        if ~meta.fallback
            h_nodes = trilinear_interp_batch(T_bwd, nodes(:, tree_idx), meta);
        else
            h_nodes = vecnorm(nodes(:, tree_idx) - x_goal, 2, 1);
        end
        keep_tree = (g_cost(tree_idx) + h_nodes) < best_goal_cost;
        keep_tree(1) = true;
        tree_idx = tree_idx(keep_tree);
    end

    n_tree = numel(tree_idx);
    n_samples = size(S, 2);
    if n_tree == 0 || n_samples == 0
        continue;
    end

    r_conn = dynamic_connection_radius(bounds, node_count + n_samples, min_r, max_r, eta_radius);
    tree_nodes = nodes(:, tree_idx);
    D = pairwise_dist(tree_nodes, S);
    edge_mask = D <= r_conn;

    [~, nearest_idx] = min(D, [], 1);
    nearest_lin = sub2ind(size(edge_mask), nearest_idx, 1:n_samples);
    edge_mask(nearest_lin) = true;

    [parent_local_tree, sample_local] = find(edge_mask);
    if isempty(parent_local_tree)
        continue;
    end

    parent_abs = tree_idx(parent_local_tree);
    P = nodes(:, parent_abs);
    C = S(:, sample_local);
    edge_len = reshape(D(sub2ind(size(D), parent_local_tree, sample_local)), [], 1);
    g_parent = reshape(g_cost(parent_abs), [], 1);

    if ~meta.fallback
        h_sample = reshape(trilinear_interp_batch(T_bwd, C, meta), [], 1);
    else
        h_sample = reshape(vecnorm(C - x_goal, 2, 1), [], 1);
    end

    if isfinite(best_goal_cost)
        edge_keep = isfinite(h_sample) & ((g_parent + edge_len + h_sample) < best_goal_cost);
        if ~any(edge_keep)
            continue;
        end
        parent_local_tree = parent_local_tree(edge_keep);
        parent_abs = parent_abs(edge_keep);
        sample_local = sample_local(edge_keep);
        P = P(:, edge_keep);
        C = C(:, edge_keep);
        edge_len = edge_len(edge_keep);
        g_parent = g_parent(edge_keep);
        h_sample = h_sample(edge_keep);
    end

    edge_key = g_parent + edge_len + h_sample;
    [~, order_edges] = sort(edge_key, 'ascend');
    parent_local_tree = parent_local_tree(order_edges);
    parent_abs = parent_abs(order_edges);
    sample_local = sample_local(order_edges);
    P = P(:, order_edges);
    C = C(:, order_edges);
    edge_len = edge_len(order_edges);
    g_parent = g_parent(order_edges);

    M = 0.5 * (P + C);
    risk_mid = reshape(risk_query(M), [], 1);
    if isinf(stealth_weight)
        edge_penalty = inf(size(edge_len));
        edge_penalty(risk_mid <= risk_limit) = 0;
        edge_w = edge_len + edge_penalty;
    else
        edge_w = edge_len + stealth_weight .* edge_len .* max(risk_mid, 0);
    end
    total_w = g_parent + edge_w;

    parent_costs = inf(n_tree, n_samples);
    lin_ps = sub2ind(size(parent_costs), parent_local_tree(:), sample_local(:));
    [lin_ps_u, iu_ps] = unique(lin_ps, 'stable');
    parent_costs(lin_ps_u) = total_w(iu_ps);
    [best_sample_cost, best_parent_local] = min(parent_costs, [], 1);
    candidate_mask = isfinite(best_sample_cost);
    if ~any(candidate_mask)
        continue;
    end

    cand_samples = S(:, candidate_mask);
    cand_parent_abs = tree_idx(best_parent_local(candidate_mask));
    cand_costs = best_sample_cost(candidate_mask);
    cand_parent_pts = nodes(:, cand_parent_abs);

    ok_edges = edge_is_safe_batch(cand_parent_pts, cand_samples, edge_check_samples, ...
        bounds, risk_query, risk_limit, terrain_query, point_clearance);
    if ~any(ok_edges)
        continue;
    end

    cand_samples = cand_samples(:, ok_edges);
    cand_parent_abs = cand_parent_abs(ok_edges);
    cand_costs = cand_costs(ok_edges);

    add_cap = max_nodes - node_count;
    if add_cap <= 0
        stop_reason = 'max_nodes';
        break;
    end

    n_add = min(add_cap, size(cand_samples, 2));
    if n_add <= 0
        continue;
    end

    old_count = node_count;
    new_indices = old_count + (1:n_add);
    nodes(:, new_indices) = cand_samples(:, 1:n_add);
    parents(new_indices) = uint32(cand_parent_abs(1:n_add));
    g_cost(new_indices) = cand_costs(1:n_add);
    node_count = old_count + n_add;

    if old_count >= 2 && n_add >= 1
        existing_idx = 2:old_count;
        new_idx = old_count + 1:node_count;

        D_rw = pairwise_dist(nodes(:, new_idx), nodes(:, existing_idx));
        rw_mask = D_rw <= (1.10 * r_conn);
        if any(rw_mask(:))
            [new_local, ex_local] = find(rw_mask);
            new_abs = new_idx(new_local);
            ex_abs = existing_idx(ex_local);

            P_rw = nodes(:, new_abs);
            C_rw = nodes(:, ex_abs);
            M_rw = 0.5 * (P_rw + C_rw);
            rw_risk = reshape(risk_query(M_rw), [], 1);
            rw_len = reshape(D_rw(sub2ind(size(D_rw), new_local, ex_local)), [], 1);
            g_new = reshape(g_cost(new_abs), [], 1);
            if isinf(stealth_weight)
                rw_penalty = inf(size(rw_len));
                rw_penalty(rw_risk <= risk_limit) = 0;
                rw_cost = g_new + rw_len + rw_penalty;
            else
                rw_cost = g_new + rw_len + stealth_weight .* rw_len .* max(rw_risk, 0);
            end

            rw_mat = inf(numel(new_idx), numel(existing_idx));
            lin_rw = sub2ind(size(rw_mat), new_local(:), ex_local(:));
            [lin_rw_u, iu_rw] = unique(lin_rw, 'stable');
            rw_mat(lin_rw_u) = rw_cost(iu_rw);
            [best_rw_cost, best_new_local] = min(rw_mat, [], 1);

            improve_mask = (best_rw_cost + 1e-9) < g_cost(existing_idx);
            if any(improve_mask)
                ex_targets = existing_idx(improve_mask);
                new_parents = new_idx(best_new_local(improve_mask));
                P_best = nodes(:, new_parents);
                C_best = nodes(:, ex_targets);
                ok_rewire = edge_is_safe_batch(P_best, C_best, edge_check_samples, ...
                    bounds, risk_query, risk_limit, terrain_query, point_clearance);
                if any(ok_rewire)
                    ex_ok = ex_targets(ok_rewire);
                    np_ok = new_parents(ok_rewire);
                    new_cost_ok = best_rw_cost(improve_mask);
                    new_cost_ok = new_cost_ok(ok_rewire);
                    parents(ex_ok) = uint32(np_ok);
                    g_cost(ex_ok) = new_cost_ok;
                end
            end
        end
    end

    D_goal = vecnorm(nodes(:, 1:node_count) - x_goal, 2, 1);
    goal_candidates = find(D_goal <= max(goal_radius, r_conn));
    if ~isempty(goal_candidates)
        P_goal = nodes(:, goal_candidates);
        G_goal = repmat(x_goal, 1, numel(goal_candidates));
        M_goal = 0.5 * (P_goal + G_goal);
        goal_risk = reshape(risk_query(M_goal), [], 1);
        goal_len = reshape(D_goal(goal_candidates), [], 1);
        g_goal = reshape(g_cost(goal_candidates), [], 1);
        if isinf(stealth_weight)
            goal_penalty = inf(size(goal_len));
            goal_penalty(goal_risk <= risk_limit) = 0;
            goal_total = g_goal + goal_len + goal_penalty;
        else
            goal_total = g_goal + goal_len + stealth_weight .* goal_len .* max(goal_risk, 0);
        end

        better = goal_total < (best_goal_cost - 1e-9);
        if any(better)
            cand_nodes = goal_candidates(better);
            cand_costs = goal_total(better);
            P_try = nodes(:, cand_nodes);
            G_try = repmat(x_goal, 1, numel(cand_nodes));
            ok_goal = edge_is_safe_batch(P_try, G_try, edge_check_samples, ...
                bounds, risk_query, risk_limit, terrain_query, point_clearance);
            if any(ok_goal)
                valid_nodes = cand_nodes(ok_goal);
                valid_costs = cand_costs(ok_goal);
                [best_local, best_k] = min(valid_costs);
                best_goal_cost = best_local;
                goal_parent = uint32(valid_nodes(best_k));
                if isnan(first_solution_batch)
                    first_solution_batch = batch_idx;
                    first_solution_time = toc(t_start);
                end
            end
        end
    end

    if goal_parent > 0 && batch_idx >= first_solution_batch + post_solution_batches
        stop_reason = 'post_solution_budget';
        break;
    end
end

planning_time = toc(t_start);
success = goal_parent > 0 || isfinite(warmstart_bound_cost);
if success
    if goal_parent > 0
        path = extract_path(nodes, parents, goal_parent, x_start, x_goal, max_nodes);
    else
        path = warmstart_world;
    end
else
    path = [];
end

if volume_draws > 0 && volume_accept_eucl > 0
    volume_ratio = (volume_accept_fsm / volume_draws) / (volume_accept_eucl / volume_draws);
else
    volume_ratio = NaN;
end

info = struct();
info.success = success;
info.iterations = batches_run;
info.tree_size = node_count;
info.planning_time = planning_time;
info.first_solution_batch = first_solution_batch;
info.first_solution_time = first_solution_time;
info.path_cost = best_goal_cost;
info.samples_generated = samples_generated;
info.samples_accepted = samples_accepted;
info.global_sampling_epsilon = eps_global;
info.connection_radius_last = dynamic_connection_radius(bounds, max(node_count, 2), min_r, max_r, eta_radius);
info.stop_reason = stop_reason;
info.max_batches = max_batches;
info.post_solution_batches = post_solution_batches;
info.risk_limit = risk_limit;
info.fsm_used = ~meta.fallback;
info.fsm_time_s = get_meta_time(meta);
info.volume_ratio = volume_ratio;
info.warmstart_cost = meta.warmstart_cost;
info.warmstart_bound_cost = warmstart_bound_cost;
info.fsm_meta = meta;
info.fsm_T_fwd = T_fwd;
info.fsm_T_bwd = T_bwd;
info.env = env_meta;
cost = best_goal_cost;
end

function [terrain_input, terrain_query, radar_los_map, risk_query, vis_threshold, bounds, point_clearance, env_meta] = ...
    parse_direct_env(env, params, x_start, x_goal)
if ~isstruct(env)
    error('fsm_bit_star_planner_direct:InvalidEnv', 'env must be a struct.');
end

point_clearance = max(get_opt(params, 'fsm_point_clearance', get_opt(params, 'min_clearance', 5.0)), 0.0);
env_meta = struct();

terrain_source = [];
terrain_fields = {'terrain', 'terrain_map', 'fsm_terrain_map'};
for i = 1:numel(terrain_fields)
    if isfield(env, terrain_fields{i}) && ~isempty(env.(terrain_fields{i}))
        terrain_source = env.(terrain_fields{i});
        break;
    end
end
if isempty(terrain_source)
    error('fsm_bit_star_planner_direct:MissingTerrain', 'env must include terrain or terrain_map.');
end

terrain_input = terrain_source_to_struct(terrain_source, env);
terrain_query = make_terrain_query(terrain_input);
env_meta.terrain_bounds = [terrain_input.N_vec(1), terrain_input.N_vec(end), ...
                           terrain_input.E_vec(1), terrain_input.E_vec(end)];

if isfield(env, 'radar_los_map') && ~isempty(env.radar_los_map)
    radar_los_map = logical(env.radar_los_map);
elseif isfield(env, 'risk_grid') && ~isempty(env.risk_grid)
    radar_los_map = [];
else
    error('fsm_bit_star_planner_direct:MissingThreat', ...
        'env must include radar_los_map or risk_grid.');
end

if isfield(env, 'visibility_threshold') && ~isempty(env.visibility_threshold)
    vis_threshold = env.visibility_threshold;
else
    vis_threshold = 0.5;
end

if ~isfield(env, 'N_vec') || ~isfield(env, 'E_vec') || ~isfield(env, 'alt_vec')
    error('fsm_bit_star_planner_direct:MissingAxes', ...
        'env must include N_vec, E_vec, and alt_vec.');
end

N_vec = env.N_vec(:)';
E_vec = env.E_vec(:)';
alt_vec = env.alt_vec(:)';

if isfield(env, 'risk_grid') && ~isempty(env.risk_grid)
    V = env.risk_grid;
    if isequal(size(V), [numel(E_vec), numel(N_vec), numel(alt_vec)])
        F = griddedInterpolant({E_vec, N_vec, alt_vec}, V, 'linear', 'none');
        risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);
        if isempty(radar_los_map)
            radar_los_map = V > vis_threshold;
        end
    else
        error('fsm_bit_star_planner_direct:RiskGridSizeMismatch', ...
            'env.risk_grid size must match [numel(E_vec) x numel(N_vec) x numel(alt_vec)].');
    end
elseif ~isempty(radar_los_map)
    if isequal(size(radar_los_map), [numel(E_vec), numel(N_vec), numel(alt_vec)])
        V = double(radar_los_map);
        F = griddedInterpolant({E_vec, N_vec, alt_vec}, V, 'nearest', 'none');
        risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);
    else
        error('fsm_bit_star_planner_direct:LOSSizeMismatch', ...
            'env.radar_los_map size must match [numel(E_vec) x numel(N_vec) x numel(alt_vec)].');
    end
end

if isfield(env, 'bounds_world') && ~isempty(env.bounds_world)
    bounds = env.bounds_world(:)';
else
    bounds = [N_vec(1), N_vec(end), E_vec(1), E_vec(end), alt_vec(1), alt_vec(end)];
end

bounds(1) = min([bounds(1), x_start(1), x_goal(1)]);
bounds(2) = max([bounds(2), x_start(1), x_goal(1)]);
bounds(3) = min([bounds(3), x_start(2), x_goal(2)]);
bounds(4) = max([bounds(4), x_start(2), x_goal(2)]);
bounds(5) = min([bounds(5), x_start(3), x_goal(3)]);
bounds(6) = max([bounds(6), x_start(3), x_goal(3)]);

env_meta.grid_size = [numel(N_vec), numel(E_vec), numel(alt_vec)];
env_meta.bounds = bounds;
end

function terrain_input = terrain_source_to_struct(src, env)
terrain_input = struct();
if isa(src, 'terrain_map')
    terrain_input.Z = src.Z;
    terrain_input.N_vec = src.N_vec(:)';
    terrain_input.E_vec = src.E_vec(:)';
    terrain_input.dx = mean(diff(src.N_vec));
    terrain_input.dy = mean(diff(src.E_vec));
elseif isstruct(src)
    terrain_input = src;
    if ~isfield(terrain_input, 'N_vec') && isfield(terrain_input, 'x_vec')
        terrain_input.N_vec = terrain_input.x_vec(:)';
    end
    if ~isfield(terrain_input, 'E_vec') && isfield(terrain_input, 'y_vec')
        terrain_input.E_vec = terrain_input.y_vec(:)';
    end
    if ~isfield(terrain_input, 'dx') && isfield(terrain_input, 'N_vec')
        terrain_input.dx = mean(diff(terrain_input.N_vec));
    end
    if ~isfield(terrain_input, 'dy') && isfield(terrain_input, 'E_vec')
        terrain_input.dy = mean(diff(terrain_input.E_vec));
    end
else
    error('fsm_bit_star_planner_direct:InvalidTerrainContext', ...
        'Unsupported terrain source type: %s', class(src));
end

if isfield(env, 'alt_vec') && ~isempty(env.alt_vec)
    terrain_input.alt_vec = env.alt_vec(:)';
end
terrain_input.x_vec = terrain_input.N_vec(:)';
terrain_input.y_vec = terrain_input.E_vec(:)';
end

function F = make_terrain_query(terrain_input)
N_vec = terrain_input.N_vec(:)';
E_vec = terrain_input.E_vec(:)';
Z = terrain_input.Z;
if isequal(size(Z), [numel(E_vec), numel(N_vec)])
    Z_NE = double(Z.');
elseif isequal(size(Z), [numel(N_vec), numel(E_vec)])
    Z_NE = double(Z);
else
    error('fsm_bit_star_planner_direct:TerrainSizeMismatch', ...
        'terrain_input.Z size does not match terrain axes.');
end
F = griddedInterpolant({N_vec, E_vec}, Z_NE, 'linear', 'nearest');
end

function tf = points_clear_of_terrain(P, terrain_query, clearance)
if isempty(P)
    tf = false(1, 0);
    return;
end
h = terrain_query(P(1, :).', P(2, :).');
h = reshape(h, 1, []);
tf = isfinite(h) & (P(3, :) >= (h + clearance));
end

function cost = evaluate_world_path_cost(path_world, risk_query, risk_limit, stealth_weight, point_clearance, terrain_query)
if isempty(path_world) || size(path_world, 2) < 2
    cost = Inf;
    return;
end

cost = 0;
for i = 1:(size(path_world, 2) - 1)
    p0 = path_world(:, i);
    p1 = path_world(:, i + 1);
    seg_mid = 0.5 * (p0 + p1);
    seg_len = norm(p1 - p0);
    if ~all(points_clear_of_terrain([p0, p1, seg_mid], terrain_query, point_clearance))
        cost = Inf;
        return;
    end
    risk_mid = risk_query(seg_mid);
    if ~isfinite(risk_mid) || risk_mid > risk_limit
        cost = Inf;
        return;
    end
    if isinf(stealth_weight)
        cost = cost + seg_len;
    else
        cost = cost + seg_len + stealth_weight * seg_len * max(risk_mid, 0);
    end
end
end

function S = sample_uniform_bounds(bounds, n)
if n <= 0
    S = zeros(3, 0);
    return;
end
S = zeros(3, n);
S(1, :) = bounds(1) + rand(1, n) .* (bounds(2) - bounds(1));
S(2, :) = bounds(3) + rand(1, n) .* (bounds(4) - bounds(3));
S(3, :) = bounds(5) + rand(1, n) .* (bounds(6) - bounds(5));
end

function tf = points_in_bounds(P, bounds)
tf = (P(1, :) >= bounds(1)) & (P(1, :) <= bounds(2)) & ...
     (P(2, :) >= bounds(3)) & (P(2, :) <= bounds(4)) & ...
     (P(3, :) >= bounds(5)) & (P(3, :) <= bounds(6));
end

function D = pairwise_dist(A, B)
AA = sum(A .^ 2, 1)';
BB = sum(B .^ 2, 1);
D2 = max(0.0, AA + BB - 2.0 .* (A' * B));
D = sqrt(D2);
end

function r = dynamic_connection_radius(bounds, n, r_min, r_max, eta)
vol = max((bounds(2) - bounds(1)) * (bounds(4) - bounds(3)) * (bounds(6) - bounds(5)), 1e-9);
zeta3 = 4.0 * pi / 3.0;
gamma = 2.0 * (1.0 + 1.0 / 3.0)^(1.0 / 3.0) * (vol / zeta3)^(1.0 / 3.0);
n_eff = max(2, n);
r_asym = eta * gamma * (log(n_eff) / n_eff)^(1.0 / 3.0);
r = min(r_max, max(r_min, r_asym));
end

function ok = edge_is_safe_batch(P0, P1, n_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance)
n_edges = size(P0, 2);
if n_edges == 0
    ok = false(1, 0);
    return;
end
t = linspace(0, 1, n_samples);
if numel(t) > 2
    t = t(2:end-1);
end
n_t = numel(t);
if n_t == 0
    ok = true(1, n_edges);
    return;
end

T = reshape(t, 1, n_t, 1);
P0r = reshape(P0, 3, 1, n_edges);
P1r = reshape(P1, 3, 1, n_edges);
Pts = P0r .* (1 - T) + P1r .* T;
Pts_flat = reshape(Pts, 3, []);
inb = points_in_bounds(Pts_flat, bounds);
risk = risk_query(Pts_flat);
terrain_ok = points_clear_of_terrain(Pts_flat, terrain_query, point_clearance);
safe = inb & isfinite(risk) & (risk <= risk_limit) & terrain_ok;
safe = reshape(safe, n_t, n_edges);
ok = all(safe, 1);
end

function path = extract_path(nodes, parents, goal_parent, start, goal, max_nodes)
chain = zeros(1, max_nodes, 'uint32');
count = 0;
idx = goal_parent;
while idx > 0 && count < max_nodes
    count = count + 1;
    chain(count) = idx;
    idx = parents(idx);
end
if count == 0
    path = [];
    return;
end
chain = fliplr(chain(1:count));
n_mid = numel(chain);
path = zeros(3, n_mid + 2);
path(:, 1) = start;
path(:, 2:n_mid+1) = nodes(:, chain);
path(:, n_mid+2) = goal;
d = vecnorm(diff(path, 1, 2), 2, 1);
keep = [true, d > 1e-9];
path = path(:, keep);
end

function vals = trilinear_interp_batch(T, P_world, meta)
if isempty(P_world)
    vals = zeros(1, 0);
    return;
end
vals = inf(1, size(P_world, 2));
for i = 1:size(P_world, 2)
    vals(i) = trilinear_interp(T, P_world(:, i), meta);
end
end

function val = trilinear_interp(T, x_world, meta)
if isempty(T) || isempty(meta) || ~isfield(meta, 'grid_x') || isempty(meta.grid_x)
    val = Inf;
    return;
end
x_world = x_world(:);
fx = (x_world(1) - meta.grid_x(1)) / meta.dx + 1;
fy = (x_world(2) - meta.grid_y(1)) / meta.dy + 1;
fz = (x_world(3) - meta.grid_z(1)) / meta.dz + 1;

sz = size(T);
if fx < 1 || fy < 1 || fz < 1 || fx > sz(1) || fy > sz(2) || fz > sz(3)
    val = Inf;
    return;
end

i0 = max(1, min(sz(1) - 1, floor(fx)));
j0 = max(1, min(sz(2) - 1, floor(fy)));
k0 = max(1, min(sz(3) - 1, floor(fz)));
i1 = min(sz(1), i0 + 1);
j1 = min(sz(2), j0 + 1);
k1 = min(sz(3), k0 + 1);

tx = max(0, min(1, fx - i0));
ty = max(0, min(1, fy - j0));
tz = max(0, min(1, fz - k0));

c000 = double(T(i0, j0, k0));
c100 = double(T(i1, j0, k0));
c010 = double(T(i0, j1, k0));
c110 = double(T(i1, j1, k0));
c001 = double(T(i0, j0, k1));
c101 = double(T(i1, j0, k1));
c011 = double(T(i0, j1, k1));
c111 = double(T(i1, j1, k1));

c00 = (1 - tx) * c000 + tx * c100;
c10 = (1 - tx) * c010 + tx * c110;
c01 = (1 - tx) * c001 + tx * c101;
c11 = (1 - tx) * c011 + tx * c111;
c0 = (1 - ty) * c00 + ty * c10;
c1 = (1 - ty) * c01 + ty * c11;
val = (1 - tz) * c0 + tz * c1;
end

function x_ned = world_to_ned(x_world)
x_ned = [x_world(1); x_world(2); -x_world(3)];
end

function t = get_meta_time(meta)
if isfield(meta, 'fsm_time_s')
    t = meta.fsm_time_s;
else
    t = 0;
end
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
