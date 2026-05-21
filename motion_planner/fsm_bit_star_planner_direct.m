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

ctx = build_direct_planner_context(x_start, x_goal, env, params);
terrain_input = ctx.terrain_input;
terrain_query = ctx.terrain_query;
radar_los_map = ctx.radar_los_map;
risk_query = ctx.risk_query;
bounds = ctx.bounds;
point_clearance = ctx.point_clearance;
env_meta = ctx.env_meta;
risk_limit = ctx.risk_limit;
altitude_cost_weight = ctx.altitude_cost_weight;
preferred_agl = ctx.preferred_agl;

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
    if get_opt(params, 'force_euclidean_heuristic', false)
        meta.fallback = true;
    else
        [T_fwd, T_bwd, ~, meta] = compute_fsm_heuristic( ...
            terrain_input, radar_los_map, ctx.world_to_ned(x_start), ctx.world_to_ned(x_goal), fsm_params);
    end
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

default_goal_radius = max(5.0, 0.02 * norm(bounds([2 4 6]) - bounds([1 3 5])));
goal_radius = get_opt(params, 'goal_radius', default_goal_radius);
min_r = get_opt(params, 'min_connection_radius', max(5.0, 0.012 * norm(bounds([2 4 6]) - bounds([1 3 5]))));
max_r = get_opt(params, 'max_connection_radius', max(40, 0.08 * norm(bounds([2 4 6]) - bounds([1 3 5]))));
eta_radius = get_opt(params, 'eta_radius', 1.8);

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
    warmstart_bound_cost = evaluate_path_cost_strict(warmstart_world, bounds, risk_query, risk_limit, ...
        terrain_query, point_clearance, edge_check_samples, stealth_weight, altitude_cost_weight, preferred_agl);
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
    safe_mask = safe_mask & ctx.points_clear_of_terrain(S_all);

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
    alt_cost = compute_altitude_edge_cost_batch(M, edge_len, terrain_query, altitude_cost_weight, preferred_agl);
    if isinf(stealth_weight)
        edge_penalty = inf(size(edge_len));
        edge_penalty(risk_mid <= risk_limit) = 0;
        edge_w = edge_len + edge_penalty + alt_cost;
    else
        edge_w = edge_len + stealth_weight .* edge_len .* max(risk_mid, 0) + alt_cost;
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
            rw_alt_cost = compute_altitude_edge_cost_batch(M_rw, rw_len, terrain_query, altitude_cost_weight, preferred_agl);
            if isinf(stealth_weight)
                rw_penalty = inf(size(rw_len));
                rw_penalty(rw_risk <= risk_limit) = 0;
                rw_cost = g_new + rw_len + rw_penalty + rw_alt_cost;
            else
                rw_cost = g_new + rw_len + stealth_weight .* rw_len .* max(rw_risk, 0) + rw_alt_cost;
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
        goal_alt_cost = compute_altitude_edge_cost_batch(M_goal, goal_len, terrain_query, altitude_cost_weight, preferred_agl);
        if isinf(stealth_weight)
            goal_penalty = inf(size(goal_len));
            goal_penalty(goal_risk <= risk_limit) = 0;
            goal_total = g_goal + goal_len + goal_penalty + goal_alt_cost;
        else
            goal_total = g_goal + goal_len + stealth_weight .* goal_len .* max(goal_risk, 0) + goal_alt_cost;
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

function total_cost = evaluate_path_cost_strict(path_world, bounds, risk_query, risk_limit, terrain_query, point_clearance, edge_check_samples, stealth_weight, altitude_cost_weight, preferred_agl)
if isempty(path_world) || size(path_world, 2) < 2
    total_cost = Inf;
    return;
end

total_cost = 0;
for i = 1:(size(path_world, 2) - 1)
    p0 = path_world(:, i);
    p1 = path_world(:, i + 1);
    if ~edge_is_safe_batch(p0, p1, edge_check_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance)
        total_cost = Inf;
        return;
    end

    seg_len = norm(p1 - p0);
    seg_mid = 0.5 * (p0 + p1);
    risk_mid = risk_query(seg_mid);
    if ~isfinite(risk_mid) || risk_mid > risk_limit
        total_cost = Inf;
        return;
    end
    cost_alt = compute_altitude_edge_cost_batch(seg_mid, seg_len, terrain_query, altitude_cost_weight, preferred_agl);

    if isinf(stealth_weight)
        total_cost = total_cost + seg_len + cost_alt;
    else
        total_cost = total_cost + seg_len + stealth_weight * seg_len * max(risk_mid, 0) + cost_alt;
    end
end
end

function cost_alt = compute_altitude_edge_cost_batch(P, seg_len, terrain_query, altitude_cost_weight, preferred_agl)
if altitude_cost_weight <= 0
    cost_alt = zeros(size(seg_len));
    return;
end
terrain_h = terrain_query(P(1, :).', P(2, :).');
agl = P(3, :)' - terrain_h;
agl_scale = max(preferred_agl, 1);
err_norm = (agl - preferred_agl) ./ agl_scale;
cost_alt = altitude_cost_weight .* seg_len .* (err_norm .^ 2);
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
h = terrain_query(Pts_flat(1, :).', Pts_flat(2, :).');
h = reshape(h, 1, []);
terrain_ok = isfinite(h) & (Pts_flat(3, :) >= (h + point_clearance));
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
