function [path, cost, info] = plain_radar_bit_star_planner(x_start, x_goal, env, params)
%PLAIN_RADAR_BIT_STAR_PLANNER Corridor-free BIT* with soft radar edge cost.
%
% Inputs
%   x_start [3x1] start [N; E; Alt]
%   x_goal  [3x1] goal  [N; E; Alt]
%   env     struct with:
%     .risk_grid    [nE x nN x nAlt] continuous cost in [0,1]
%     .N_vec        North axis for risk_grid
%     .E_vec        East axis for risk_grid
%     .alt_vec      altitude axis for risk_grid
%     .bounds_world [Nmin Nmax Emin Emax Altmin Altmax] optional
%     .terrain_map  terrain_map object or .terrain struct optional
%   params optional fields:
%     .stealth_weight          radar cost multiplier W (default 4)
%     .point_clearance         hard terrain clearance [m] (default 10)
%     .edge_sample_spacing     edge cost/check spacing [m] (default 40)
%     .max_edge_samples        cap per edge check (default 45)
%     .max_batches             default 100
%     .batch_size              default 320
%     .max_nodes               default 18000
%     .post_solution_batches   default 12
%     .max_time                max planning time [s] (default inf)
%
% The only hard feasibility checks are map bounds and terrain clearance.
% Radar exposure is never used as a strict rejection rule; it only changes
% edge cost as edge_length * (1 + W * mean(edge_risk)).

if nargin < 4
    params = struct();
end

validateattributes(x_start, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_start', 1);
validateattributes(x_goal, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_goal', 2);
x_start = x_start(:);
x_goal = x_goal(:);

ctx = build_context(env, params);
bounds = ctx.bounds;
risk_query = ctx.risk_query;
terrain_query = ctx.terrain_query;

max_batches = max(1, round(get_opt(params, 'max_batches', 100)));
batch_size = max(32, round(get_opt(params, 'batch_size', 320)));
max_nodes = max(batch_size + 2, round(get_opt(params, 'max_nodes', 18000)));
post_solution_batches = max(0, round(get_opt(params, 'post_solution_batches', 12)));
max_time = get_opt(params, 'max_time', inf);
stealth_weight = max(0, get_opt(params, 'stealth_weight', 4.0));
point_clearance = max(0, get_opt(params, 'point_clearance', 10.0));
edge_sample_spacing = max(1, get_opt(params, 'edge_sample_spacing', 40.0));
max_edge_samples = max(3, round(get_opt(params, 'max_edge_samples', 45)));

diag_len = norm(bounds([2 4 6]) - bounds([1 3 5]));
goal_radius = get_opt(params, 'goal_radius', max(20.0, 0.018 * diag_len));
min_r = get_opt(params, 'min_connection_radius', max(20.0, 0.012 * diag_len));
max_r = get_opt(params, 'max_connection_radius', max(80.0, 0.075 * diag_len));
eta_radius = get_opt(params, 'eta_radius', 1.8);

info = struct();
info.success = false;
info.stop_reason = 'not_started';
info.planning_time = 0;
info.samples_generated = 0;
info.samples_accepted = 0;
info.tree_size = 0;
info.iterations = 0;
info.first_solution_batch = NaN;
info.first_solution_time = NaN;
info.path_cost = Inf;
info.mean_path_risk = NaN;
info.planner_name = 'plain_radar_bit_star_planner';

if ~points_in_bounds(x_start, bounds) || ~points_in_bounds(x_goal, bounds)
    path = [];
    cost = Inf;
    info.stop_reason = 'endpoint_out_of_bounds';
    return;
end
if ~points_clear_of_terrain(x_start, terrain_query, point_clearance) || ...
        ~points_clear_of_terrain(x_goal, terrain_query, point_clearance)
    path = [];
    cost = Inf;
    info.stop_reason = 'endpoint_terrain_clearance';
    return;
end

nodes = nan(3, max_nodes);
parents = zeros(1, max_nodes, 'uint32');
g_cost = inf(1, max_nodes);

nodes(:, 1) = x_start;
parents(1) = uint32(0);
g_cost(1) = 0.0;
node_count = 1;

best_goal_cost = inf;
goal_parent = uint32(0);
first_solution_batch = NaN;
first_solution_time = NaN;
samples_generated = 0;
samples_accepted = 0;

t_start = tic;
batches_run = 0;
stop_reason = 'max_batches_reached';

for batch_idx = 1:max_batches
    if toc(t_start) > max_time
        stop_reason = 'max_time';
        break;
    end

    batches_run = batch_idx;
    if isfinite(best_goal_cost)
        S_all = sample_informed_ellipsoid(x_start, x_goal, best_goal_cost, bounds, batch_size);
    else
        S_all = sample_uniform_bounds(bounds, batch_size);
    end
    samples_generated = samples_generated + size(S_all, 2);

    keep_samples = points_in_bounds(S_all, bounds) & ...
        points_clear_of_terrain(S_all, terrain_query, point_clearance);
    S = S_all(:, keep_samples);
    if isempty(S)
        continue;
    end
    samples_accepted = samples_accepted + size(S, 2);

    tree_idx = 1:node_count;
    if isfinite(best_goal_cost)
        h_nodes = vecnorm(nodes(:, tree_idx) - x_goal, 2, 1);
        tree_keep = (g_cost(tree_idx) + h_nodes) < best_goal_cost;
        tree_keep(1) = true;
        tree_idx = tree_idx(tree_keep);
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

    if isfinite(best_goal_cost)
        h_sample = reshape(vecnorm(C - x_goal, 2, 1), [], 1);
        admissible_keep = (g_parent + edge_len + h_sample) < best_goal_cost;
        if ~any(admissible_keep)
            continue;
        end
        parent_local_tree = parent_local_tree(admissible_keep);
        sample_local = sample_local(admissible_keep);
        P = P(:, admissible_keep);
        C = C(:, admissible_keep);
        g_parent = g_parent(admissible_keep);
    end

    [ok_edges, edge_w] = evaluate_edges(P, C, bounds, risk_query, terrain_query, ...
        point_clearance, stealth_weight, edge_sample_spacing, max_edge_samples);
    if ~any(ok_edges)
        continue;
    end

    parent_local_tree = parent_local_tree(ok_edges);
    sample_local = sample_local(ok_edges);
    edge_w = edge_w(ok_edges);
    g_parent = g_parent(ok_edges);
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
    cand_parents = tree_idx(best_parent_local(candidate_mask));
    cand_costs = best_sample_cost(candidate_mask);

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
    parents(new_indices) = uint32(cand_parents(1:n_add));
    g_cost(new_indices) = cand_costs(1:n_add);
    node_count = old_count + n_add;

    if old_count >= 2 && n_add >= 1
        existing_idx = 2:old_count;
        new_idx = old_count + 1:node_count;
        D_rw = pairwise_dist(nodes(:, new_idx), nodes(:, existing_idx));
        rw_mask = D_rw <= 1.10 * r_conn;
        if any(rw_mask(:))
            [new_local, ex_local] = find(rw_mask);
            new_abs = new_idx(new_local);
            ex_abs = existing_idx(ex_local);
            P_rw = nodes(:, new_abs);
            C_rw = nodes(:, ex_abs);
            rw_len = reshape(D_rw(sub2ind(size(D_rw), new_local, ex_local)), [], 1);
            g_new = reshape(g_cost(new_abs), [], 1);

            if isfinite(best_goal_cost)
                h_ex = reshape(vecnorm(C_rw - x_goal, 2, 1), [], 1);
                rw_keep = (g_new + rw_len + h_ex) < best_goal_cost;
                new_local = new_local(rw_keep);
                ex_local = ex_local(rw_keep);
                new_abs = new_abs(rw_keep);
                P_rw = P_rw(:, rw_keep);
                C_rw = C_rw(:, rw_keep);
                g_new = g_new(rw_keep);
            end

            if ~isempty(new_abs)
                [ok_rw, rw_edge_w] = evaluate_edges(P_rw, C_rw, bounds, risk_query, terrain_query, ...
                    point_clearance, stealth_weight, edge_sample_spacing, max_edge_samples);
                if any(ok_rw)
                    rw_cost = g_new(ok_rw) + rw_edge_w(ok_rw);
                    rw_mat = inf(numel(new_idx), numel(existing_idx));
                    lin_rw = sub2ind(size(rw_mat), new_local(ok_rw), ex_local(ok_rw));
                    [lin_rw_u, iu_rw] = unique(lin_rw, 'stable');
                    rw_mat(lin_rw_u) = rw_cost(iu_rw);
                    [best_rw_cost, best_new_local] = min(rw_mat, [], 1);
                    improve_mask = (best_rw_cost + 1e-9) < g_cost(existing_idx);
                    if any(improve_mask)
                        ex_targets = existing_idx(improve_mask);
                        new_parents = new_idx(best_new_local(improve_mask));
                        parents(ex_targets) = uint32(new_parents);
                        g_cost(ex_targets) = best_rw_cost(improve_mask);
                    end
                end
            end
        end
    end

    D_goal = vecnorm(nodes(:, 1:node_count) - x_goal, 2, 1);
    goal_candidates = find(D_goal <= max(goal_radius, r_conn));
    if ~isempty(goal_candidates)
        P_goal = nodes(:, goal_candidates);
        G_goal = repmat(x_goal, 1, numel(goal_candidates));
        goal_len = reshape(D_goal(goal_candidates), [], 1);
        g_goal = reshape(g_cost(goal_candidates), [], 1);
        if isfinite(best_goal_cost)
            goal_keep = (g_goal + goal_len) < best_goal_cost;
            P_goal = P_goal(:, goal_keep);
            G_goal = G_goal(:, goal_keep);
            goal_candidates = goal_candidates(goal_keep);
            g_goal = g_goal(goal_keep);
        end
        if ~isempty(goal_candidates)
            [ok_goal, goal_edge_w] = evaluate_edges(P_goal, G_goal, bounds, risk_query, terrain_query, ...
                point_clearance, stealth_weight, edge_sample_spacing, max_edge_samples);
            if any(ok_goal)
                goal_total = g_goal(ok_goal) + goal_edge_w(ok_goal);
                valid_nodes = goal_candidates(ok_goal);
                [best_local, best_k] = min(goal_total);
                if best_local < best_goal_cost
                    best_goal_cost = best_local;
                    goal_parent = uint32(valid_nodes(best_k));
                    if isnan(first_solution_batch)
                        first_solution_batch = batch_idx;
                        first_solution_time = toc(t_start);
                    end
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
success = goal_parent > 0;
if success
    path = extract_path(nodes, parents, goal_parent, x_goal, max_nodes);
    cost = best_goal_cost;
else
    path = [];
    cost = Inf;
end

info.success = success;
info.stop_reason = stop_reason;
info.planning_time = planning_time;
info.samples_generated = samples_generated;
info.samples_accepted = samples_accepted;
info.tree_size = node_count;
info.iterations = batches_run;
info.first_solution_batch = first_solution_batch;
info.first_solution_time = first_solution_time;
info.path_cost = best_goal_cost;
if success
    info.mean_path_risk = mean_path_risk(path, risk_query, edge_sample_spacing, max_edge_samples);
end
end

function ctx = build_context(env, params)
required = {'risk_grid', 'N_vec', 'E_vec', 'alt_vec'};
for i = 1:numel(required)
    if ~isfield(env, required{i}) || isempty(env.(required{i}))
        error('plain_radar_bit_star_planner:InvalidEnv', 'env.%s is required.', required{i});
    end
end

N_vec = env.N_vec(:)';
E_vec = env.E_vec(:)';
alt_vec = env.alt_vec(:)';
V = double(env.risk_grid);
if ~isequal(size(V), [numel(E_vec), numel(N_vec), numel(alt_vec)])
    error('plain_radar_bit_star_planner:RiskSizeMismatch', ...
        'env.risk_grid must be [numel(E_vec) x numel(N_vec) x numel(alt_vec)].');
end
F = griddedInterpolant({E_vec, N_vec, alt_vec}, V, 'linear', 'nearest');
risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);

if isfield(env, 'bounds_world') && ~isempty(env.bounds_world)
    bounds = env.bounds_world(:)';
else
    bounds = [N_vec(1), N_vec(end), E_vec(1), E_vec(end), alt_vec(1), alt_vec(end)];
end

if isfield(params, 'terrain_query') && isa(params.terrain_query, 'function_handle')
    terrain_query = params.terrain_query;
elseif isfield(env, 'terrain_query') && isa(env.terrain_query, 'function_handle')
    terrain_query = env.terrain_query;
elseif isfield(env, 'terrain_map') && ~isempty(env.terrain_map)
    tm = env.terrain_map;
    terrain_query = @(N, E) tm.get_height(N, E);
elseif isfield(env, 'terrain') && ~isempty(env.terrain)
    terrain_query = make_terrain_query(env.terrain);
else
    terrain_query = @(N, E) -inf(size(N));
end

ctx = struct();
ctx.bounds = bounds;
ctx.risk_query = risk_query;
ctx.terrain_query = terrain_query;
end

function terrain_query = make_terrain_query(src)
N_vec = src.N_vec(:)';
E_vec = src.E_vec(:)';
Z = double(src.Z);
if isequal(size(Z), [numel(N_vec), numel(E_vec)])
    Z = Z';
end
F = griddedInterpolant({E_vec, N_vec}, Z, 'linear', 'nearest');
terrain_query = @(N, E) F(E, N);
end

function S = sample_uniform_bounds(bounds, n)
S = zeros(3, n);
S(1, :) = bounds(1) + rand(1, n) .* (bounds(2) - bounds(1));
S(2, :) = bounds(3) + rand(1, n) .* (bounds(4) - bounds(3));
S(3, :) = bounds(5) + rand(1, n) .* (bounds(6) - bounds(5));
end

function S = sample_informed_ellipsoid(x_start, x_goal, c_best, bounds, n)
c_min = norm(x_goal - x_start);
if ~isfinite(c_best) || c_best <= c_min + 1e-9
    S = sample_uniform_bounds(bounds, n);
    return;
end

center = 0.5 * (x_start + x_goal);
a1 = (x_goal - x_start) / max(c_min, 1e-12);
C = basis_from_axis(a1);
r1 = 0.5 * c_best;
r2 = 0.5 * sqrt(max(c_best^2 - c_min^2, 0));
L = diag([r1, r2, r2]);

S = zeros(3, n);
count = 0;
draws = 0;
max_draws = max(10 * n, 200);
while count < n && draws < max_draws
    draws = draws + n;
    U = sample_unit_ball(n);
    X = center + C * L * U;
    inb = points_in_bounds(X, bounds);
    n_take = min(n - count, sum(inb));
    if n_take > 0
        X_keep = X(:, inb);
        S(:, count + (1:n_take)) = X_keep(:, 1:n_take);
        count = count + n_take;
    end
end
if count < n
    S(:, count+1:n) = sample_uniform_bounds(bounds, n - count);
end
end

function U = sample_unit_ball(n)
dirs = randn(3, n);
dirs = dirs ./ max(vecnorm(dirs, 2, 1), 1e-12);
rho = rand(1, n) .^ (1 / 3);
U = dirs .* rho;
end

function C = basis_from_axis(a1)
a1 = a1(:) / max(norm(a1), 1e-12);
if abs(dot(a1, [1; 0; 0])) < 0.9
    tmp = [1; 0; 0];
else
    tmp = [0; 1; 0];
end
a2 = cross(a1, tmp);
a2 = a2 / max(norm(a2), 1e-12);
a3 = cross(a1, a2);
a3 = a3 / max(norm(a3), 1e-12);
C = [a1, a2, a3];
end

function tf = points_in_bounds(P, bounds)
tf = (P(1, :) >= bounds(1)) & (P(1, :) <= bounds(2)) & ...
     (P(2, :) >= bounds(3)) & (P(2, :) <= bounds(4)) & ...
     (P(3, :) >= bounds(5)) & (P(3, :) <= bounds(6));
end

function tf = points_clear_of_terrain(P, terrain_query, clearance)
h = terrain_query(P(1, :).', P(2, :).');
h = reshape(h, 1, []);
tf = isfinite(h) & (P(3, :) >= h + clearance);
end

function [ok, edge_w, avg_risk] = evaluate_edges(P0, P1, bounds, risk_query, terrain_query, ...
    point_clearance, stealth_weight, edge_sample_spacing, max_edge_samples)
n_edges = size(P0, 2);
ok = false(n_edges, 1);
edge_w = inf(n_edges, 1);
avg_risk = nan(n_edges, 1);
if n_edges == 0
    return;
end

seg = P1 - P0;
seg_len = vecnorm(seg, 2, 1)';
max_len = max(seg_len);
n_t = max(3, ceil(max_len / edge_sample_spacing) + 1);
n_t = min(n_t, max_edge_samples);
t = linspace(0, 1, n_t);

T = reshape(t, 1, n_t, 1);
P0r = reshape(P0, 3, 1, n_edges);
P1r = reshape(P1, 3, 1, n_edges);
Pts = P0r .* (1 - T) + P1r .* T;
Pts_flat = reshape(Pts, 3, []);

inb = points_in_bounds(Pts_flat, bounds);
clear = points_clear_of_terrain(Pts_flat, terrain_query, point_clearance);
risk = risk_query(Pts_flat);
risk = max(0, min(1, reshape(risk, 1, [])));
finite_risk = isfinite(risk);

valid = inb & clear & finite_risk;
valid = reshape(valid, n_t, n_edges);
risk = reshape(risk, n_t, n_edges);

ok_row = all(valid, 1);
ok = ok_row(:);
if any(ok)
    avg = mean(risk(:, ok_row), 1);
    avg_risk(ok) = avg(:);
    edge_w(ok) = seg_len(ok) .* (1 + stealth_weight .* avg_risk(ok));
end
end

function D = pairwise_dist(A, B)
AA = sum(A .^ 2, 1)';
BB = sum(B .^ 2, 1);
D2 = max(0.0, AA + BB - 2.0 .* (A' * B));
D = sqrt(D2);
end

function r = dynamic_connection_radius(bounds, n, r_min, r_max, eta)
vol = max((bounds(2) - bounds(1)) * (bounds(4) - bounds(3)) * ...
          (bounds(6) - bounds(5)), 1e-9);
zeta3 = 4.0 * pi / 3.0;
gamma = 2.0 * (1.0 + 1.0 / 3.0)^(1.0 / 3.0) * (vol / zeta3)^(1.0 / 3.0);
n_eff = max(2, n);
r_asym = eta * gamma * (log(n_eff) / n_eff)^(1.0 / 3.0);
r = min(r_max, max(r_min, r_asym));
end

function path = extract_path(nodes, parents, goal_parent, x_goal, max_nodes)
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
path = [nodes(:, chain), x_goal(:)];
d = vecnorm(diff(path, 1, 2), 2, 1);
path = path(:, [true, d > 1e-9]);
end

function mr = mean_path_risk(path, risk_query, sample_spacing, max_samples)
if isempty(path) || size(path, 2) < 2
    mr = NaN;
    return;
end
risks = [];
for k = 1:size(path, 2)-1
    p0 = path(:, k);
    p1 = path(:, k+1);
    len = norm(p1 - p0);
    n = min(max(2, ceil(len / sample_spacing) + 1), max_samples);
    t = linspace(0, 1, n);
    P = p0 .* (1 - t) + p1 .* t;
    risks = [risks, risk_query(P)]; %#ok<AGROW>
end
risks = risks(isfinite(risks));
if isempty(risks)
    mr = NaN;
else
    mr = mean(max(0, min(1, risks)));
end
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
