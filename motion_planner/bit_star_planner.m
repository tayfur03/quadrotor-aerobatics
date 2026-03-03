function [path, info] = bit_star_planner(start, goal, corridor, params)
%BIT_STAR_PLANNER Corridor-guided batch informed tree planner.
%
% Inputs
%   start     [3x1] start point [N;E;Alt] (or voxel coords if no world axes)
%   goal      [3x1] goal point  [N;E;Alt]
%   corridor  struct from compute_stealth_corridor
%   params    optional fields:
%     .max_batches             (default 80)
%     .batch_size              (default 240)
%     .max_nodes               (default 20000)
%     .eps_global              epsilon global-sampling probability (default 0.08)
%     .stealth_weight          risk penalty multiplier (default 2.0)
%                              set to inf for strict no-radar-violation mode
%     .edge_check_samples      edge safety samples (default 7)
%     .post_solution_batches   improvement batches after first hit (default 12)
%     .goal_radius             goal connection radius (default auto)
%     .min_connection_radius   minimum dynamic radius (default auto)
%     .max_connection_radius   maximum dynamic radius (default auto)
%     .eta_radius              asymptotic radius scaling (default 1.8)
%     .max_time                max planning time [s] (default inf)
%
% Outputs
%   path [3xM] planned path [N;E;Alt] (empty if failure)
%   info struct with timing and convergence metrics
%
% Notes
% - Sampling is restricted to the union of corridor spheres for
%   (1-eps_global) of samples. The remaining eps_global samples are global,
%   preserving probabilistic completeness.
% - Deterministic geometric corridor guidance is explainable and auditable,
%   unlike black-box learned initializers, which is useful for safety cases.

if nargin < 4
    params = struct();
end

validateattributes(start, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'start', 1);
validateattributes(goal, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'goal', 2);
start = start(:);
goal = goal(:);

if ~isstruct(corridor) || ~isfield(corridor, 'skeleton') || ~isfield(corridor, 'radii')
    error('bit_star_planner:InvalidCorridor', ...
          'corridor must include .skeleton and .radii fields.');
end

centers = corridor.skeleton;
radii = corridor.radii(:)';
if size(centers, 1) ~= 3
    error('bit_star_planner:InvalidCorridor', 'corridor.skeleton must be [3xM].');
end
if numel(radii) ~= size(centers, 2)
    error('bit_star_planner:InvalidCorridor', 'corridor.radii must match skeleton point count.');
end
if isempty(centers)
    error('bit_star_planner:InvalidCorridor', 'corridor.skeleton is empty.');
end

max_batches = round(get_opt(params, 'max_batches', 80));
batch_size = round(get_opt(params, 'batch_size', 240));
max_nodes = round(get_opt(params, 'max_nodes', 20000));
eps_global = get_opt(params, 'eps_global', 0.08);
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

bounds = get_bounds_with_endpoints(corridor, start, goal);

default_goal_radius = max(5.0, 0.35 * median(max(radii, 1e-3)));
goal_radius = get_opt(params, 'goal_radius', default_goal_radius);
min_r = get_opt(params, 'min_connection_radius', max(5.0, 0.20 * median(max(radii, 1e-3))));
max_r = get_opt(params, 'max_connection_radius', max(max(radii), min_r * 2.0));
eta_radius = get_opt(params, 'eta_radius', 1.8);

[risk_query, vis_threshold] = make_risk_query(corridor);
if isinf(stealth_weight)
    risk_limit = strict_zero_risk_tol;
else
    risk_limit = vis_threshold;
end

if ~isfinite(risk_query(start)) || risk_query(start) > vis_threshold
    warning('bit_star_planner:UnsafeStart', ...
        'Start appears outside safe corridor mask. Planner will continue.');
end
if ~isfinite(risk_query(goal)) || risk_query(goal) > vis_threshold
    warning('bit_star_planner:UnsafeGoal', ...
        'Goal appears outside safe corridor mask. Planner will continue.');
end

nodes = nan(3, max_nodes);
parents = zeros(1, max_nodes, 'uint32');
g_cost = inf(1, max_nodes);

nodes(:, 1) = start;
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
    n_global = round(eps_global * batch_size);
    n_corr = batch_size - n_global;

    S = zeros(3, batch_size);
    if n_corr > 0
        S(:, 1:n_corr) = sample_from_corridor_union(centers, radii, n_corr);
    end
    if n_global > 0
        S(:, n_corr+1:batch_size) = sample_uniform_bounds(bounds, n_global);
    end
    samples_generated = samples_generated + batch_size;

    in_bounds = points_in_bounds(S, bounds);
    risk_S = risk_query(S);
    safe_samples = in_bounds & isfinite(risk_S) & (risk_S <= risk_limit);
    S = S(:, safe_samples);

    if isempty(S)
        continue;
    end

    samples_accepted = samples_accepted + size(S, 2);

    n_tree = node_count;
    n_samples = size(S, 2);
    if n_samples == 0
        continue;
    end

    r_conn = dynamic_connection_radius(bounds, n_tree + n_samples, min_r, max_r, eta_radius);
    D = pairwise_dist(nodes(:, 1:n_tree), S);  % [n_tree x n_samples]
    edge_mask = D <= r_conn;

    [~, nearest_idx] = min(D, [], 1);
    nearest_lin = sub2ind(size(edge_mask), nearest_idx, 1:n_samples);
    edge_mask(nearest_lin) = true;

    [parent_local, sample_local] = find(edge_mask);
    edge_count = numel(parent_local);
    if edge_count == 0
        continue;
    end

    P = nodes(:, parent_local);
    C = S(:, sample_local);
    M = 0.5 * (P + C);

    risk_mid = risk_query(M);
    edge_len = D(sub2ind(size(D), parent_local, sample_local));
    edge_len = reshape(edge_len, [], 1);
    risk_mid = reshape(risk_mid, [], 1);
    g_parent = reshape(g_cost(parent_local), [], 1);
    if isinf(stealth_weight)
        edge_penalty = inf(size(edge_len));
        edge_penalty(risk_mid <= risk_limit) = 0;
        edge_w = edge_len + edge_penalty;
    else
        edge_w = edge_len + stealth_weight .* edge_len .* max(risk_mid, 0);
    end
    total_w = g_parent + edge_w;

    parent_costs = inf(n_tree, n_samples);
    lin_ps = sub2ind(size(parent_costs), parent_local(:), sample_local(:));
    [lin_ps_u, iu_ps] = unique(lin_ps, 'stable');
    parent_costs(lin_ps_u) = total_w(iu_ps);
    [best_sample_cost, best_parent_local] = min(parent_costs, [], 1);
    candidate_mask = isfinite(best_sample_cost);

    if ~any(candidate_mask)
        continue;
    end

    cand_samples = S(:, candidate_mask);
    cand_parents = best_parent_local(candidate_mask);
    cand_costs = best_sample_cost(candidate_mask);
    cand_parent_pts = nodes(:, cand_parents);

    ok_edges = edge_is_safe_batch(cand_parent_pts, cand_samples, edge_check_samples, bounds, risk_query, risk_limit);
    if ~any(ok_edges)
        continue;
    end

    cand_samples = cand_samples(:, ok_edges);
    cand_parents = cand_parents(ok_edges);
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
    parents(new_indices) = uint32(cand_parents(1:n_add));
    g_cost(new_indices) = cand_costs(1:n_add);
    node_count = old_count + n_add;

    % Batch rewiring (vectorized candidate generation + filtered edge checks).
    if old_count >= 2 && n_add >= 1
        existing_idx = 2:old_count;  % do not rewire root
        new_idx = old_count + 1:node_count;

        D_rw = pairwise_dist(nodes(:, new_idx), nodes(:, existing_idx));  % [n_new x n_old]
        r_rewire = 1.10 * r_conn;
        rw_mask = D_rw <= r_rewire;

        if any(rw_mask(:))
            [new_local, ex_local] = find(rw_mask);
            new_abs = new_idx(new_local);
            ex_abs = existing_idx(ex_local);

            P_rw = nodes(:, new_abs);
            C_rw = nodes(:, ex_abs);
            M_rw = 0.5 * (P_rw + C_rw);

            rw_risk = risk_query(M_rw);
            rw_len = D_rw(sub2ind(size(D_rw), new_local, ex_local));
            rw_len = reshape(rw_len, [], 1);
            rw_risk = reshape(rw_risk, [], 1);
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
                ok_rewire = edge_is_safe_batch(P_best, C_best, edge_check_samples, bounds, risk_query, risk_limit);

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

    % Goal connection.
    D_goal = vecnorm(nodes(:, 1:node_count) - goal, 2, 1);
    goal_candidates = find(D_goal <= max(goal_radius, r_conn));
    if ~isempty(goal_candidates)
        P_goal = nodes(:, goal_candidates);
        G_goal = repmat(goal, 1, numel(goal_candidates));
        M_goal = 0.5 * (P_goal + G_goal);
        goal_risk = risk_query(M_goal);
        goal_len = D_goal(goal_candidates);
        goal_len = reshape(goal_len, [], 1);
        goal_risk = reshape(goal_risk, [], 1);
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
            G_try = repmat(goal, 1, numel(cand_nodes));
            ok_goal = edge_is_safe_batch(P_try, G_try, edge_check_samples, bounds, risk_query, risk_limit);

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
success = goal_parent > 0;

if success
    path = extract_path(nodes, parents, goal_parent, start, goal, max_nodes);
else
    path = [];
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
end

function bounds = get_bounds_with_endpoints(corridor, start, goal)
if isfield(corridor, 'bounds_world') && ~isempty(corridor.bounds_world)
    bounds = corridor.bounds_world(:)';
else
    c = corridor.skeleton;
    r = corridor.radii;
    bounds = [min(c(1, :) - r), max(c(1, :) + r), ...
              min(c(2, :) - r), max(c(2, :) + r), ...
              min(c(3, :) - r), max(c(3, :) + r)];
end

bounds(1) = min([bounds(1), start(1), goal(1)]);
bounds(2) = max([bounds(2), start(1), goal(1)]);
bounds(3) = min([bounds(3), start(2), goal(2)]);
bounds(4) = max([bounds(4), start(2), goal(2)]);
bounds(5) = min([bounds(5), start(3), goal(3)]);
bounds(6) = max([bounds(6), start(3), goal(3)]);
end

function [risk_query, vis_threshold] = make_risk_query(corridor)
vis_threshold = corridor.visibility_threshold;

if isfield(corridor, 'grid') && isfield(corridor.grid, 'has_world_axes') && corridor.grid.has_world_axes
    E_vec = corridor.grid.E_vec;
    N_vec = corridor.grid.N_vec;
    alt_vec = corridor.grid.alt_vec;
    V = corridor.V;
    F = griddedInterpolant({E_vec, N_vec, alt_vec}, V, 'linear', 'none');
    risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);
else
    V = corridor.V;
    risk_query = @(P) reshape(interpn(V, P(1, :), P(2, :), P(3, :), 'linear', NaN), 1, []);
end
end

function S = sample_from_corridor_union(centers, radii, n)
if n <= 0
    S = zeros(3, 0);
    return;
end

weights = max(radii, 1e-6) .^ 3;
cdf = cumsum(weights);
cdf = cdf ./ cdf(end);

u = rand(1, n);
sphere_idx = sum(u(:) > cdf, 2) + 1;
sphere_idx = sphere_idx(:)';
sphere_idx = min(max(sphere_idx, 1), numel(radii));

C = centers(:, sphere_idx);
R = radii(sphere_idx);

dirs = randn(3, n);
dir_norm = vecnorm(dirs, 2, 1);
dir_norm = max(dir_norm, 1e-12);
dirs = dirs ./ dir_norm;

rho = R .* (rand(1, n) .^ (1.0 / 3.0));
S = C + dirs .* rho;
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
% A: [3 x nA], B: [3 x nB] -> D: [nA x nB]
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

function ok = edge_is_safe_batch(P0, P1, n_samples, bounds, risk_query, risk_limit)
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
Pts = P0r .* (1 - T) + P1r .* T;  % [3 x n_t x n_edges]

Pts_flat = reshape(Pts, 3, []);
inb = points_in_bounds(Pts_flat, bounds);
risk = risk_query(Pts_flat);
safe = inb & isfinite(risk) & (risk <= risk_limit);
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

chain = chain(1:count);
chain = fliplr(chain);

n_mid = numel(chain);
path = zeros(3, n_mid + 2);
path(:, 1) = start;
path(:, 2:n_mid+1) = nodes(:, chain);
path(:, n_mid+2) = goal;

% Remove consecutive duplicates.
d = vecnorm(diff(path, 1, 2), 2, 1);
keep = [true, d > 1e-9];
path = path(:, keep);
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
