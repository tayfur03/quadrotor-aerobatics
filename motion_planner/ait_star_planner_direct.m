function [path, cost, info] = ait_star_planner_direct(x_start, x_goal, env, params)
%AIT_STAR_PLANNER_DIRECT Direct AIT*-style planner for radar-aware environments.
%
% Inputs
%   x_start   [3x1] start point [N;E;Alt]
%   x_goal    [3x1] goal point  [N;E;Alt]
%   env       direct planning environment, same contract as FSM direct planner
%   params    planner options
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
terrain_query = ctx.terrain_query;
risk_query = ctx.risk_query;
bounds = ctx.bounds;
point_clearance = ctx.point_clearance;
max_clearance = ctx.max_clearance;
risk_limit = ctx.risk_limit;
altitude_cost_weight = ctx.altitude_cost_weight;
preferred_agl = ctx.preferred_agl;

max_batches = round(get_opt(params, 'max_batches', 80));
batch_size = round(get_opt(params, 'batch_size', 240));
max_nodes = round(get_opt(params, 'max_nodes', 20000));
eps_global = get_opt(params, 'eps_global', 0.12);
stealth_weight = get_opt(params, 'stealth_weight', 2.0);
edge_check_samples = round(get_opt(params, 'edge_check_samples', 7));
post_solution_batches = round(get_opt(params, 'post_solution_batches', 12));
max_time = get_opt(params, 'max_time', inf);
reverse_k_neighbors = round(get_opt(params, 'reverse_k_neighbors', 24));
forward_k_neighbors = round(get_opt(params, 'forward_k_neighbors', 16));
reverse_debug = logical(get_opt(params, 'reverse_debug', false));
goal_sample_fraction = get_opt(params, 'goal_sample_fraction', 0.20);
line_sample_fraction = get_opt(params, 'line_sample_fraction', 0.25);
goal_probe_count = round(get_opt(params, 'goal_probe_count', 10));

profile_ait = logical(get_opt(params, 'profile_ait', false));
neighbor_precompute_mode = get_opt(params, 'neighbor_precompute_mode', 'batch');
if isstring(neighbor_precompute_mode)
    neighbor_precompute_mode = char(neighbor_precompute_mode);
end
if ~ischar(neighbor_precompute_mode)
    neighbor_precompute_mode = 'batch';
end
neighbor_precompute_mode = lower(neighbor_precompute_mode);
use_parallel_reverse = logical(get_opt(params, 'use_parallel_reverse', false));
use_parallel_edge_checks = logical(get_opt(params, 'use_parallel_edge_checks', false));
parallel_threshold_neighbors = round(get_opt(params, 'parallel_threshold_neighbors', 128));
parallel_threshold_edges = round(get_opt(params, 'parallel_threshold_edges', 64));
reverse_batch_repair_limit = round(get_opt(params, 'reverse_batch_repair_limit', 256));
forward_expand_chunk = round(get_opt(params, 'forward_expand_chunk', 32));
kd_rebuild_interval = round(get_opt(params, 'kd_rebuild_interval', 1)); %#ok<NASGU>
graph_rebuild_interval = round(get_opt(params, 'graph_rebuild_interval', 1));
reverse_full_rebuild_interval = round(get_opt(params, 'reverse_full_rebuild_interval', graph_rebuild_interval));

max_batches = max(1, max_batches);
batch_size = max(16, batch_size);
max_nodes = max(batch_size + 2, max_nodes);
eps_global = min(max(eps_global, 0.0), 0.5);
edge_check_samples = max(3, edge_check_samples);
post_solution_batches = max(0, post_solution_batches);
reverse_k_neighbors = max(4, reverse_k_neighbors);
forward_k_neighbors = max(4, forward_k_neighbors);
goal_sample_fraction = min(max(goal_sample_fraction, 0.0), 0.7);
line_sample_fraction = min(max(line_sample_fraction, 0.0), 0.7);
goal_probe_count = max(0, goal_probe_count);
parallel_threshold_neighbors = max(2, parallel_threshold_neighbors);
parallel_threshold_edges = max(2, parallel_threshold_edges);
reverse_batch_repair_limit = max(8, reverse_batch_repair_limit);
forward_expand_chunk = max(1, forward_expand_chunk);
graph_rebuild_interval = max(1, graph_rebuild_interval);
reverse_full_rebuild_interval = max(1, reverse_full_rebuild_interval);
if ~strcmp(neighbor_precompute_mode, 'batch')
    warning('AIT*: Unsupported neighbor_precompute_mode="%s". Falling back to "batch".', neighbor_precompute_mode);
    neighbor_precompute_mode = 'batch';
end
if use_parallel_reverse && ~can_use_parfor()
    warning('AIT*: parfor unavailable. Falling back to serial reverse precompute.');
    use_parallel_reverse = false;
end
if use_parallel_edge_checks && ~can_use_parfor()
    warning('AIT*: parfor unavailable. Falling back to serial edge checks.');
    use_parallel_edge_checks = false;
end

default_goal_radius = max(5.0, 0.02 * norm(bounds([2 4 6]) - bounds([1 3 5])));
goal_radius = get_opt(params, 'goal_radius', default_goal_radius);
min_r = get_opt(params, 'min_connection_radius', max(5.0, 0.012 * norm(bounds([2 4 6]) - bounds([1 3 5]))));
max_r = get_opt(params, 'max_connection_radius', max(40, 0.08 * norm(bounds([2 4 6]) - bounds([1 3 5]))));
eta_radius = get_opt(params, 'eta_radius', 1.8);

nodes = nan(3, max_nodes);
nodes(:, 1) = x_start;
nodes(:, 2) = x_goal;
node_count = 2;

in_tree_fwd = false(1, max_nodes);
parent_fwd = zeros(1, max_nodes, 'uint32');
g_fwd = inf(1, max_nodes);
in_tree_fwd(1) = true;
g_fwd(1) = 0.0;

rev_state = init_reverse_state(max_nodes);
invalid_adj = spalloc(max_nodes, max_nodes, max_nodes * 6);
fwd_state = init_forward_edge_state(max_nodes);

best_goal_cost = Inf;
first_solution_batch = NaN;
first_solution_time = NaN;
samples_generated = 0;
samples_accepted = 0;
batches_run = 0;
stop_reason = 'max_batches_reached';

reverse_rebuild_time_s = 0;
reverse_repair_time_s = 0;
forward_time_s = 0;
forward_queue_build_time_s = 0;
goal_probe_time_s = 0;
neighbor_precompute_time_s = 0;
edge_validation_time_s = 0;
reverse_repairs = 0;
invalid_edge_count = 0;
reverse_queue_pops = 0;
forward_edge_checks = 0;
adaptive_heuristic_ready_count = 0;
valid_edge_cache_hits = 0;
forward_stale_pops = 0;
goal_probe_checks = 0;
goal_probe_successes = 0;
parallel_used_reverse = false;
parallel_used_edge_checks = false;

t_start = tic;
last_r_graph = max_r;
graph_state = [];
batch_profile = repmat(empty_batch_profile(), 1, max_batches);

for batch_idx = 1:max_batches
    if toc(t_start) > max_time
        stop_reason = 'max_time';
        break;
    end

    batches_run = batch_idx;
    batch_edge_checks_start = forward_edge_checks;
    batch_edge_eval_start = edge_validation_time_s;
    batch_accepted_start = samples_accepted;
    graph_elapsed = 0;
    rev_elapsed = 0;
    queue_elapsed = 0;
    probe_elapsed = 0;
    forward_elapsed = 0;
    S_all = sample_mixed_bounds(bounds, x_start, x_goal, batch_size, goal_sample_fraction, line_sample_fraction);
    samples_generated = samples_generated + batch_size;

    in_bounds = points_in_bounds(S_all, bounds);
    risk_all = risk_query(S_all);
    safe_mask = in_bounds & isfinite(risk_all) & (risk_all <= risk_limit);
    safe_mask = safe_mask & ctx.points_clear_of_terrain(S_all);

    if isfinite(best_goal_cost)
        eucl_informed = (vecnorm(S_all - x_start, 2, 1) + vecnorm(S_all - x_goal, 2, 1)) < best_goal_cost;
        n_global = round(eps_global * batch_size);
        global_mask = false(1, batch_size);
        if n_global > 0
            global_mask(randperm(batch_size, min(n_global, batch_size))) = true;
        end
        keep_mask = safe_mask & (global_mask | eucl_informed);
    else
        keep_mask = safe_mask;
    end

    S = S_all(:, keep_mask);
    new_idx = zeros(1, 0);
    if ~isempty(S)
        add_cap = max_nodes - node_count;
        n_add = min(add_cap, size(S, 2));
        if n_add <= 0
            stop_reason = 'max_nodes';
            break;
        end
        add_idx = node_count + (1:n_add);
        nodes(:, add_idx) = S(:, 1:n_add);
        node_count = node_count + n_add;
        samples_accepted = samples_accepted + n_add;
        new_idx = add_idx;
    end

    last_r_graph = dynamic_connection_radius(bounds, node_count, min_r, max_r, eta_radius);

    t_graph = tic;
    force_graph_rebuild = isempty(graph_state) || batch_idx == 1 || ...
        graph_rebuild_interval <= 1 || mod(batch_idx - 1, graph_rebuild_interval) == 0 || ...
        isempty(new_idx);
    if force_graph_rebuild
        [graph_state, graph_parallel_used] = build_graph_state( ...
            nodes, node_count, last_r_graph, max(goal_radius, last_r_graph), ...
            reverse_k_neighbors, forward_k_neighbors, 1, 2, neighbor_precompute_mode, ...
            use_parallel_reverse, parallel_threshold_neighbors);
    else
        [graph_state, graph_parallel_used] = append_graph_state( ...
            graph_state, nodes, node_count, new_idx, last_r_graph, max(goal_radius, last_r_graph), ...
            reverse_k_neighbors, forward_k_neighbors, 2, use_parallel_reverse, parallel_threshold_neighbors);
    end
    graph_elapsed = toc(t_graph);
    neighbor_precompute_time_s = neighbor_precompute_time_s + graph_elapsed;
    parallel_used_reverse = parallel_used_reverse || graph_parallel_used;

    t_rev = tic;
    force_reverse_rebuild = force_graph_rebuild || reverse_full_rebuild_interval <= 1 || ...
        mod(batch_idx - 1, reverse_full_rebuild_interval) == 0;
    if force_reverse_rebuild
        [rev_state, pops_batch] = rebuild_reverse_search(rev_state, graph_state, invalid_adj, 1, 2);
    else
        [rev_state, pops_batch] = repair_reverse_after_new_nodes( ...
            rev_state, graph_state, invalid_adj, new_idx, 1, 2, reverse_batch_repair_limit);
    end
    rev_elapsed = toc(t_rev);
    reverse_rebuild_time_s = reverse_rebuild_time_s + rev_elapsed;
    reverse_queue_pops = reverse_queue_pops + pops_batch;
    if reverse_debug || profile_ait
        adaptive_heuristic_ready_count = sum(isfinite(rev_state.g(1:node_count)));
    end
    if reverse_debug
        fprintf(['AIT* batch %d: nodes=%d r=%.1f graph_t=%.2fs reverse_t=%.2fs pops=%d ', ...
                 'avg_deg=%.1f max_deg=%d ready=%d\n'], ...
            batch_idx, node_count, last_r_graph, graph_elapsed, rev_elapsed, pops_batch, ...
            graph_state.neighbor_avg_degree, graph_state.neighbor_max_degree, adaptive_heuristic_ready_count);
    end

    fwd_state = reset_forward_open_state(fwd_state);
    t_queue = tic;
    edge_heap = init_edge_heap(max(4096, node_count * 8));
    [edge_heap, fwd_state] = rebuild_forward_queue( ...
        edge_heap, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost);
    queue_elapsed = toc(t_queue);
    forward_queue_build_time_s = forward_queue_build_time_s + queue_elapsed;

    t_probe = tic;
    [best_goal_cost, in_tree_fwd, parent_fwd, g_fwd, invalid_adj, fwd_state, ...
        probe_checks_batch, probe_success_batch, edge_heap, probe_eval_time, probe_parallel_used, ...
        queue_rebuild_time_probe] = attempt_goal_probes( ...
        best_goal_cost, in_tree_fwd, parent_fwd, g_fwd, nodes, node_count, ...
        rev_state, graph_state, invalid_adj, fwd_state, x_goal, edge_check_samples, ...
        bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance, stealth_weight, altitude_cost_weight, preferred_agl, ...
        goal_probe_count, edge_heap, use_parallel_edge_checks, parallel_threshold_edges);
    probe_elapsed = toc(t_probe);
    goal_probe_time_s = goal_probe_time_s + probe_elapsed;
    edge_validation_time_s = edge_validation_time_s + probe_eval_time;
    parallel_used_edge_checks = parallel_used_edge_checks || probe_parallel_used;
    forward_queue_build_time_s = forward_queue_build_time_s + queue_rebuild_time_probe;
    goal_probe_checks = goal_probe_checks + probe_checks_batch;
    goal_probe_successes = goal_probe_successes + probe_success_batch;
    if probe_success_batch > 0 && isnan(first_solution_batch) && isfinite(best_goal_cost)
        first_solution_batch = batch_idx;
        first_solution_time = toc(t_start);
    end

    t_fwd = tic;
    while edge_heap.size > 0
        if toc(t_start) > max_time
            stop_reason = 'max_time';
            break;
        end

        [cand, edge_heap, fwd_state, stale_inc, stop_by_bound] = collect_forward_chunk( ...
            edge_heap, fwd_state, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, ...
            best_goal_cost, forward_expand_chunk);
        forward_stale_pops = forward_stale_pops + stale_inc;
        if stop_by_bound
            break;
        end
        if cand.count == 0
            continue;
        end

        cached_mask = get_edge_flags(fwd_state.valid_adj, cand.u, cand.v);
        valid_edge_cache_hits = valid_edge_cache_hits + nnz(cached_mask);
        unknown_mask = ~cached_mask;
        ok_mask = cached_mask;
        eval_time = 0;
        used_parallel_batch = false;

        if any(unknown_mask)
            forward_edge_checks = forward_edge_checks + nnz(unknown_mask);
            t_eval = tic;
            [ok_unknown, used_parallel_batch] = evaluate_edge_safety_chunk( ...
                nodes(:, cand.u(unknown_mask)), nodes(:, cand.v(unknown_mask)), edge_check_samples, ...
                bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance, ...
                use_parallel_edge_checks, parallel_threshold_edges);
            eval_time = toc(t_eval);
            edge_validation_time_s = edge_validation_time_s + eval_time;
            parallel_used_edge_checks = parallel_used_edge_checks || used_parallel_batch;
            ok_mask(unknown_mask) = ok_unknown;
        end

        if any(unknown_mask)
            fwd_state.valid_adj = set_edge_flags(fwd_state.valid_adj, cand.u(unknown_mask & ok_mask), cand.v(unknown_mask & ok_mask), true);
            fwd_state.valid_adj = set_edge_flags(fwd_state.valid_adj, cand.u(unknown_mask & ~ok_mask), cand.v(unknown_mask & ~ok_mask), false);
        end

        invalid_mask = ~ok_mask;
        if any(invalid_mask)
            bad_u = cand.u(invalid_mask);
            bad_v = cand.v(invalid_mask);
            invalid_adj = set_edge_flags(invalid_adj, bad_u, bad_v, true);
            invalid_edge_count = invalid_edge_count + numel(bad_u);
        end

        valid_mask = ok_mask;
        if any(valid_mask)
            edge_costs = compute_true_edge_cost_batch( ...
                nodes(:, cand.u(valid_mask)), nodes(:, cand.v(valid_mask)), risk_query, risk_limit, stealth_weight, terrain_query, altitude_cost_weight, preferred_agl);
            good_cost_mask = isfinite(edge_costs);
            if any(good_cost_mask)
                valid_u = cand.u(valid_mask);
                valid_v = cand.v(valid_mask);
                tentative = g_fwd(valid_u(good_cost_mask)) + edge_costs(good_cost_mask);
                proc_u = valid_u(good_cost_mask);
                proc_v = valid_v(good_cost_mask);
                [~, order] = sort(tentative, 'ascend');
                proc_u = proc_u(order);
                proc_v = proc_v(order);
                tentative = tentative(order);
                for i_edge = 1:numel(proc_u)
                    u = proc_u(i_edge);
                    v = proc_v(i_edge);
                    if tentative(i_edge) + 1e-9 < g_fwd(v) || ~in_tree_fwd(v)
                        in_tree_fwd(v) = true;
                        parent_fwd(v) = uint32(u);
                        g_fwd(v) = tentative(i_edge);
                        [edge_heap, fwd_state] = enqueue_forward_neighbors( ...
                            edge_heap, v, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost);
                        if v == 2
                            best_goal_cost = tentative(i_edge);
                            if isnan(first_solution_batch)
                                first_solution_batch = batch_idx;
                                first_solution_time = toc(t_start);
                            end
                        end
                    end
                end
            end
        end

        if any(invalid_mask)
            t_repair = tic;
            [rev_state, pops_repair] = repair_reverse_after_invalid_edges( ...
                rev_state, graph_state, invalid_adj, cand.u(invalid_mask), cand.v(invalid_mask), ...
                1, 2, reverse_batch_repair_limit);
            repair_elapsed = toc(t_repair);
            reverse_repair_time_s = reverse_repair_time_s + repair_elapsed;
            reverse_queue_pops = reverse_queue_pops + pops_repair;
            reverse_repairs = reverse_repairs + 1;
            if reverse_debug || profile_ait
                adaptive_heuristic_ready_count = sum(isfinite(rev_state.g(1:node_count)));
            end

            fwd_state = reset_forward_open_state(fwd_state);
            t_queue = tic;
            edge_heap = init_edge_heap(max(4096, node_count * 8));
            [edge_heap, fwd_state] = rebuild_forward_queue( ...
                edge_heap, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost);
            forward_queue_build_time_s = forward_queue_build_time_s + toc(t_queue);
        end
    end
    forward_elapsed = toc(t_fwd);
    forward_time_s = forward_time_s + forward_elapsed;
    batch_profile(batch_idx) = make_batch_profile(batch_idx, node_count, graph_elapsed, rev_elapsed, ...
        queue_elapsed, probe_elapsed, forward_elapsed, edge_validation_time_s - batch_edge_eval_start, ...
        samples_accepted - batch_accepted_start, forward_edge_checks - batch_edge_checks_start, ...
        force_graph_rebuild, force_reverse_rebuild);

    if strcmp(stop_reason, 'max_time')
        break;
    end
    if in_tree_fwd(2) && batch_idx >= first_solution_batch + post_solution_batches
        stop_reason = 'post_solution_budget';
        break;
    end
    if node_count >= max_nodes
        stop_reason = 'max_nodes';
        break;
    end
end

planning_time = toc(t_start);
success = in_tree_fwd(2) && isfinite(best_goal_cost);
if ~success
    best_goal_cost = Inf;
end

if success
    path = extract_tree_path(nodes, parent_fwd, 2, max_nodes);
else
    path = [];
end

reverse_time_s = reverse_rebuild_time_s + reverse_repair_time_s + neighbor_precompute_time_s;
if ~reverse_debug && ~profile_ait
    adaptive_heuristic_ready_count = sum(isfinite(rev_state.g(1:node_count)));
end

info = struct();
info.success = success;
info.iterations = batches_run;
info.tree_size = sum(in_tree_fwd(1:node_count));
info.planning_time = planning_time;
info.first_solution_batch = first_solution_batch;
info.first_solution_time = first_solution_time;
info.path_cost = best_goal_cost;
info.samples_generated = samples_generated;
info.samples_accepted = samples_accepted;
info.connection_radius_last = last_r_graph;
info.stop_reason = stop_reason;
info.risk_limit = risk_limit;
info.ait_used = true;
info.reverse_time_s = reverse_time_s;
info.forward_time_s = forward_time_s;
info.reverse_rebuild_time_s = reverse_rebuild_time_s;
info.reverse_repair_time_s = reverse_repair_time_s;
info.forward_queue_build_time_s = forward_queue_build_time_s;
info.goal_probe_time_s = goal_probe_time_s;
info.neighbor_precompute_time_s = neighbor_precompute_time_s;
info.edge_validation_time_s = edge_validation_time_s;
info.parallel_used_reverse = parallel_used_reverse;
info.parallel_used_edge_checks = parallel_used_edge_checks;
info.reverse_repairs = reverse_repairs;
info.invalid_edge_count = invalid_edge_count;
info.reverse_queue_pops = reverse_queue_pops;
info.forward_edge_checks = forward_edge_checks;
info.adaptive_heuristic_ready_count = adaptive_heuristic_ready_count;
info.valid_edge_cache_hits = valid_edge_cache_hits;
info.forward_stale_pops = forward_stale_pops;
info.goal_probe_checks = goal_probe_checks;
info.goal_probe_successes = goal_probe_successes;
info.profile_ait = profile_ait;
info.neighbor_precompute_mode = neighbor_precompute_mode;
info.env = ctx.env_meta;
info.sampled_nodes = nodes(:, 1:node_count);
info.in_tree_fwd = in_tree_fwd(1:node_count);
info.parent_fwd = parent_fwd(1:node_count);
info.g_fwd = g_fwd(1:node_count);
info.g_rev = rev_state.g(1:node_count);
info.batch_profile = batch_profile(1:max(batches_run, 0));
cost = best_goal_cost;
end

function p = empty_batch_profile()
p = struct( ...
    'batch', 0, ...
    'nodes', 0, ...
    'graph_t', 0, ...
    'reverse_t', 0, ...
    'queue_t', 0, ...
    'probe_t', 0, ...
    'forward_t', 0, ...
    'edge_eval_t', 0, ...
    'accepted_samples', 0, ...
    'edge_checks', 0, ...
    'full_graph_rebuild', false, ...
    'full_reverse_rebuild', false);
end

function p = make_batch_profile(batch_idx, node_count, graph_t, reverse_t, queue_t, probe_t, forward_t, edge_eval_t, accepted_samples, edge_checks, full_graph_rebuild, full_reverse_rebuild)
p = empty_batch_profile();
p.batch = batch_idx;
p.nodes = node_count;
p.graph_t = graph_t;
p.reverse_t = reverse_t;
p.queue_t = queue_t;
p.probe_t = probe_t;
p.forward_t = forward_t;
p.edge_eval_t = edge_eval_t;
p.accepted_samples = accepted_samples;
p.edge_checks = edge_checks;
p.full_graph_rebuild = full_graph_rebuild;
p.full_reverse_rebuild = full_reverse_rebuild;
end

function rev_state = init_reverse_state(max_nodes)
rev_state = struct();
rev_state.g = inf(1, max_nodes);
rev_state.rhs = inf(1, max_nodes);
rev_state.parent = zeros(1, max_nodes, 'uint32');
rev_state.children = cell(1, max_nodes);
rev_state.open_pending = false(1, max_nodes);
rev_state.open_key1 = inf(1, max_nodes);
rev_state.open_key2 = inf(1, max_nodes);
end

function fwd_state = init_forward_edge_state(max_nodes)
cap = max(max_nodes * 12, 1024);
fwd_state = struct();
fwd_state.open_pending = spalloc(max_nodes, max_nodes, cap);
fwd_state.open_key1 = spalloc(max_nodes, max_nodes, cap);
fwd_state.open_key2 = spalloc(max_nodes, max_nodes, cap);
fwd_state.valid_adj = spalloc(max_nodes, max_nodes, cap);
end

function fwd_state = reset_forward_open_state(fwd_state)
sz = size(fwd_state.open_pending);
fwd_state.open_pending = spalloc(sz(1), sz(2), nnz(fwd_state.open_pending));
fwd_state.open_key1 = spalloc(sz(1), sz(2), nnz(fwd_state.open_key1));
fwd_state.open_key2 = spalloc(sz(1), sz(2), nnz(fwd_state.open_key2));
end

function fwd_state = clear_forward_open_edge(fwd_state, u, v)
fwd_state.open_pending(u, v) = false;
fwd_state.open_key1(u, v) = 0;
fwd_state.open_key2(u, v) = 0;
end

function [graph_state, parallel_used] = build_graph_state(nodes, node_count, r_graph, goal_radius, reverse_cap, forward_cap, start_idx, goal_idx, neighbor_mode, use_parallel_reverse, parallel_threshold_neighbors) %#ok<INUSD>
graph_state = struct();
graph_state.nodes = nodes;
graph_state.node_count = node_count;
graph_state.r_graph = r_graph;
graph_state.goal_radius = goal_radius;
graph_state.reverse_k_neighbor_cap = reverse_cap;
graph_state.forward_k_neighbor_cap = forward_cap;
graph_state.node_to_start = vecnorm(nodes(:, 1:node_count) - nodes(:, start_idx), 2, 1);
graph_state.node_to_goal = vecnorm(nodes(:, 1:node_count) - nodes(:, goal_idx), 2, 1);
graph_state.reverse_neighbors = cell(1, node_count);
graph_state.reverse_lengths = cell(1, node_count);
graph_state.forward_neighbors = cell(1, node_count);
graph_state.forward_lengths = cell(1, node_count);
graph_state.neighbor_query_count = node_count;
graph_state.neighbor_total_degree = 0;
graph_state.neighbor_max_degree = 0;
graph_state.neighbor_avg_degree = 0;
parallel_used = false;

data = nodes(:, 1:node_count).';
raw_idx = cell(1, node_count);
raw_dist = cell(1, node_count);
if can_use_kdtree() && node_count >= 2
    try
        kd_tree = createns(data, 'NSMethod', 'kdtree');
        [raw_idx, raw_dist] = rangesearch(kd_tree, data, r_graph);
    catch
        [raw_idx, raw_dist] = brute_force_rangesearch(data, r_graph);
    end
else
    [raw_idx, raw_dist] = brute_force_rangesearch(data, r_graph);
end

do_parallel = use_parallel_reverse && node_count >= parallel_threshold_neighbors;
if do_parallel
    try
        rev_n = cell(1, node_count);
        rev_l = cell(1, node_count);
        fwd_n = cell(1, node_count);
        fwd_l = cell(1, node_count);
        deg = zeros(1, node_count);
        parfor idx = 1:node_count
            [rev_n{idx}, rev_l{idx}, fwd_n{idx}, fwd_l{idx}, deg(idx)] = cap_neighbor_lists( ...
                raw_idx{idx}, raw_dist{idx}, idx, goal_idx, graph_state.node_to_goal(idx), ...
                goal_radius, reverse_cap, forward_cap);
        end
        graph_state.reverse_neighbors = rev_n;
        graph_state.reverse_lengths = rev_l;
        graph_state.forward_neighbors = fwd_n;
        graph_state.forward_lengths = fwd_l;
        graph_state.neighbor_total_degree = sum(deg);
        graph_state.neighbor_max_degree = max(deg);
        graph_state.neighbor_avg_degree = mean(deg);
        parallel_used = true;
        return;
    catch
        parallel_used = false;
    end
end

deg = zeros(1, node_count);
for idx = 1:node_count
    [graph_state.reverse_neighbors{idx}, graph_state.reverse_lengths{idx}, ...
        graph_state.forward_neighbors{idx}, graph_state.forward_lengths{idx}, deg(idx)] = ...
        cap_neighbor_lists(raw_idx{idx}, raw_dist{idx}, idx, goal_idx, graph_state.node_to_goal(idx), ...
        goal_radius, reverse_cap, forward_cap);
end
graph_state.neighbor_total_degree = sum(deg);
graph_state.neighbor_max_degree = max(deg);
graph_state.neighbor_avg_degree = mean(deg);
end

function [graph_state, parallel_used] = append_graph_state(graph_state, nodes, node_count, new_idx, r_graph, goal_radius, reverse_cap, forward_cap, goal_idx, use_parallel_reverse, parallel_threshold_neighbors)
parallel_used = false;
old_count = graph_state.node_count;
new_idx = reshape(new_idx, 1, []);
new_idx = new_idx(new_idx >= 1 & new_idx <= node_count);
if isempty(new_idx) || node_count <= old_count
    graph_state.node_count = node_count;
    return;
end

graph_state.nodes = nodes;
graph_state.node_count = node_count;
graph_state.r_graph = r_graph;
graph_state.goal_radius = goal_radius;
graph_state.reverse_k_neighbor_cap = reverse_cap;
graph_state.forward_k_neighbor_cap = forward_cap;
graph_state.node_to_start = vecnorm(nodes(:, 1:node_count) - nodes(:, 1), 2, 1);
graph_state.node_to_goal = vecnorm(nodes(:, 1:node_count) - nodes(:, goal_idx), 2, 1);
graph_state.reverse_neighbors(old_count + 1:node_count) = {[]};
graph_state.reverse_lengths(old_count + 1:node_count) = {[]};
graph_state.forward_neighbors(old_count + 1:node_count) = {[]};
graph_state.forward_lengths(old_count + 1:node_count) = {[]};
graph_state.neighbor_query_count = graph_state.neighbor_query_count + numel(new_idx);

data = nodes(:, 1:node_count).';
do_parallel = use_parallel_reverse && numel(new_idx) >= parallel_threshold_neighbors;
if do_parallel
    try
        new_rev_n = cell(1, numel(new_idx));
        new_rev_l = cell(1, numel(new_idx));
        new_fwd_n = cell(1, numel(new_idx));
        new_fwd_l = cell(1, numel(new_idx));
        raw_idx_new = cell(1, numel(new_idx));
        raw_dist_new = cell(1, numel(new_idx));
        parfor ii = 1:numel(new_idx)
            idx = new_idx(ii);
            delta = data - data(idx, :);
            dist = sqrt(sum(delta .^ 2, 2));
            keep = dist <= r_graph;
            raw_idx_new{ii} = find(keep).';
            raw_dist_new{ii} = dist(keep).';
            [new_rev_n{ii}, new_rev_l{ii}, new_fwd_n{ii}, new_fwd_l{ii}, ~] = cap_neighbor_lists( ...
                raw_idx_new{ii}, raw_dist_new{ii}, idx, goal_idx, graph_state.node_to_goal(idx), ...
                goal_radius, reverse_cap, forward_cap);
        end
        for ii = 1:numel(new_idx)
            idx = new_idx(ii);
            graph_state.reverse_neighbors{idx} = new_rev_n{ii};
            graph_state.reverse_lengths{idx} = new_rev_l{ii};
            graph_state.forward_neighbors{idx} = new_fwd_n{ii};
            graph_state.forward_lengths{idx} = new_fwd_l{ii};
        end
        parallel_used = true;
    catch
        parallel_used = false;
    end
end

if ~parallel_used
    for ii = 1:numel(new_idx)
        idx = new_idx(ii);
        delta = data - data(idx, :);
        dist = sqrt(sum(delta .^ 2, 2));
        keep = dist <= r_graph;
        raw_idx = find(keep).';
        raw_dist = dist(keep).';
        [graph_state.reverse_neighbors{idx}, graph_state.reverse_lengths{idx}, ...
            graph_state.forward_neighbors{idx}, graph_state.forward_lengths{idx}, ~] = ...
            cap_neighbor_lists(raw_idx, raw_dist, idx, goal_idx, graph_state.node_to_goal(idx), ...
            goal_radius, reverse_cap, forward_cap);
    end
end

for ii = 1:numel(new_idx)
    idx = new_idx(ii);
    old_nodes = 1:(idx - 1);
    if isempty(old_nodes)
        continue;
    end
    delta_old = nodes(:, old_nodes) - nodes(:, idx);
    dist_old = sqrt(sum(delta_old .^ 2, 1));
    near_old = old_nodes(dist_old <= r_graph);
    dist_near = dist_old(dist_old <= r_graph);
    for jj = 1:numel(near_old)
        old = near_old(jj);
        d = dist_near(jj);
        [graph_state.reverse_neighbors{old}, graph_state.reverse_lengths{old}] = ...
            merge_capped_neighbor(graph_state.reverse_neighbors{old}, graph_state.reverse_lengths{old}, ...
            idx, d, old, goal_idx, graph_state.node_to_goal(old), goal_radius, reverse_cap, false);
        [graph_state.forward_neighbors{old}, graph_state.forward_lengths{old}] = ...
            merge_capped_neighbor(graph_state.forward_neighbors{old}, graph_state.forward_lengths{old}, ...
            idx, d, old, goal_idx, graph_state.node_to_goal(old), goal_radius, forward_cap, true);
    end
end

deg = cellfun(@numel, graph_state.forward_neighbors(1:node_count));
graph_state.neighbor_total_degree = sum(deg);
graph_state.neighbor_max_degree = max([0, deg]);
graph_state.neighbor_avg_degree = mean(deg);
end

function [idx_out, len_out] = merge_capped_neighbor(idx_in, len_in, new_idx, new_len, self_idx, goal_idx, goal_dist, goal_radius, cap, keep_goal)
idx_out = reshape(idx_in, 1, []);
len_out = reshape(len_in, 1, []);
if new_idx == self_idx || any(idx_out == new_idx)
    return;
end
idx_out = [idx_out, new_idx];
len_out = [len_out, new_len];
[len_out, order] = sort(len_out, 'ascend');
idx_out = idx_out(order);
if numel(idx_out) > cap
    keep = false(1, numel(idx_out));
    keep(1:min(cap, numel(idx_out))) = true;
    if keep_goal
        goal_loc = find(idx_out == goal_idx, 1);
        if isempty(goal_loc) && self_idx ~= goal_idx && goal_dist <= goal_radius
            idx_out = [idx_out, goal_idx];
            len_out = [len_out, goal_dist];
            [len_out, order] = sort(len_out, 'ascend');
            idx_out = idx_out(order);
            goal_loc = find(idx_out == goal_idx, 1);
            keep = false(1, numel(idx_out));
            keep(1:min(cap, numel(idx_out))) = true;
        end
        if ~isempty(goal_loc) && ~keep(goal_loc)
            last_kept = find(keep, 1, 'last');
            keep(last_kept) = false;
            keep(goal_loc) = true;
        end
    end
    idx_out = idx_out(keep);
    len_out = len_out(keep);
    [len_out, order] = sort(len_out, 'ascend');
    idx_out = idx_out(order);
end
end

function [raw_idx, raw_dist] = brute_force_rangesearch(data, radius)
node_count = size(data, 1);
raw_idx = cell(1, node_count);
raw_dist = cell(1, node_count);
for idx = 1:node_count
    delta = data - data(idx, :);
    dist = sqrt(sum(delta .^ 2, 2));
    keep = dist <= radius;
    raw_idx{idx} = find(keep).';
    raw_dist{idx} = dist(keep).';
end
end

function [rev_idx, rev_len, fwd_idx, fwd_len, degree] = cap_neighbor_lists(raw_idx, raw_dist, self_idx, goal_idx, goal_dist, goal_radius, reverse_cap, forward_cap)
if isempty(raw_idx)
    rev_idx = zeros(1, 0);
    rev_len = zeros(1, 0);
    fwd_idx = zeros(1, 0);
    fwd_len = zeros(1, 0);
    degree = 0;
    return;
end

raw_idx = reshape(raw_idx, 1, []);
raw_dist = reshape(raw_dist, 1, []);
keep = isfinite(raw_dist) & raw_idx ~= self_idx;
raw_idx = raw_idx(keep);
raw_dist = raw_dist(keep);
if isempty(raw_idx)
    rev_idx = zeros(1, 0);
    rev_len = zeros(1, 0);
    fwd_idx = zeros(1, 0);
    fwd_len = zeros(1, 0);
    degree = 0;
    return;
end

[raw_dist, order] = sort(raw_dist, 'ascend');
raw_idx = raw_idx(order);
degree = numel(raw_idx);

n_rev = min(degree, reverse_cap);
rev_idx = raw_idx(1:n_rev);
rev_len = raw_dist(1:n_rev);

n_fwd = min(degree, forward_cap);
fwd_idx = raw_idx(1:n_fwd);
fwd_len = raw_dist(1:n_fwd);

if self_idx ~= goal_idx && goal_dist <= goal_radius && ~any(fwd_idx == goal_idx)
    fwd_idx = [fwd_idx, goal_idx];
    fwd_len = [fwd_len, goal_dist];
    [fwd_len, order] = sort(fwd_len, 'ascend');
    fwd_idx = fwd_idx(order);
end
end

function [rev_state, pops] = rebuild_reverse_search(rev_state, graph_state, invalid_adj, start_idx, goal_idx)
node_count = graph_state.node_count;
rev_state.g(1:node_count) = inf;
rev_state.rhs(1:node_count) = inf;
rev_state.parent(1:node_count) = uint32(0);
rev_state.children(1:node_count) = {[]};
rev_state.open_pending(1:node_count) = false;
rev_state.open_key1(1:node_count) = inf;
rev_state.open_key2(1:node_count) = inf;
rev_state.rhs(goal_idx) = 0;

queue = init_vertex_heap(max(4096, node_count * 8));
[queue, rev_state] = queue_reverse_vertex(queue, rev_state, reverse_key(rev_state, graph_state, goal_idx), goal_idx);
[rev_state, ~, pops] = compute_reverse_shortest_path(rev_state, graph_state, invalid_adj, queue, start_idx, goal_idx, inf);
rev_state = rebuild_reverse_children(rev_state, node_count, goal_idx);
end

function [rev_state, pops] = repair_reverse_after_invalid_edges(rev_state, graph_state, invalid_adj, bad_u, bad_v, start_idx, goal_idx, pop_limit)
affected = unique([bad_u(:); bad_v(:)]).';
if isempty(affected)
    pops = 0;
    return;
end

for i = 1:numel(bad_u)
    u = bad_u(i);
    v = bad_v(i);
    if rev_state.parent(u) == v
        [rev_state, ~] = invalidate_rev_branch(rev_state, u, goal_idx);
    elseif rev_state.parent(v) == u
        [rev_state, ~] = invalidate_rev_branch(rev_state, v, goal_idx);
    end
end

expanded = affected;
for i = 1:numel(affected)
    expanded = [expanded, graph_state.reverse_neighbors{affected(i)}]; %#ok<AGROW>
end
expanded = unique(expanded);

queue = init_vertex_heap(max(4096, graph_state.node_count * 8));
for i = 1:numel(expanded)
    [rev_state, queue] = update_reverse_state(rev_state, graph_state, invalid_adj, expanded(i), start_idx, goal_idx, queue);
end

[rev_state, ~, pops] = compute_reverse_shortest_path(rev_state, graph_state, invalid_adj, queue, start_idx, goal_idx, pop_limit);
rev_state = rebuild_reverse_children(rev_state, graph_state.node_count, goal_idx);
end

function [rev_state, pops] = repair_reverse_after_new_nodes(rev_state, graph_state, invalid_adj, new_idx, start_idx, goal_idx, pop_limit)
new_idx = reshape(new_idx, 1, []);
new_idx = new_idx(new_idx >= 1 & new_idx <= graph_state.node_count);
if isempty(new_idx)
    pops = 0;
    return;
end

affected = new_idx;
for i = 1:numel(new_idx)
    idx = new_idx(i);
    affected = [affected, graph_state.reverse_neighbors{idx}, graph_state.forward_neighbors{idx}]; %#ok<AGROW>
end
affected = unique(affected);

queue = init_vertex_heap(max(4096, graph_state.node_count * 8));
rev_state.rhs(goal_idx) = 0;
for i = 1:numel(affected)
    [rev_state, queue] = update_reverse_state(rev_state, graph_state, invalid_adj, affected(i), start_idx, goal_idx, queue);
end
[queue, rev_state] = queue_reverse_vertex(queue, rev_state, reverse_key(rev_state, graph_state, goal_idx), goal_idx);
[rev_state, ~, pops] = compute_reverse_shortest_path(rev_state, graph_state, invalid_adj, queue, start_idx, goal_idx, pop_limit);
rev_state = rebuild_reverse_children(rev_state, graph_state.node_count, goal_idx);
end

function [rev_state, branch_nodes] = invalidate_rev_branch(rev_state, root_idx, goal_idx)
branch_nodes = zeros(1, 0);
if root_idx == goal_idx
    return;
end

stack = root_idx;
while ~isempty(stack)
    x = stack(end);
    stack(end) = [];
    branch_nodes(end + 1) = x; %#ok<AGROW>
    kids = rev_state.children{x};
    if ~isempty(kids)
        stack = [stack, kids(:).']; %#ok<AGROW>
    end
    rev_state.g(x) = inf;
    rev_state.rhs(x) = inf;
    rev_state.parent(x) = uint32(0);
    rev_state.children{x} = [];
end
end

function rev_state = rebuild_reverse_children(rev_state, node_count, goal_idx)
rev_state.children(1:node_count) = {[]};
for idx = 1:node_count
    p = double(rev_state.parent(idx));
    if idx ~= goal_idx && p >= 1 && p <= node_count
        rev_state.children{p} = [rev_state.children{p}, idx];
    end
end
end

function [rev_state, queue, pops] = compute_reverse_shortest_path(rev_state, graph_state, invalid_adj, queue, start_idx, goal_idx, pop_limit) %#ok<INUSD>
pops = 0;
while queue.size > 0 && pops < pop_limit
    [entry, queue] = pop_vertex(queue);
    if ~rev_state.open_pending(entry.idx)
        continue;
    end
    if ~lex_eq([rev_state.open_key1(entry.idx), rev_state.open_key2(entry.idx)], [entry.key1, entry.key2], 1e-9)
        continue;
    end
    rev_state.open_pending(entry.idx) = false;
    rev_state.open_key1(entry.idx) = inf;
    rev_state.open_key2(entry.idx) = inf;

    curr_key = reverse_key(rev_state, graph_state, entry.idx);
    if ~lex_eq(curr_key, [entry.key1, entry.key2], 1e-9)
        if abs(rev_state.g(entry.idx) - rev_state.rhs(entry.idx)) > 1e-9 || ...
                xor(isinf(rev_state.g(entry.idx)), isinf(rev_state.rhs(entry.idx)))
            [queue, rev_state] = queue_reverse_vertex(queue, rev_state, curr_key, entry.idx);
        end
        continue;
    end

    pops = pops + 1;
    v = entry.idx;
    neigh = graph_state.reverse_neighbors{v};
    if rev_state.g(v) > rev_state.rhs(v)
        rev_state.g(v) = rev_state.rhs(v);
        for k = 1:numel(neigh)
            [rev_state, queue] = update_reverse_state(rev_state, graph_state, invalid_adj, neigh(k), start_idx, goal_idx, queue);
        end
    else
        rev_state.g(v) = inf;
        [rev_state, queue] = update_reverse_state(rev_state, graph_state, invalid_adj, v, start_idx, goal_idx, queue);
        for k = 1:numel(neigh)
            [rev_state, queue] = update_reverse_state(rev_state, graph_state, invalid_adj, neigh(k), start_idx, goal_idx, queue);
        end
    end
end
end

function [rev_state, queue] = update_reverse_state(rev_state, graph_state, invalid_adj, idx, start_idx, goal_idx, queue) %#ok<INUSD>
if idx == goal_idx
    rev_state.rhs(idx) = 0;
    rev_state.parent(idx) = uint32(0);
else
    neigh = graph_state.reverse_neighbors{idx};
    lens = graph_state.reverse_lengths{idx};
    if isempty(neigh)
        best_rhs = inf;
        best_parent = uint32(0);
    else
        bad = logical(full(invalid_adj(idx, neigh)));
        cand = rev_state.g(neigh) + lens;
        cand(bad) = inf;
        [best_rhs, loc] = min(cand);
        if isfinite(best_rhs)
            best_parent = uint32(neigh(loc));
        else
            best_parent = uint32(0);
        end
    end
    rev_state.rhs(idx) = best_rhs;
    rev_state.parent(idx) = best_parent;
end

if abs(rev_state.g(idx) - rev_state.rhs(idx)) > 1e-9 || xor(isinf(rev_state.g(idx)), isinf(rev_state.rhs(idx)))
    [queue, rev_state] = queue_reverse_vertex(queue, rev_state, reverse_key(rev_state, graph_state, idx), idx);
else
    rev_state.open_pending(idx) = false;
    rev_state.open_key1(idx) = inf;
    rev_state.open_key2(idx) = inf;
end
end

function key = reverse_key(rev_state, graph_state, idx)
min_val = min(rev_state.g(idx), rev_state.rhs(idx));
key = [min_val + graph_state.node_to_start(idx), min_val];
end

function [edge_heap, fwd_state] = rebuild_forward_queue(edge_heap, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost)
tree_idx = find(in_tree_fwd(1:graph_state.node_count));
for i = 1:numel(tree_idx)
    [edge_heap, fwd_state] = enqueue_forward_neighbors( ...
        edge_heap, tree_idx(i), in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost);
end
end

function [edge_heap, fwd_state] = enqueue_forward_neighbors(edge_heap, src_idx, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost)
if ~in_tree_fwd(src_idx)
    return;
end

neigh = graph_state.forward_neighbors{src_idx};
lens = graph_state.forward_lengths{src_idx};
for k = 1:numel(neigh)
    dst_idx = neigh(k);
    if src_idx == dst_idx || invalid_adj(src_idx, dst_idx)
        continue;
    end
    edge_len = lens(k);
    key2 = g_fwd(src_idx) + edge_len;
    key1 = key2 + adaptive_goal_heuristic(rev_state, graph_state, dst_idx);
    if isfinite(key1) && key1 < best_goal_cost - 1e-9
        [edge_heap, fwd_state] = queue_forward_edge(edge_heap, fwd_state, key1, key2, src_idx, dst_idx, edge_len);
    end
end
end

function h_val = adaptive_goal_heuristic(rev_state, graph_state, idx)
eucl_goal = graph_state.node_to_goal(idx);
if isfinite(rev_state.g(idx))
    h_val = min(rev_state.g(idx), eucl_goal);
else
    h_val = eucl_goal;
end
end

function h_val = adaptive_goal_heuristic_batch(rev_state, graph_state, idx)
eucl_goal = graph_state.node_to_goal(idx);
h_val = eucl_goal;
finite_mask = isfinite(rev_state.g(idx));
h_val(finite_mask) = min(rev_state.g(idx(finite_mask)), eucl_goal(finite_mask));
end

function [cand, edge_heap, fwd_state, stale_pops, stop_by_bound] = collect_forward_chunk(edge_heap, fwd_state, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, best_goal_cost, chunk_size)
cand = struct();
cand.count = 0;
cand.u = zeros(1, chunk_size);
cand.v = zeros(1, chunk_size);
cand.len = zeros(1, chunk_size);
cand.key1 = inf(1, chunk_size);
cand.key2 = inf(1, chunk_size);
stale_pops = 0;
stop_by_bound = false;

while edge_heap.size > 0 && cand.count < chunk_size
    [entry, edge_heap] = pop_edge(edge_heap);
    if ~full(fwd_state.open_pending(entry.u, entry.v))
        stale_pops = stale_pops + 1;
        continue;
    end
    queued_key = [full(fwd_state.open_key1(entry.u, entry.v)), full(fwd_state.open_key2(entry.u, entry.v))];
    if ~lex_eq(queued_key, [entry.key1, entry.key2], 1e-9)
        stale_pops = stale_pops + 1;
        continue;
    end
    fwd_state = clear_forward_open_edge(fwd_state, entry.u, entry.v);

    if entry.key1 >= best_goal_cost - 1e-9
        stop_by_bound = true;
        break;
    end

    u = entry.u;
    v = entry.v;
    if u < 1 || u > graph_state.node_count || v < 1 || v > graph_state.node_count || u == v
        continue;
    end
    if ~in_tree_fwd(u) || invalid_adj(u, v)
        continue;
    end

    curr_key2 = g_fwd(u) + entry.len;
    curr_key1 = curr_key2 + adaptive_goal_heuristic(rev_state, graph_state, v);
    if ~isfinite(curr_key1) || lex_gt([entry.key1, entry.key2], [curr_key1, curr_key2], 1e-9)
        continue;
    end
    if curr_key1 >= best_goal_cost - 1e-9
        continue;
    end

    cand.count = cand.count + 1;
    cand.u(cand.count) = u;
    cand.v(cand.count) = v;
    cand.len(cand.count) = entry.len;
    cand.key1(cand.count) = entry.key1;
    cand.key2(cand.count) = entry.key2;
end

if cand.count > 0
    cand.u = cand.u(1:cand.count);
    cand.v = cand.v(1:cand.count);
    cand.len = cand.len(1:cand.count);
    cand.key1 = cand.key1(1:cand.count);
    cand.key2 = cand.key2(1:cand.count);
else
    cand.u = zeros(1, 0);
    cand.v = zeros(1, 0);
    cand.len = zeros(1, 0);
    cand.key1 = zeros(1, 0);
    cand.key2 = zeros(1, 0);
end
end

function edge_cost = compute_true_edge_cost_batch(p0, p1, risk_query, risk_limit, stealth_weight, terrain_query, altitude_cost_weight, preferred_agl)
seg_mid = 0.5 * (p0 + p1);
seg_len = vecnorm(p1 - p0, 2, 1);
risk_mid = risk_query(seg_mid);
ok = isfinite(seg_len) & isfinite(risk_mid) & (risk_mid <= risk_limit);
edge_cost = inf(1, size(p0, 2));
alt_cost = compute_altitude_edge_cost_batch(seg_mid, seg_len, terrain_query, altitude_cost_weight, preferred_agl);
if isinf(stealth_weight)
    edge_cost(ok) = seg_len(ok) + alt_cost(ok);
else
    edge_cost(ok) = seg_len(ok) + stealth_weight .* seg_len(ok) .* max(risk_mid(ok), 0) + alt_cost(ok);
end
edge_cost(seg_len < 1e-9) = 0;
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
cost_alt = (altitude_cost_weight .* seg_len(:) .* (err_norm .^ 2)).';
end

function [best_goal_cost, in_tree_fwd, parent_fwd, g_fwd, invalid_adj, fwd_state, checks, successes, edge_heap, eval_time, parallel_used, queue_rebuild_time] = ...
        attempt_goal_probes(best_goal_cost, in_tree_fwd, parent_fwd, g_fwd, nodes, node_count, ...
        rev_state, graph_state, invalid_adj, fwd_state, x_goal, edge_check_samples, ...
        bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance, stealth_weight, altitude_cost_weight, preferred_agl, ...
        goal_probe_count, edge_heap, use_parallel_edge_checks, parallel_threshold_edges)
checks = 0;
successes = 0;
eval_time = 0;
parallel_used = false;
queue_rebuild_time = 0;
if goal_probe_count <= 0
    return;
end

tree_idx = find(in_tree_fwd(1:node_count));
tree_idx = tree_idx(tree_idx ~= 2);
if isempty(tree_idx)
    return;
end

goal_d = graph_state.node_to_goal(tree_idx);
goal_h = adaptive_goal_heuristic_batch(rev_state, graph_state, tree_idx);
rank_key = g_fwd(tree_idx) + goal_h + 0.25 * goal_d;
[~, order] = sort(rank_key, 'ascend');
probe_idx = tree_idx(order(1:min(goal_probe_count, numel(order))));

invalid_mask = get_edge_flags(invalid_adj, probe_idx, repmat(2, 1, numel(probe_idx)));
probe_idx = probe_idx(~invalid_mask);
if isempty(probe_idx)
    return;
end

goal_idx = repmat(2, 1, numel(probe_idx));
cached_mask = get_edge_flags(fwd_state.valid_adj, probe_idx, goal_idx);
ok_mask = cached_mask;
unknown_mask = ~cached_mask;

if any(unknown_mask)
    checks = nnz(unknown_mask);
    t_eval = tic;
    [ok_unknown, parallel_used] = evaluate_edge_safety_chunk( ...
        nodes(:, probe_idx(unknown_mask)), repmat(x_goal(:), 1, checks), edge_check_samples, ...
        bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance, ...
        use_parallel_edge_checks, parallel_threshold_edges);
    eval_time = toc(t_eval);
    ok_mask(unknown_mask) = ok_unknown;
    fwd_state.valid_adj = set_edge_flags(fwd_state.valid_adj, probe_idx(unknown_mask & ok_mask), goal_idx(unknown_mask & ok_mask), true);
end

bad_mask = ~ok_mask;
if any(bad_mask)
    invalid_adj = set_edge_flags(invalid_adj, probe_idx(bad_mask), goal_idx(bad_mask), true);
end

valid_idx = probe_idx(ok_mask);
if isempty(valid_idx)
    return;
end

edge_cost = compute_true_edge_cost_batch(nodes(:, valid_idx), repmat(x_goal(:), 1, numel(valid_idx)), risk_query, risk_limit, stealth_weight, terrain_query, altitude_cost_weight, preferred_agl);
tentative = g_fwd(valid_idx) + edge_cost;
[tentative, order] = sort(tentative, 'ascend');
valid_idx = valid_idx(order);
for i = 1:numel(valid_idx)
    if isfinite(tentative(i)) && tentative(i) + 1e-9 < g_fwd(2)
        in_tree_fwd(2) = true;
        parent_fwd(2) = uint32(valid_idx(i));
        g_fwd(2) = tentative(i);
        best_goal_cost = min(best_goal_cost, tentative(i));
        successes = successes + 1;
    end
end

if successes > 0
    fwd_state = reset_forward_open_state(fwd_state);
    t_queue = tic;
    edge_heap = init_edge_heap(max(4096, node_count * 8));
    [edge_heap, fwd_state] = rebuild_forward_queue( ...
        edge_heap, in_tree_fwd, g_fwd, rev_state, graph_state, invalid_adj, fwd_state, best_goal_cost);
    queue_rebuild_time = toc(t_queue);
end
end

function [ok, parallel_used] = evaluate_edge_safety_chunk(P0, P1, n_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance, use_parallel_edge_checks, parallel_threshold_edges)
n_edges = size(P0, 2);
parallel_used = false;
if n_edges == 0
    ok = false(1, 0);
    return;
end
if ~use_parallel_edge_checks || n_edges < parallel_threshold_edges || ~can_use_parfor()
    ok = edge_is_safe_batch(P0, P1, n_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance);
    return;
end

chunk_size = max(8, ceil(n_edges / min(8, n_edges)));
starts = 1:chunk_size:n_edges;
num_chunks = numel(starts);
local_ok = cell(1, num_chunks);
try
    parfor i = 1:num_chunks
        cols = starts(i):min(n_edges, starts(i) + chunk_size - 1);
        local_ok{i} = edge_is_safe_batch(P0(:, cols), P1(:, cols), n_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance);
    end
    ok = false(1, n_edges);
    for i = 1:num_chunks
        cols = starts(i):min(n_edges, starts(i) + chunk_size - 1);
        ok(cols) = local_ok{i};
    end
    parallel_used = true;
catch
    ok = edge_is_safe_batch(P0, P1, n_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance);
    parallel_used = false;
end
end

function ok = edge_is_safe_batch(P0, P1, n_samples, bounds, risk_query, risk_limit, terrain_query, point_clearance, max_clearance)
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
agl = Pts_flat(3, :) - h;
terrain_ok = isfinite(h) & (agl >= point_clearance) & (agl <= max_clearance);
safe = inb & isfinite(risk) & (risk <= risk_limit) & terrain_ok;
safe = reshape(safe, n_t, n_edges);
ok = all(safe, 1);
end

function flags = get_edge_flags(mat, u, v)
u = double(reshape(u, 1, []));
v = double(reshape(v, 1, []));
if isempty(u)
    flags = false(1, 0);
    return;
end
lin = sub2ind(size(mat), u, v);
flags = full(mat(lin)) ~= 0;
flags = reshape(flags, 1, []);
end

function mat = set_edge_flags(mat, u, v, val)
if isempty(u)
    return;
end
u = double(reshape(u, 1, []));
v = double(reshape(v, 1, []));
lin = sub2ind(size(mat), [u, v], [v, u]);
mat(lin) = val;
end

function path = extract_tree_path(nodes, parents, goal_idx, max_nodes)
chain = zeros(1, max_nodes, 'uint32');
count = 0;
idx = goal_idx;
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
path = nodes(:, chain);
d = vecnorm(diff(path, 1, 2), 2, 1);
keep = [true, d > 1e-9];
path = path(:, keep);
end

function S = sample_mixed_bounds(bounds, x_start, x_goal, n, goal_frac, line_frac)
if n <= 0
    S = zeros(3, 0);
    return;
end

n_goal = min(n, round(goal_frac * n));
n_line = min(n - n_goal, round(line_frac * n));
n_uniform = n - n_goal - n_line;

S = zeros(3, n);
cursor = 0;
if n_uniform > 0
    S(:, cursor + (1:n_uniform)) = sample_uniform_bounds(bounds, n_uniform);
    cursor = cursor + n_uniform;
end
if n_line > 0
    S(:, cursor + (1:n_line)) = sample_line_biased(bounds, x_start, x_goal, n_line);
    cursor = cursor + n_line;
end
if n_goal > 0
    S(:, cursor + (1:n_goal)) = sample_goal_biased(bounds, x_goal, n_goal);
end
end

function S = sample_goal_biased(bounds, x_goal, n)
if n <= 0
    S = zeros(3, 0);
    return;
end

span = [bounds(2) - bounds(1); bounds(4) - bounds(3); bounds(6) - bounds(5)];
sigma = max([0.12 * span(1); 0.12 * span(2); 0.08 * span(3)], [40; 40; 20]);
S = x_goal(:) + sigma .* randn(3, n);
S = clamp_points_to_bounds(S, bounds);
end

function S = sample_line_biased(bounds, x_start, x_goal, n)
if n <= 0
    S = zeros(3, 0);
    return;
end

t = rand(1, n);
base = x_start(:) .* (1 - t) + x_goal(:) .* t;
span = [bounds(2) - bounds(1); bounds(4) - bounds(3); bounds(6) - bounds(5)];
noise = [max(0.08 * span(1), 35); max(0.08 * span(2), 35); max(0.05 * span(3), 15)] .* randn(3, n);
S = base + noise;
S = clamp_points_to_bounds(S, bounds);
end

function P = clamp_points_to_bounds(P, bounds)
P(1, :) = min(max(P(1, :), bounds(1)), bounds(2));
P(2, :) = min(max(P(2, :), bounds(3)), bounds(4));
P(3, :) = min(max(P(3, :), bounds(5)), bounds(6));
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

function r = dynamic_connection_radius(bounds, n, r_min, r_max, eta)
vol = max((bounds(2) - bounds(1)) * (bounds(4) - bounds(3)) * (bounds(6) - bounds(5)), 1e-9);
zeta3 = 4.0 * pi / 3.0;
gamma = 2.0 * (1.0 + 1.0 / 3.0)^(1.0 / 3.0) * (vol / zeta3)^(1.0 / 3.0);
n_eff = max(2, n);
r_asym = eta * gamma * (log(n_eff) / n_eff)^(1.0 / 3.0);
r = min(r_max, max(r_min, r_asym));
end

function tf = can_use_kdtree()
tf = (exist('createns', 'builtin') > 0 || exist('createns', 'file') > 0) && ...
     (exist('rangesearch', 'builtin') > 0 || exist('rangesearch', 'file') > 0);
end

function tf = can_use_parfor()
tf = license('test', 'Distrib_Computing_Toolbox') && ...
     (exist('parfor', 'builtin') > 0 || exist('parfor', 'file') > 0);
end

function tf = lex_eq(a, b, tol)
tf = abs(a(1) - b(1)) <= tol && abs(a(2) - b(2)) <= tol;
end

function tf = lex_gt(a, b, tol)
if a(1) > b(1) + tol
    tf = true;
elseif a(1) < b(1) - tol
    tf = false;
else
    tf = a(2) > b(2) + tol;
end
end

function heap = init_vertex_heap(capacity)
heap = struct();
heap.key1 = inf(capacity, 1);
heap.key2 = inf(capacity, 1);
heap.idx = zeros(capacity, 1, 'uint32');
heap.size = 0;
end

function heap = push_vertex(heap, key, idx)
heap.size = heap.size + 1;
if heap.size > numel(heap.key1)
    new_cap = max(heap.size, 2 * numel(heap.key1));
    heap.key1(new_cap, 1) = inf;
    heap.key2(new_cap, 1) = inf;
    heap.idx(new_cap, 1) = uint32(0);
end

pos = heap.size;
heap.key1(pos) = key(1);
heap.key2(pos) = key(2);
heap.idx(pos) = uint32(idx);
while pos > 1
    parent = floor(pos / 2);
    if ~lex_gt([heap.key1(parent), heap.key2(parent)], [heap.key1(pos), heap.key2(pos)], 0)
        break;
    end
    heap = swap_vertex_heap(heap, parent, pos);
    pos = parent;
end
end

function [heap, rev_state] = queue_reverse_vertex(heap, rev_state, key, idx)
if rev_state.open_pending(idx)
    curr = [rev_state.open_key1(idx), rev_state.open_key2(idx)];
    if ~lex_gt(curr, key, 1e-9)
        return;
    end
end
heap = push_vertex(heap, key, idx);
rev_state.open_pending(idx) = true;
rev_state.open_key1(idx) = key(1);
rev_state.open_key2(idx) = key(2);
end

function [entry, heap] = pop_vertex(heap)
entry = struct('key1', heap.key1(1), 'key2', heap.key2(1), 'idx', double(heap.idx(1)));
heap.key1(1) = heap.key1(heap.size);
heap.key2(1) = heap.key2(heap.size);
heap.idx(1) = heap.idx(heap.size);
heap.size = heap.size - 1;

pos = 1;
while true
    left = 2 * pos;
    right = left + 1;
    if left > heap.size
        break;
    end
    best = left;
    if right <= heap.size && lex_gt([heap.key1(left), heap.key2(left)], [heap.key1(right), heap.key2(right)], 0)
        best = right;
    end
    if ~lex_gt([heap.key1(pos), heap.key2(pos)], [heap.key1(best), heap.key2(best)], 0)
        break;
    end
    heap = swap_vertex_heap(heap, pos, best);
    pos = best;
end
end

function heap = swap_vertex_heap(heap, a, b)
[heap.key1(a), heap.key1(b)] = deal(heap.key1(b), heap.key1(a));
[heap.key2(a), heap.key2(b)] = deal(heap.key2(b), heap.key2(a));
[heap.idx(a), heap.idx(b)] = deal(heap.idx(b), heap.idx(a));
end

function heap = init_edge_heap(capacity)
heap = struct();
heap.key1 = inf(capacity, 1);
heap.key2 = inf(capacity, 1);
heap.len = zeros(capacity, 1);
heap.u = zeros(capacity, 1, 'uint32');
heap.v = zeros(capacity, 1, 'uint32');
heap.size = 0;
end

function heap = push_edge(heap, key1, key2, u, v, edge_len)
heap.size = heap.size + 1;
if heap.size > numel(heap.key1)
    new_cap = max(heap.size, 2 * numel(heap.key1));
    heap.key1(new_cap, 1) = inf;
    heap.key2(new_cap, 1) = inf;
    heap.len(new_cap, 1) = 0;
    heap.u(new_cap, 1) = uint32(0);
    heap.v(new_cap, 1) = uint32(0);
end

pos = heap.size;
heap.key1(pos) = key1;
heap.key2(pos) = key2;
heap.len(pos) = edge_len;
heap.u(pos) = uint32(u);
heap.v(pos) = uint32(v);
while pos > 1
    parent = floor(pos / 2);
    if ~lex_gt([heap.key1(parent), heap.key2(parent)], [heap.key1(pos), heap.key2(pos)], 0)
        break;
    end
    heap = swap_edge_heap(heap, parent, pos);
    pos = parent;
end
end

function [heap, fwd_state] = queue_forward_edge(heap, fwd_state, key1, key2, u, v, edge_len)
if full(fwd_state.open_pending(u, v))
    curr = [full(fwd_state.open_key1(u, v)), full(fwd_state.open_key2(u, v))];
    if ~lex_gt(curr, [key1, key2], 1e-9)
        return;
    end
end
heap = push_edge(heap, key1, key2, u, v, edge_len);
fwd_state.open_pending(u, v) = true;
fwd_state.open_key1(u, v) = key1;
fwd_state.open_key2(u, v) = key2;
end

function [entry, heap] = pop_edge(heap)
entry = struct( ...
    'key1', heap.key1(1), ...
    'key2', heap.key2(1), ...
    'len', heap.len(1), ...
    'u', double(heap.u(1)), ...
    'v', double(heap.v(1)));
heap.key1(1) = heap.key1(heap.size);
heap.key2(1) = heap.key2(heap.size);
heap.len(1) = heap.len(heap.size);
heap.u(1) = heap.u(heap.size);
heap.v(1) = heap.v(heap.size);
heap.size = heap.size - 1;

pos = 1;
while true
    left = 2 * pos;
    right = left + 1;
    if left > heap.size
        break;
    end
    best = left;
    if right <= heap.size && lex_gt([heap.key1(left), heap.key2(left)], [heap.key1(right), heap.key2(right)], 0)
        best = right;
    end
    if ~lex_gt([heap.key1(pos), heap.key2(pos)], [heap.key1(best), heap.key2(best)], 0)
        break;
    end
    heap = swap_edge_heap(heap, pos, best);
    pos = best;
end
end

function heap = swap_edge_heap(heap, a, b)
[heap.key1(a), heap.key1(b)] = deal(heap.key1(b), heap.key1(a));
[heap.key2(a), heap.key2(b)] = deal(heap.key2(b), heap.key2(a));
[heap.len(a), heap.len(b)] = deal(heap.len(b), heap.len(a));
[heap.u(a), heap.u(b)] = deal(heap.u(b), heap.u(a));
[heap.v(a), heap.v(b)] = deal(heap.v(b), heap.v(a));
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
