%% TEST_AIT_STAR_DIRECT_SMOKE
% Smoke test for direct AIT* planner integration.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
setup_project_paths();

rng(11, 'twister');

[X, Y] = meshgrid(linspace(0, 2500, 55), linspace(0, 2500, 55));
Z = 180 + 90 * sin(X / 420) .* cos(Y / 360);
terrain_data = struct();
terrain_data.N_vec = X(1, :);
terrain_data.E_vec = Y(:, 1)';
terrain_data.Z = Z;
terrain_data.bounds = [min(terrain_data.N_vec), max(terrain_data.N_vec), ...
                       min(terrain_data.E_vec), max(terrain_data.E_vec)];
terrain_data.resolution = mean(diff(terrain_data.N_vec));
terrain_data.type = 'synthetic_ait_test';

tm = terrain_map(terrain_data);
los = los_checker(tm);

diag_len = hypot(diff(tm.bounds(1:2)), diff(tm.bounds(3:4)));
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));

threat_params = struct();
threat_params.resolution = [max(35, diag_len / 90), 35];
threat_params.alt_range = [terrain_min + 1, terrain_max + 250];
threat = threat_map(tm, los, 0.1, threat_params);

radar_pos = [tm.bounds(1) + 0.55 * diff(tm.bounds(1:2)); ...
             tm.bounds(3) + 0.48 * diff(tm.bounds(3:4))];
rH = tm.get_height(radar_pos(1), radar_pos(2)) + 15;
radar_obj = radar_site([radar_pos; rH], 'SmokeRadar', 'tracking');
radar_obj.R_max = 1400;
radar_obj.P_t = 120e3;
threat.add_radar(radar_obj);
threat.compute_map('binary', struct('use_parallel', false, 'show_progress', false));

start_NE = [tm.bounds(1) + 0.08 * diff(tm.bounds(1:2)); tm.bounds(3) + 0.10 * diff(tm.bounds(3:4))];
goal_NE = [tm.bounds(2) - 0.10 * diff(tm.bounds(1:2)); tm.bounds(4) - 0.12 * diff(tm.bounds(3:4))];
start_h = tm.get_height(start_NE(1), start_NE(2));
goal_h = tm.get_height(goal_NE(1), goal_NE(2));
start_world = [start_NE; start_h + 70];
goal_world = [goal_NE; goal_h + 70];

[Nq, Eq] = meshgrid(threat.N_vec, threat.E_vec);
Zq = tm.get_height(Nq, Eq);
env = struct();
env.terrain = struct( ...
    'Z', Zq, ...
    'N_vec', threat.N_vec, ...
    'E_vec', threat.E_vec, ...
    'x_vec', threat.N_vec, ...
    'y_vec', threat.E_vec, ...
    'dx', mean(diff(threat.N_vec)), ...
    'dy', mean(diff(threat.E_vec)), ...
    'alt_vec', threat.alt_vec);
env.risk_grid = threat.risk_grid;
env.radar_los_map = threat.risk_grid >= 0.5;
env.N_vec = threat.N_vec;
env.E_vec = threat.E_vec;
env.alt_vec = threat.alt_vec;
env.visibility_threshold = 0.5;
env.bounds_world = [threat.N_vec(1), threat.N_vec(end), threat.E_vec(1), threat.E_vec(end), ...
                    threat.alt_vec(1), threat.alt_vec(end)];

params = struct();
params.max_batches = 18;
params.batch_size = 90;
params.max_nodes = 2200;
params.eps_global = 0.18;
params.stealth_weight = inf;
params.strict_zero_risk_tol = 1e-6;
params.edge_check_samples = 7;
params.post_solution_batches = 4;
params.goal_radius = 90;
params.min_connection_radius = 70;
params.max_connection_radius = 220;
params.eta_radius = 1.8;
params.max_time = 12;
params.min_clearance = 10;
params.reverse_k_neighbors = 20;
params.forward_k_neighbors = 14;
params.goal_sample_fraction = 0.22;
params.line_sample_fraction = 0.28;
params.goal_probe_count = 8;
params.profile_ait = false;
params.neighbor_precompute_mode = 'batch';
params.use_parallel_reverse = false;
params.use_parallel_edge_checks = false;
params.parallel_threshold_neighbors = 64;
params.parallel_threshold_edges = 32;
params.reverse_batch_repair_limit = 128;
params.forward_expand_chunk = 16;
params.kd_rebuild_interval = 1;
params.reverse_debug = false;

[path_world, cost, info] = ait_star_planner_direct(start_world, goal_world, env, params);

assert(isstruct(info), 'AIT* smoke test expected info struct output.');
assert(islogical(info.success) && isscalar(info.success), 'AIT* smoke test expected logical scalar success flag.');

if info.success
    assert(~isempty(path_world) && size(path_world, 2) >= 2, 'Successful AIT* run must return a path with at least 2 waypoints.');
    n_seg_samples = max(2, 8 * (size(path_world, 2) - 1));
    path_samples = interpolate_path(path_world, n_seg_samples);
    terrain_h = tm.get_height(path_samples(1, :), path_samples(2, :));
    assert(all(isfinite(terrain_h)), 'Path samples must have finite terrain heights.');
    assert(all(path_samples(3, :) >= terrain_h + params.min_clearance - 1e-6), ...
        'Successful AIT* path violates terrain clearance.');
    risk_vals = arrayfun(@(k) threat.get_risk(path_samples(1, k), path_samples(2, k), path_samples(3, k)), 1:size(path_samples, 2));
    assert(all(isfinite(risk_vals)), 'Path samples must have finite risk values.');
    assert(all(risk_vals <= params.strict_zero_risk_tol + 1e-6), ...
        'Strict zero-risk AIT* path violates radar threshold.');
    assert(isfinite(cost), 'Successful AIT* run must report finite cost.');
else
    assert(isempty(path_world), 'Failed AIT* run should not return a path.');
end

fprintf('AIT* smoke test completed. success=%d, iterations=%d, invalid_edges=%d\n', ...
    info.success, info.iterations, info.invalid_edge_count);

function samples = interpolate_path(path_world, min_samples)
if size(path_world, 2) < 2
    samples = path_world;
    return;
end

segments = size(path_world, 2) - 1;
samples = path_world(:, 1);
for i = 1:segments
    p0 = path_world(:, i);
    p1 = path_world(:, i + 1);
    n_local = max(2, ceil(min_samples / segments));
    lam = linspace(0, 1, n_local);
    seg = p0 .* (1 - lam) + p1 .* lam;
    if i > 1
        seg = seg(:, 2:end);
    end
    samples = [samples, seg]; %#ok<AGROW>
end
end
