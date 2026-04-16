%% DEMO_RADAR_PLANNER_BENCHMARK
% Unified benchmark for FSM/BIT*, FMM/BIT*, and radar-aware RRT*.

clear; clc; close all;
clear fmm_bit_star_planner fsm_bit_star_planner_direct rrt_star_radar compute_stealth_corridor;
rehash;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
setup_project_paths();

fprintf('=== Unified Radar Planner Benchmark ===\n');

cfg = struct();
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 30;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 2500;
cfg.visibility_threshold = 0.5;
cfg.start_goal_agl = 80;
cfg.threat_vert_res = 40;
cfg.use_parallel_threat = false;
cfg.radar_alt_offset = 20;
cfg.radar_range = 3000;
cfg.radar_power = 150e3;
cfg.radar_fraction_ne = [0.22, 0.78; 0.25, 0.72];
cfg.corridor_alpha = 2.3;
cfg.corridor_lambda = 1.8;
cfg.corridor_r_min = 120;
cfg.corridor_r_max = 450;
cfg.corridor_max_sweeps = 60;
cfg.common_stealth_alpha = 2.0;
cfg.strict_zero_risk_tol = 1e-6;
cfg.bit_strict_zero_risk = false;
cfg.bit_stealth_weight = 3.0;
cfg.fsm_max_sweeps = 96;
cfg.fsm_tol = 1e-4;
cfg.rng_seed = 7;
cfg.colors = struct( ...
    'fsm_bit', [0.90, 0.20, 0.18], ...
    'fmm_bit', [0.12, 0.55, 0.86], ...
    'rrt_star', [0.15, 0.68, 0.36], ...
    'radar', [1.00, 0.52, 0.10], ...
    'terrain', [0.67, 0.61, 0.54]);

rng(cfg.rng_seed, 'twister');

setup_tic = tic;

%% 1) Terrain
t_dem = tic;
[terrain_data, dem_source] = load_benchmark_terrain(cfg);
t_dem_s = toc(t_dem);

tm = terrain_map(terrain_data);
los = los_checker(tm);

terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);

%% 2) Radar / threat volume
threat_params = struct();
threat_params.resolution = [max(40, diag_len / 120), cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];

threat = threat_map(tm, los, 0.1, threat_params);
radar_ne = [tm.bounds(1) + cfg.radar_fraction_ne(1, :) * N_span; ...
            tm.bounds(3) + cfg.radar_fraction_ne(2, :) * E_span];

for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + cfg.radar_alt_offset;
    radar_obj = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar_obj.R_max = cfg.radar_range;
    radar_obj.P_t = cfg.radar_power;
    threat.add_radar(radar_obj);
end

fprintf('Computing binary threat map...\n');
t_threat = tic;
threat.compute_map('binary', struct('use_parallel', cfg.use_parallel_threat, 'show_progress', false));
t_threat_s = toc(t_threat);
fprintf('Threat map compute time: %.2f s\n', t_threat_s);

radar_los_map = threat.risk_grid >= cfg.visibility_threshold;

%% 3) Shared start / goal
[start_NE, start_h, start_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(1) + 0.08 * N_span; tm.bounds(3) + 0.10 * E_span], ...
    cfg.start_goal_agl, cfg.visibility_threshold);
[goal_NE, goal_h, goal_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(2) - 0.08 * N_span; tm.bounds(4) - 0.10 * E_span], ...
    cfg.start_goal_agl, cfg.visibility_threshold);

start_world = [start_NE; start_h + cfg.start_goal_agl];
goal_world = [goal_NE; goal_h + cfg.start_goal_agl];
start_ned = world_to_ned(start_world);
goal_ned = world_to_ned(goal_world);

fprintf('Start risk: %.2f | Goal risk: %.2f\n', start_risk, goal_risk);

shared_times = struct();
shared_times.dem_load_s = t_dem_s;
shared_times.threat_compute_s = t_threat_s;
shared_times.total_setup_s = toc(setup_tic);

%% 4) Shared benchmark metadata
planner_benchmark = struct();
planner_benchmark.cfg = cfg;
planner_benchmark.environment = build_environment_summary( ...
    cfg, tm, threat, dem_source, radar_los_map, radar_ne, ...
    start_world, goal_world, start_ned, goal_ned, start_risk, goal_risk);
planner_benchmark.shared_times = shared_times;
planner_benchmark.planners = struct();

%% 5) FSM/BIT*
fprintf('\n[1/3] Running FSM/BIT*...\n');
planner_benchmark.planners.fsm_bit = run_fsm_bit_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len);

%% 6) FMM/BIT*
fprintf('\n[2/3] Running FMM/BIT*...\n');
planner_benchmark.planners.fmm_bit = run_fmm_bit_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len);

%% 7) RRT*
fprintf('\n[3/3] Running RRT*...\n');
planner_benchmark.planners.rrt_star = run_rrt_star_benchmark(cfg, tm, threat, start_ned, goal_ned, terrain_min, terrain_max, diag_len);

%% 8) Comparison arrays
planner_benchmark.comparison = build_comparison_struct(planner_benchmark);

%% 9) Figures
plot_benchmark_metrics(planner_benchmark);
plot_overlay_figure(planner_benchmark, tm, threat);
plot_planner_panels(planner_benchmark, tm, threat);

%% 10) Console summary and workspace export
print_benchmark_summary(planner_benchmark);
assignin('base', 'planner_benchmark', planner_benchmark);

%% ---------------- Local helpers ----------------
function [terrain_data, dem_source] = load_benchmark_terrain(cfg)
terrain_data = [];
dem_source = 'synthetic_fallback';

try
    if exist(cfg.dem_file, 'file')
        dem_params = struct();
        dem_params.target_resolution = cfg.dem_target_resolution;
        dem_params.fill_nodata = cfg.dem_fill_nodata;
        if ~isempty(cfg.dem_crop_half_size) && cfg.dem_crop_half_size > 0
            h = cfg.dem_crop_half_size;
            dem_params.crop_bounds = [-h, h, -h, h];
        end
        terrain_data = dem_loader(cfg.dem_file, dem_params);
        dem_source = cfg.dem_file;
    end
catch ME
    warning('DEM load failed (%s). Using synthetic terrain fallback.', ME.message);
    terrain_data = [];
end

if isempty(terrain_data)
    [X, Y] = meshgrid(linspace(0, 5000, 100), linspace(0, 5000, 100));
    Z = 200 + 150 * sin(X / 800) .* cos(Y / 600) + 80 * randn(size(X)) * 0.05;
    terrain_data = struct();
    terrain_data.N_vec = X(1, :);
    terrain_data.E_vec = Y(:, 1)';
    terrain_data.Z = Z;
    terrain_data.bounds = [min(terrain_data.N_vec), max(terrain_data.N_vec), ...
                           min(terrain_data.E_vec), max(terrain_data.E_vec)];
    terrain_data.resolution = mean(diff(terrain_data.N_vec));
    terrain_data.type = 'synthetic_benchmark';
end
end

function env_summary = build_environment_summary(cfg, tm, threat, dem_source, radar_los_map, radar_ne, start_world, goal_world, start_ned, goal_ned, start_risk, goal_risk)
radar_sites = repmat(struct('name', '', 'position', [], 'range_m', NaN), 1, numel(threat.radars));
for i = 1:numel(threat.radars)
    rr = threat.radars{i};
    radar_sites(i).name = rr.name;
    radar_sites(i).position = rr.position;
    radar_sites(i).range_m = rr.R_max;
end

env_summary = struct();
env_summary.dem_source = dem_source;
env_summary.terrain_type = getfield_with_default(tm, 'type', '');
env_summary.terrain_bounds = tm.bounds;
env_summary.terrain_resolution = getfield_with_default(tm, 'resolution', mean(diff(tm.N_vec)));
env_summary.threat_resolution = threat.resolution;
env_summary.threat_alt_range = [threat.alt_vec(1), threat.alt_vec(end)];
env_summary.threat_grid_size = size(threat.risk_grid);
env_summary.visibility_threshold = cfg.visibility_threshold;
env_summary.radar_fraction_ne = cfg.radar_fraction_ne;
env_summary.radar_positions_ne = radar_ne;
env_summary.radar_sites = radar_sites;
env_summary.radar_los_map_size = size(radar_los_map);
env_summary.start_world = start_world;
env_summary.goal_world = goal_world;
env_summary.start_ned = start_ned;
env_summary.goal_ned = goal_ned;
env_summary.start_risk = start_risk;
env_summary.goal_risk = goal_risk;
end

function result = run_fsm_bit_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len)
[Nq_fsm, Eq_fsm] = meshgrid(threat.N_vec, threat.E_vec);
Z_fsm = tm.get_height(Nq_fsm, Eq_fsm);

env = struct();
env.terrain = struct( ...
    'Z', Z_fsm, ...
    'N_vec', threat.N_vec, ...
    'E_vec', threat.E_vec, ...
    'x_vec', threat.N_vec, ...
    'y_vec', threat.E_vec, ...
    'dx', mean(diff(threat.N_vec)), ...
    'dy', mean(diff(threat.E_vec)), ...
    'alt_vec', threat.alt_vec);
env.risk_grid = threat.risk_grid;
env.radar_los_map = radar_los_map;
env.N_vec = threat.N_vec;
env.E_vec = threat.E_vec;
env.alt_vec = threat.alt_vec;
env.visibility_threshold = cfg.visibility_threshold;
env.bounds_world = [threat.N_vec(1), threat.N_vec(end), threat.E_vec(1), threat.E_vec(end), ...
                    threat.alt_vec(1), threat.alt_vec(end)];

planner_params = struct();
planner_params.max_batches = 120;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.12;
if cfg.bit_strict_zero_risk
    planner_params.stealth_weight = inf;
else
    planner_params.stealth_weight = cfg.bit_stealth_weight;
end
planner_params.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
planner_params.edge_check_samples = 9;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.eta_radius = 1.8;
planner_params.max_time = 45;
planner_params.min_clearance = 10;
planner_params.w_mask = 1.0;
planner_params.w_alt = 0.08;
planner_params.h_opt = cfg.start_goal_agl;
planner_params.h_scale = 200;
planner_params.dz = mean(diff(threat.alt_vec));
planner_params.mem_limit_gb = 1.5;
planner_params.fsm_point_clearance = 10;
planner_params.max_sweeps = cfg.fsm_max_sweeps;
planner_params.tol = cfg.fsm_tol;
planner_params.fsm_bounds_pad_xy = 180;
planner_params.fsm_bounds_pad_z = 80;

result = init_planner_result('FSM/BIT*');
result.strict_zero_risk = cfg.bit_strict_zero_risk;
t_total = tic;
try
    [path_world, native_cost, info] = fsm_bit_star_planner_direct(start_world, goal_world, env, planner_params);
    result = finalize_world_path_result(result, path_world, threat, cfg.common_stealth_alpha);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = finite_or_nan(info.fsm_time_s);
    result.heuristic_time_s = finite_or_nan(info.fsm_time_s);
    result.search_time_s = finite_or_nan(info.planning_time);
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = read_info_text(info, 'stop_reason');
    result.native_path_cost = finite_or_nan(native_cost);
    result.warmstart_cost = read_info_field(info, 'warmstart_cost');
    result.volume_ratio = read_info_field(info, 'volume_ratio');
    result.fsm_used = read_info_field(info, 'fsm_used');
    result.samples_generated = read_info_field(info, 'samples_generated');
    result.samples_accepted = read_info_field(info, 'samples_accepted');
    result.connection_radius_last = read_info_field(info, 'connection_radius_last');
    result.first_solution_batch = read_info_field(info, 'first_solution_batch');
    result.native_summary = extract_native_summary(info, {'risk_limit', 'max_batches', 'post_solution_batches'});
    if isfield(info, 'fsm_meta') && ~isempty(info.fsm_meta)
        result.native_summary.fsm_sweeps_fwd = read_info_field(info.fsm_meta, 'sweeps_fwd');
        result.native_summary.fsm_sweeps_bwd = read_info_field(info.fsm_meta, 'sweeps_bwd');
        result.native_summary.fsm_delta_fwd = read_info_field(info.fsm_meta, 'final_delta_fwd');
        result.native_summary.fsm_delta_bwd = read_info_field(info.fsm_meta, 'final_delta_bwd');
    end
    if ~result.success
        result.error_message = derive_failure_message(info, result);
    end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = run_fmm_bit_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len)
result = init_planner_result('FMM/BIT*');
result.strict_zero_risk = cfg.bit_strict_zero_risk;

corr_params = struct();
corr_params.N_vec = threat.N_vec;
corr_params.E_vec = threat.E_vec;
corr_params.alt_vec = threat.alt_vec;
corr_params.max_sweeps = cfg.corridor_max_sweeps;
corr_params.lambda = cfg.corridor_lambda;
corr_params.r_min = cfg.corridor_r_min;
corr_params.r_max = cfg.corridor_r_max;
corr_params.max_backtrack_points = 5000;

[Nq_fmm, Eq_fmm] = meshgrid(threat.N_vec, threat.E_vec);
Z_fmm = tm.get_height(Nq_fmm, Eq_fmm);

planner_params = struct();
planner_params.max_batches = 120;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.08;
if cfg.bit_strict_zero_risk
    planner_params.stealth_weight = inf;
else
    planner_params.stealth_weight = cfg.bit_stealth_weight;
end
planner_params.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
planner_params.edge_check_samples = 9;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.eta_radius = 1.8;
planner_params.max_time = 45;
planner_params.min_clearance = 10;
planner_params.w_mask = 1.0;
planner_params.w_alt = 0.08;
planner_params.h_opt = cfg.start_goal_agl;
planner_params.h_scale = 200;
planner_params.dz = mean(diff(threat.alt_vec));
planner_params.mem_limit_gb = 1.5;
planner_params.fmm_point_clearance = 10;
planner_params.fmm_bounds_pad_xy = 180;
planner_params.fmm_bounds_pad_z = 80;
planner_params.fmm_terrain_map = struct( ...
    'Z', Z_fmm, ...
    'N_vec', threat.N_vec, ...
    'E_vec', threat.E_vec, ...
    'x_vec', threat.N_vec, ...
    'y_vec', threat.E_vec, ...
    'dx', mean(diff(threat.N_vec)), ...
    'dy', mean(diff(threat.E_vec)), ...
    'alt_vec', threat.alt_vec);
planner_params.fmm_radar_los_map = radar_los_map;

t_total = tic;
try
    t_corr = tic;
    corridor = compute_stealth_corridor( ...
        threat.risk_grid, start_world, goal_world, cfg.visibility_threshold, cfg.corridor_alpha, corr_params);
    result.corridor_time_s = toc(t_corr);
    goal_idx = round(corridor.grid.goal_idx);
    result.corridor_goal_connected = isfinite(corridor.T(goal_idx(1), goal_idx(2), goal_idx(3)));
    result.corridor_skeleton_points = size(corridor.skeleton, 2);

    [path_world, native_cost, info] = fmm_bit_star_planner(start_world, goal_world, corridor, planner_params);
    result = finalize_world_path_result(result, path_world, threat, cfg.common_stealth_alpha);
    result.success = logical(info.success);
    result.error_message = '';
    result.heuristic_time_s = finite_or_nan(info.fmm_time_s);
    result.preprocess_time_s = finite_or_nan(result.corridor_time_s + finite_or_zero(info.fmm_time_s));
    result.search_time_s = finite_or_nan(info.planning_time);
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = read_info_text(info, 'stop_reason');
    result.native_path_cost = finite_or_nan(native_cost);
    result.warmstart_cost = read_info_field(info, 'warmstart_cost');
    result.volume_ratio = read_info_field(info, 'volume_ratio');
    result.fmm_used = read_info_field(info, 'fmm_used');
    result.samples_generated = read_info_field(info, 'samples_generated');
    result.samples_accepted = read_info_field(info, 'samples_accepted');
    result.connection_radius_last = read_info_field(info, 'connection_radius_last');
    result.first_solution_batch = read_info_field(info, 'first_solution_batch');
    result.native_summary = extract_native_summary(info, {'risk_limit', 'max_batches', 'post_solution_batches'});
    if ~result.success
        result.error_message = derive_failure_message(info, result);
    end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = run_rrt_star_benchmark(cfg, tm, threat, start_ned, goal_ned, terrain_min, terrain_max, diag_len)
planner_params = struct();
planner_params.max_iter = cfg.dem_crop_half_size * 6;
planner_params.base_step_size = max(40, min(75, diag_len / 85));
planner_params.max_step_size = planner_params.base_step_size * 3.5;
planner_params.rewire_radius = planner_params.base_step_size * 3;
planner_params.rewire_mode = 'fixed';
planner_params.rewire_k_const = 3.2 * exp(1) * (1 + 1 / 3);
planner_params.rewire_k_min = 50;
planner_params.rewire_k_max = 350;
planner_params.goal_bias = 0.15;
planner_params.goal_bias_after_goal = 0.04;
planner_params.min_clearance = 20;
planner_params.max_flight_alt = 500;
planner_params.alpha = 1.0;
planner_params.beta = 0.0;
planner_params.gamma = 0.2;
planner_params.preferred_agl = 40;
planner_params.shadow_bias = 0.85;
planner_params.use_parallel_rewire = false;
planner_params.animate = false;
planner_params.plot_interval = 500;
planner_params.height_query_mode = 'lut';
planner_params.radar_hard_constraint = true;
planner_params.radar_visibility_threshold = cfg.visibility_threshold;

min_alt = terrain_min + planner_params.min_clearance;
max_alt = terrain_max + planner_params.max_flight_alt;
planner_params.bounds = [tm.bounds(1), tm.bounds(2), tm.bounds(3), tm.bounds(4), -max_alt, -min_alt];

result = init_planner_result('RRT*');
t_total = tic;
try
    [path_ned, info] = rrt_star_radar(start_ned, goal_ned, tm, threat, planner_params);
    result = finalize_ned_path_result(result, path_ned, threat, cfg.common_stealth_alpha);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = 0.0;
    result.heuristic_time_s = 0.0;
    result.search_time_s = read_info_field(info, 'planning_time');
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = infer_rrt_stop_reason(info);
    result.native_path_cost = read_info_field(info, 'path_cost');
    result.terminated_by_time_limit = read_info_field(info, 'terminated_by_time_limit');
    result.native_summary = extract_native_summary(info, {'first_solution_iter', 'path_risk'});
    if ~result.success
        result.error_message = derive_failure_message(info, result);
    end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = init_planner_result(label)
result = struct();
result.label = label;
result.success = false;
result.error_message = 'Planner did not run.';
result.path_world = [];
result.path_ned = [];
result.num_waypoints = 0;
result.preprocess_time_s = NaN;
result.heuristic_time_s = NaN;
result.search_time_s = NaN;
result.total_time_s = NaN;
result.first_solution_time_s = NaN;
result.iterations = NaN;
result.tree_size = NaN;
result.stop_reason = '';
result.path_length_m = NaN;
result.risk_max = NaN;
result.risk_avg = NaN;
result.risk_integrated = NaN;
result.common_stealth_cost = NaN;
result.strict_zero_risk = NaN;
result.native_path_cost = NaN;
result.warmstart_cost = NaN;
result.volume_ratio = NaN;
result.fsm_used = NaN;
result.fmm_used = NaN;
result.corridor_time_s = NaN;
result.corridor_goal_connected = NaN;
result.corridor_skeleton_points = NaN;
result.samples_generated = NaN;
result.samples_accepted = NaN;
result.connection_radius_last = NaN;
result.first_solution_batch = NaN;
result.terminated_by_time_limit = NaN;
result.native_summary = struct();
end

function result = finalize_world_path_result(result, path_world, threat, alpha_cost)
result.path_world = path_world;
if isempty(path_world)
    result.path_ned = [];
    result.num_waypoints = 0;
    return;
end

path_ned = [path_world(1, :); path_world(2, :); -path_world(3, :)];
result = finalize_ned_path_result(result, path_ned, threat, alpha_cost);
end

function result = finalize_ned_path_result(result, path_ned, threat, alpha_cost)
result.path_ned = path_ned;
if isempty(path_ned)
    result.path_world = [];
    result.num_waypoints = 0;
    return;
end

result.path_world = [path_ned(1, :); path_ned(2, :); -path_ned(3, :)];
result.num_waypoints = size(path_ned, 2);
result.path_length_m = compute_path_length(path_ned);
[result.risk_max, result.risk_avg, result.risk_integrated] = threat.evaluate_path(path_ned, 1);
result.common_stealth_cost = compute_common_stealth_cost(path_ned, threat, alpha_cost);
end

function comparison = build_comparison_struct(planner_benchmark)
keys = {'fsm_bit', 'fmm_bit', 'rrt_star'};
labels = {'FSM/BIT*', 'FMM/BIT*', 'RRT*'};
colors = [planner_benchmark.cfg.colors.fsm_bit; ...
          planner_benchmark.cfg.colors.fmm_bit; ...
          planner_benchmark.cfg.colors.rrt_star];

n = numel(keys);
success = false(1, n);
total_time_s = nan(1, n);
preprocess_time_s = nan(1, n);
search_time_s = nan(1, n);
path_length_m = nan(1, n);
risk_integrated = nan(1, n);
common_stealth_cost = nan(1, n);

for i = 1:n
    item = planner_benchmark.planners.(keys{i});
    success(i) = logical(item.success);
    total_time_s(i) = item.total_time_s;
    preprocess_time_s(i) = item.preprocess_time_s;
    search_time_s(i) = item.search_time_s;
    path_length_m(i) = item.path_length_m;
    risk_integrated(i) = item.risk_integrated;
    common_stealth_cost(i) = item.common_stealth_cost;
end

comparison = struct();
comparison.keys = keys;
comparison.labels = labels;
comparison.colors = colors;
comparison.success = success;
comparison.total_time_s = total_time_s;
comparison.preprocess_time_s = preprocess_time_s;
comparison.search_time_s = search_time_s;
comparison.path_length_m = path_length_m;
comparison.risk_integrated = risk_integrated;
comparison.common_stealth_cost = common_stealth_cost;
comparison.time_stack_s = [preprocess_time_s(:), search_time_s(:)];
end

function plot_benchmark_metrics(planner_benchmark)
cmp = planner_benchmark.comparison;

fig = figure('Name', 'Planner Benchmark Metrics', ...
    'Position', [90, 70, 1400, 820], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Unified Radar Planner Benchmark');

ax1 = nexttile(tl, 1);
plot_scalar_metric(ax1, cmp.total_time_s, cmp.labels, cmp.colors, 'Total Planner Time [s]', 'Total Time');

ax2 = nexttile(tl, 2);
plot_stacked_time_metric(ax2, cmp.time_stack_s, cmp.labels, cmp.colors, cmp.success);

ax3 = nexttile(tl, 3);
plot_scalar_metric(ax3, cmp.path_length_m, cmp.labels, cmp.colors, 'Path Length [m]', 'Path Length');

ax4 = nexttile(tl, 4);
plot_scalar_metric(ax4, cmp.risk_integrated, cmp.labels, cmp.colors, 'Integrated Risk', 'Integrated Risk');

ax5 = nexttile(tl, 5);
plot_scalar_metric(ax5, cmp.common_stealth_cost, cmp.labels, cmp.colors, 'Common Stealth Cost', 'Shared Cost');

ax6 = nexttile(tl, 6);
axis(ax6, 'off');
text(ax6, 0.02, 0.92, sprintf('DEM source: %s', planner_benchmark.environment.dem_source), 'FontSize', 11);
text(ax6, 0.02, 0.80, sprintf('Threat resolution: [%.1f, %.1f] m', ...
    planner_benchmark.environment.threat_resolution(1), planner_benchmark.environment.threat_resolution(2)), 'FontSize', 11);
text(ax6, 0.02, 0.68, sprintf('Setup time: %.2f s', planner_benchmark.shared_times.total_setup_s), 'FontSize', 11);
text(ax6, 0.02, 0.56, sprintf('Start risk: %.3f', planner_benchmark.environment.start_risk), 'FontSize', 11);
text(ax6, 0.02, 0.44, sprintf('Goal risk: %.3f', planner_benchmark.environment.goal_risk), 'FontSize', 11);
for i = 1:numel(cmp.labels)
    planner = planner_benchmark.planners.(cmp.keys{i});
    y = 0.30 - 0.12 * (i - 1);
    if planner.success
        status_text = sprintf('%s: success', cmp.labels{i});
    else
        status_text = sprintf('%s: failed', cmp.labels{i});
    end
    text(ax6, 0.02, y, status_text, 'FontSize', 11, 'Color', cmp.colors(i, :));
end
end

function plot_overlay_figure(planner_benchmark, tm, threat)
cmp = planner_benchmark.comparison;
fig = figure('Name', 'Planner Trajectory Overlay', ...
    'Position', [100, 90, 1200, 860], 'Color', 'w');
ax = axes('Parent', fig);
plot_terrain_surface(ax, tm, planner_benchmark.cfg.colors.terrain);
hold(ax, 'on');

for r = 1:numel(threat.radars)
    rr = threat.radars{r};
    draw_radar_sphere(ax, rr, planner_benchmark.cfg.colors.radar);
    plot3(ax, rr.position(1), rr.position(2), rr.position(3), 'kd', ...
        'MarkerFaceColor', planner_benchmark.cfg.colors.radar, 'MarkerSize', 8);
end

for i = 1:numel(cmp.keys)
    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success && ~isempty(planner.path_world)
        plot3(ax, planner.path_world(1, :), planner.path_world(2, :), planner.path_world(3, :), ...
            '-', 'Color', cmp.colors(i, :), 'LineWidth', 2.4, 'DisplayName', cmp.labels{i});
    end
end

plot3(ax, planner_benchmark.environment.start_world(1), planner_benchmark.environment.start_world(2), planner_benchmark.environment.start_world(3), ...
    'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9, 'DisplayName', 'Start');
plot3(ax, planner_benchmark.environment.goal_world(1), planner_benchmark.environment.goal_world(2), planner_benchmark.environment.goal_world(3), ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9, 'DisplayName', 'Goal');

style_3d_axes(ax, tm, planner_benchmark);
title(ax, '3D Terrain, Radars, and Planner Trajectories');
legend(ax, 'Location', 'northeastoutside');
end

function plot_planner_panels(planner_benchmark, tm, threat)
cmp = planner_benchmark.comparison;
fig = figure('Name', 'Per-Planner 3D Views', ...
    'Position', [110, 100, 1500, 560], 'Color', 'w');
tl = tiledlayout(fig, 1, numel(cmp.keys), 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(cmp.keys)
    ax = nexttile(tl, i);
    plot_terrain_surface(ax, tm, planner_benchmark.cfg.colors.terrain);
    hold(ax, 'on');

    for r = 1:numel(threat.radars)
        rr = threat.radars{r};
        draw_radar_sphere(ax, rr, planner_benchmark.cfg.colors.radar);
        plot3(ax, rr.position(1), rr.position(2), rr.position(3), 'kd', ...
            'MarkerFaceColor', planner_benchmark.cfg.colors.radar, 'MarkerSize', 7);
    end

    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success && ~isempty(planner.path_world)
        plot3(ax, planner.path_world(1, :), planner.path_world(2, :), planner.path_world(3, :), ...
            '-', 'Color', cmp.colors(i, :), 'LineWidth', 2.5);
    else
        text(ax, mean(tm.bounds(1:2)), mean(tm.bounds(3:4)), max(tm.Z(:)) + 60, ...
            'FAILED', 'Color', [0.7, 0.1, 0.1], 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    plot3(ax, planner_benchmark.environment.start_world(1), planner_benchmark.environment.start_world(2), planner_benchmark.environment.start_world(3), ...
        'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot3(ax, planner_benchmark.environment.goal_world(1), planner_benchmark.environment.goal_world(2), planner_benchmark.environment.goal_world(3), ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

    style_3d_axes(ax, tm, planner_benchmark);
    title(ax, cmp.labels{i});
end
end

function plot_scalar_metric(ax, values, labels, colors, y_label, ttl)
axes(ax);
bar_data = values;
h = bar(ax, bar_data, 0.65, 'FaceColor', 'flat');
for i = 1:numel(values)
    h.CData(i, :) = colors(i, :);
end
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels);
ylabel(ax, y_label);
title(ax, ttl);
grid(ax, 'on');

yl = ylim(ax);
for i = 1:numel(values)
    if isfinite(values(i))
        text(ax, i, values(i) + 0.03 * max(yl(2) - yl(1), 1), sprintf('%.2f', values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    else
        text(ax, i, yl(1) + 0.85 * (yl(2) - yl(1)), 'FAILED', ...
            'HorizontalAlignment', 'center', 'Color', [0.7, 0.1, 0.1], 'FontSize', 10, 'FontWeight', 'bold');
    end
end
end

function plot_stacked_time_metric(ax, time_stack_s, labels, colors, success)
axes(ax);
h = bar(ax, time_stack_s, 'stacked', 'BarWidth', 0.65);
h(1).FaceColor = [0.78, 0.83, 0.89];
h(2).FaceColor = [0.27, 0.48, 0.69];
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels);
ylabel(ax, 'Time [s]');
title(ax, 'Preprocess vs Search');
legend(ax, {'Preprocess', 'Search'}, 'Location', 'best');
grid(ax, 'on');

totals = sum(time_stack_s, 2, 'omitnan');
yl = ylim(ax);
for i = 1:numel(labels)
    if success(i) && isfinite(totals(i))
        text(ax, i, totals(i) + 0.03 * max(yl(2) - yl(1), 1), sprintf('%.2f', totals(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    else
        text(ax, i, yl(1) + 0.85 * (yl(2) - yl(1)), 'FAILED', ...
            'HorizontalAlignment', 'center', 'Color', [0.7, 0.1, 0.1], 'FontSize', 10, 'FontWeight', 'bold');
    end
end
end

function plot_terrain_surface(ax, tm, face_color)
[N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
surf(ax, N_grid, E_grid, tm.Z, ...
    'EdgeColor', 'none', ...
    'FaceColor', face_color, ...
    'FaceAlpha', 0.92);
colormap(ax, turbo);
end

function style_3d_axes(ax, tm, planner_benchmark)
all_alt = [planner_benchmark.environment.start_world(3), planner_benchmark.environment.goal_world(3)];
cmp = planner_benchmark.comparison;
for i = 1:numel(cmp.keys)
    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success && ~isempty(planner.path_world)
        all_alt = [all_alt, planner.path_world(3, :)]; %#ok<AGROW>
    end
end

terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));
alt_max = max([terrain_max + 100, all_alt + 50]);

xlim(ax, [tm.bounds(1), tm.bounds(2)]);
ylim(ax, [tm.bounds(3), tm.bounds(4)]);
zlim(ax, [terrain_min - 20, alt_max]);
xlabel(ax, 'North [m]');
ylabel(ax, 'East [m]');
zlabel(ax, 'Altitude [m]');
view(ax, 38, 30);
grid(ax, 'on');
axis(ax, 'tight');
end

function draw_radar_sphere(ax, radar_obj, color_rgb)
[sx, sy, sz] = sphere(24);
surf(ax, ...
    sx * radar_obj.R_max + radar_obj.position(1), ...
    sy * radar_obj.R_max + radar_obj.position(2), ...
    sz * radar_obj.R_max + radar_obj.position(3), ...
    'FaceColor', color_rgb, ...
    'FaceAlpha', 0.05, ...
    'EdgeColor', color_rgb, ...
    'EdgeAlpha', 0.08, ...
    'HandleVisibility', 'off');
end

function print_benchmark_summary(planner_benchmark)
cmp = planner_benchmark.comparison;
fprintf('\n=== Benchmark Summary ===\n');
fprintf('Shared setup: %.2f s (DEM %.2f s, threat %.2f s)\n', ...
    planner_benchmark.shared_times.total_setup_s, ...
    planner_benchmark.shared_times.dem_load_s, ...
    planner_benchmark.shared_times.threat_compute_s);

for i = 1:numel(cmp.keys)
    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success
        fprintf(['%-10s | total=%.2f s | preprocess=%.2f s | search=%.2f s | ', ...
                 'length=%.1f m | int_risk=%.3f | stealth=%.3f\n'], ...
            cmp.labels{i}, planner.total_time_s, planner.preprocess_time_s, planner.search_time_s, ...
            planner.path_length_m, planner.risk_integrated, planner.common_stealth_cost);
    else
        fprintf('%-10s | FAILED | total=%.2f s | preprocess=%.2f s | search=%.2f s | %s\n', ...
            cmp.labels{i}, planner.total_time_s, planner.preprocess_time_s, planner.search_time_s, planner.error_message);
    end
end

fprintf('planner_benchmark exported to workspace.\n');
end

function [NE, h, risk] = choose_endpoint(threat, tm, target_NE, agl, visibility_threshold)
offset_vals = -700:140:700;
offsets = zeros(numel(offset_vals)^2, 2);
cursor = 1;
for dn = offset_vals
    for de = offset_vals
        offsets(cursor, :) = [dn, de];
        cursor = cursor + 1;
    end
end

best_cost = inf;
best_NE = target_NE(:);
best_h = tm.get_height(best_NE(1), best_NE(2));
best_risk = 1.0;

for i = 1:size(offsets, 1)
    cand_NE = target_NE(:) + offsets(i, :)';
    cand_NE(1) = min(max(cand_NE(1), tm.bounds(1)), tm.bounds(2));
    cand_NE(2) = min(max(cand_NE(2), tm.bounds(3)), tm.bounds(4));
    cand_h = tm.get_height(cand_NE(1), cand_NE(2));
    if ~isfinite(cand_h)
        continue;
    end
    cand_alt = cand_h + agl;
    cand_risk = threat.get_risk(cand_NE(1), cand_NE(2), cand_alt);
    score = 1000 * max(0, cand_risk - visibility_threshold) + norm(offsets(i, :));
    if score < best_cost
        best_cost = score;
        best_NE = cand_NE;
        best_h = cand_h;
        best_risk = cand_risk;
    end
end

NE = best_NE;
h = best_h;
risk = best_risk;
end

function cost_val = compute_common_stealth_cost(path_ned, threat, alpha)
cost_val = NaN;
if isempty(path_ned) || size(path_ned, 2) < 2
    return;
end

cost_val = 0.0;
for k = 2:size(path_ned, 2)
    p0 = path_ned(:, k - 1);
    p1 = path_ned(:, k);
    seg_len = norm(p1 - p0);
    if seg_len < 1e-9
        continue;
    end
    n_samples = max(2, ceil(seg_len / 30));
    lambda = linspace(0, 1, n_samples);
    seg_cost = 0.0;
    for s = 1:n_samples
        ps = (1 - lambda(s)) * p0 + lambda(s) * p1;
        risk_s = threat.get_risk(ps(1), ps(2), -ps(3));
        seg_cost = seg_cost + exp(alpha * risk_s);
    end
    cost_val = cost_val + seg_len * seg_cost / n_samples;
end
end

function len = compute_path_length(path_xyz)
if isempty(path_xyz) || size(path_xyz, 2) < 2
    len = NaN;
    return;
end
len = sum(vecnorm(diff(path_xyz, 1, 2), 2, 1));
end

function x_ned = world_to_ned(x_world)
x_ned = [x_world(1); x_world(2); -x_world(3)];
end

function val = read_info_field(info, field_name)
val = NaN;
if isstruct(info) && isfield(info, field_name) && ~isempty(info.(field_name))
    val = info.(field_name);
end
end

function txt = read_info_text(info, field_name)
txt = '';
if isstruct(info) && isfield(info, field_name) && ~isempty(info.(field_name))
    raw = info.(field_name);
    if isstring(raw)
        txt = char(raw);
    elseif ischar(raw)
        txt = raw;
    end
end
end

function stop_reason = infer_rrt_stop_reason(info)
stop_reason = '';
if isfield(info, 'terminated_by_time_limit') && logical(info.terminated_by_time_limit)
    stop_reason = 'time_limit';
elseif isfield(info, 'success') && logical(info.success)
    stop_reason = 'goal_reached';
else
    stop_reason = 'max_iter_or_failure';
end
end

function msg = derive_failure_message(info, result)
msg = read_info_text(info, 'stop_reason');
if isempty(msg) && isfield(result, 'corridor_goal_connected') && islogical(result.corridor_goal_connected) && ~result.corridor_goal_connected
    msg = 'Strict corridor disconnected under the binary masking threshold.';
end
if isempty(msg) && isfield(info, 'terminated_by_time_limit') && logical(info.terminated_by_time_limit)
    msg = 'Planner hit the time limit before finding a path.';
end
if isempty(msg)
    msg = 'Planner finished without returning a feasible path.';
end
end

function summary = extract_native_summary(info, field_names)
summary = struct();
for i = 1:numel(field_names)
    name = field_names{i};
    if isfield(info, name) && ~isempty(info.(name))
        summary.(name) = info.(name);
    end
end
end

function val = getfield_with_default(s, field_name, default_val)
val = default_val;
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    val = s.(field_name);
end
end

function val = finite_or_nan(val)
if ~isfinite(val)
    val = NaN;
end
end

function val = finite_or_zero(val)
if ~isfinite(val)
    val = 0;
end
end
