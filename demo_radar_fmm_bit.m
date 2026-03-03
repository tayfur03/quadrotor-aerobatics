%% DEMO_RADAR_FMM_BIT
% FMM-guided BIT* benchmark for deterministic binary radar masking.
%
% This script reuses the terrain/radar/threat setup style from
% demo_radar_binary_masking.m, then compares:
%   1) FMM-guided RRT* (corridor-bounded search)
%   2) FMM-guided BIT* (corridor-union informed sampler)
%
% Certification rationale:
% The "geometric tube" initializer (skeleton + adaptive radii) is fully
% deterministic, interpretable, and auditable. This is preferable to
% black-box GAN initializers when flight-safety traceability is required.

clear; clc; close all;
clear bit_star_planner compute_stealth_corridor;
rehash;

addpath('terrain');
addpath('radar');
addpath('motion_planner');

fprintf('=== FMM Corridor + RRT* + BIT* Radar Demo ===\n');
fprintf('Using bit_star_planner: %s\n', which('bit_star_planner'));
fprintf('Using compute_stealth_corridor: %s\n', which('compute_stealth_corridor'));

%% 0) Configuration (aligned with demo_radar_binary_masking.m workflow)
cfg = struct();
cfg.profile = 'FAST';   % 'FAST' | 'ACCURATE'
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 50;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 4000;
cfg.use_mesh_los_threat = false;
cfg.use_mesh_los_verify = true;
cfg.mesh_los_eps = 0.75;
cfg.threat_horiz_res = [];
cfg.threat_vert_res = 40;
cfg.use_parallel_threat = true;
cfg.skymap_gui = true;                    % Sliding SkyMap-style GUI
cfg.skymap_mode = 'app';                  % 'app' | 'viewer'

cfg.radar_position_mode = 'manual_ne';  % 'manual_ne' | 'manual_fraction' | 'figure_click'
cfg.num_radars = 1;
cfg.radar_alt_offset = [];
cfg.radar_manual_ne = [-600; 1500];
cfg.radar_manual_fraction = [0.50; 0.50];
cfg.figure_click_use_contour = true;
cfg.figure_click_close_after = true;
cfg.radar_min_separation = 800;
cfg.radar_edge_margin_frac = 0.08;

cfg.height_query_mode = 'LUT';  % 'LUT' | 'Mesh'
cfg.height_lut_resolution = [];

cfg.visibility_threshold = 0.5;
cfg.start_agl = 40;
cfg.goal_agl = 40;
cfg.stealth_cost_alpha = 2.0;  % exp(alpha * risk) path integral
cfg.run_fmm_rrt_trial = false;
cfg.bit_force_full_budget = true;  % true => keep optimizing until max_batches/max_time
cfg.bit_strict_no_radar_violation = true; % true => BIT* uses stealth_weight = inf
cfg.bit_strict_zero_risk_tol = 1e-6;      % strict zero-risk threshold

if isempty(which('bwdist'))
    error(['Image Processing Toolbox function bwdist is missing. ', ...
           'compute_stealth_corridor requires bwdist for adaptive radii.']);
end

profile_name = upper(char(cfg.profile));
switch profile_name
    case 'FAST'
        cfg.threat_horiz_res = [];
        cfg.threat_vert_res = 40;
        cfg.corridor_alpha = 2.2;
        cfg.corridor_lambda = 1.8;
        cfg.corridor_r_min = 300;
        cfg.corridor_r_max = 800;
        cfg.corridor_max_sweeps = 56;

        bit_tuning = struct();
        bit_tuning.max_batches = 150;
        bit_tuning.batch_size = 300;
        bit_tuning.max_nodes = 20000;
        bit_tuning.eps_global = 0.1;
        bit_tuning.stealth_weight = 3;
        bit_tuning.edge_check_samples = 9;
        bit_tuning.post_solution_batches = 10;
        bit_tuning.eta_radius = 1.8;

        rrt_tuning = struct();
        rrt_tuning.max_iter = 10000;
        rrt_tuning.max_planning_time = 45;   % [s]
        rrt_tuning.goal_bias = 0.20;
        rrt_tuning.shadow_bias = 0.92;
        rrt_tuning.rewire_mode = 'knearest';
        rrt_tuning.rewire_k_min = 24;
        rrt_tuning.rewire_k_max = 96;
        rrt_tuning.edge_sample_spacing = 20;

    case 'ACCURATE'
        % Memory-safe accurate profile: higher fidelity than FAST without
        % exploding voxel count on large DEMs.
        cfg.threat_horiz_res = 30;
        cfg.threat_vert_res = 30;
        cfg.corridor_alpha = 2.6;
        cfg.corridor_lambda = 1.4;
        cfg.corridor_r_min = 70;
        cfg.corridor_r_max = 300;
        cfg.corridor_max_sweeps = 72;

        bit_tuning = struct();
        bit_tuning.max_batches = 110;
        bit_tuning.batch_size = 300;
        bit_tuning.max_nodes = 24000;
        bit_tuning.eps_global = 0.04;
        bit_tuning.stealth_weight = 3.2;
        bit_tuning.edge_check_samples = 13;
        bit_tuning.post_solution_batches = 14;
        bit_tuning.eta_radius = 1.7;

        rrt_tuning = struct();
        rrt_tuning.max_iter = 5000;
        rrt_tuning.max_planning_time = 75;   % [s]
        rrt_tuning.goal_bias = 0.22;
        rrt_tuning.shadow_bias = 0.94;
        rrt_tuning.rewire_mode = 'knearest';
        rrt_tuning.rewire_k_min = 32;
        rrt_tuning.rewire_k_max = 128;
        rrt_tuning.edge_sample_spacing = 16;

    otherwise
        error('Unknown cfg.profile "%s". Use ''FAST'' or ''ACCURATE''.', cfg.profile);
end
fprintf('Planner profile: %s\n', profile_name);

%% 1) Terrain and LOS setup
use_dem = true;
if ~exist(cfg.dem_file, 'file')
    warning('DEM file not found: %s. Falling back to synthetic terrain.', cfg.dem_file);
    use_dem = false;
end

if use_dem
    dem_params = struct();
    dem_params.target_resolution = cfg.dem_target_resolution;
    dem_params.fill_nodata = cfg.dem_fill_nodata;
    if ~isempty(cfg.dem_crop_half_size) && cfg.dem_crop_half_size > 0
        h = cfg.dem_crop_half_size;
        dem_params.crop_bounds = [-h, h, -h, h];
    end
    terrain_data = dem_loader(cfg.dem_file, dem_params);
else
    terrain_params = struct();
    terrain_params.bounds = [0, 1200, -600, 600];
    terrain_params.resolution = 10;
    terrain_params.type = 'mountain';
    terrain_params.amplitude = 260;
    terrain_params.wavelength = 380;
    terrain_params.base_height = 5;
    terrain_params.seed = 42;
    terrain_data = terrain_generator(terrain_params);
end

tm = terrain_map(terrain_data);
mesh = terrain_mesh(terrain_data);
los_threat = los_checker(tm, struct('use_mesh', cfg.use_mesh_los_threat, 'mesh', mesh, 'los_eps', cfg.mesh_los_eps));
los_verify = los_checker(tm, struct('use_mesh', cfg.use_mesh_los_verify, 'mesh', mesh, 'los_eps', cfg.mesh_los_eps)); %#ok<NASGU>

%% 2) Radar and binary threat map
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));
center_E = mean(tm.bounds(3:4));

if isempty(cfg.radar_alt_offset)
    radar_alt_offset = max(10, 0.02 * (terrain_max - terrain_min + 1));
else
    radar_alt_offset = cfg.radar_alt_offset;
end

threat_params = struct();
if isempty(cfg.threat_horiz_res)
    threat_horiz_res = max(35, diag_len / 120);
else
    threat_horiz_res = cfg.threat_horiz_res;
end
threat_params.resolution = [threat_horiz_res, cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + 1, terrain_max + 300];

threat = threat_map(tm, los_threat, 0.1, threat_params);
radar_ne = select_radar_positions_demo(tm, cfg, cfg.num_radars);
radar_range = max(500, 0.35 * diag_len);

for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + radar_alt_offset;
    radar_obj = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar_obj.R_max = radar_range;
    radar_obj.P_t = 250e3;
    threat.add_radar(radar_obj);
end

fprintf('Computing deterministic visibility map...\n');
tic;
threat.compute_map('binary', struct('use_parallel', cfg.use_parallel_threat, 'show_progress', true));
t_threat = toc;
fprintf('Threat map compute time: %.2f s\n', t_threat);

%% 3) Mission endpoints
start_N_fixed = tm.bounds(1) + 0.05 * N_span;
goal_N_fixed = tm.bounds(2) - 0.05 * N_span;

[start_NE, start_h, start_risk] = find_hidden_endpoint_demo( ...
    tm, threat, start_N_fixed, center_E, cfg.start_agl, cfg.visibility_threshold, E_span);
[goal_NE, goal_h, goal_risk] = find_hidden_endpoint_demo( ...
    tm, threat, goal_N_fixed, center_E, cfg.goal_agl, cfg.visibility_threshold, E_span);

fprintf('Start risk (binary): %.2f | Goal risk (binary): %.2f\n', start_risk, goal_risk);

start_pos = [start_NE; -(start_h + cfg.start_agl)]; % NED
goal_pos = [goal_NE; -(goal_h + cfg.goal_agl)];     % NED

start_world = [start_pos(1); start_pos(2); -start_pos(3)]; % [N;E;Alt]
goal_world = [goal_pos(1); goal_pos(2); -goal_pos(3)];

%% 4) Base planner params (same family as demo_radar_binary_masking.m)
planner_params = struct();
planner_params.max_iter = cfg.dem_crop_half_size * 6;
planner_params.base_step_size = max(40, min(75, diag_len / 85));
planner_params.max_step_size = planner_params.base_step_size * 3.5;
planner_params.rewire_radius = planner_params.base_step_size * 3;
planner_params.rewire_mode = 'fixed';
planner_params.rewire_k_const = 3.2 * exp(1) * (1 + 1/3);
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
planner_params.height_query_mode = lower(cfg.height_query_mode);
planner_params.height_lut_resolution = cfg.height_lut_resolution;
planner_params.height_mesh_mode = 'raycast';
planner_params.terrain_mesh = mesh;
planner_params.radar_hard_constraint = true;
planner_params.radar_visibility_threshold = cfg.visibility_threshold;

min_alt = terrain_min + planner_params.min_clearance;
max_alt = terrain_max + planner_params.max_flight_alt;
base_bounds = [tm.bounds(1), tm.bounds(2), tm.bounds(3), tm.bounds(4), -max_alt, -min_alt];
planner_params.bounds = base_bounds;

%% 5) Strategic initializer: deterministic stealth corridor
corr_params = struct();
corr_params.N_vec = threat.N_vec;
corr_params.E_vec = threat.E_vec;
corr_params.alt_vec = threat.alt_vec;
corr_params.max_sweeps = cfg.corridor_max_sweeps;
corr_params.lambda = cfg.corridor_lambda;
corr_params.r_min = cfg.corridor_r_min;
corr_params.r_max = cfg.corridor_r_max;
corr_params.max_backtrack_points = 8000;

fprintf('\nComputing FMM/FSM stealth corridor...\n');
t_corr = tic;
corridor = compute_stealth_corridor( ...
    threat.risk_grid, start_world, goal_world, cfg.visibility_threshold, cfg.corridor_alpha, corr_params);
fprintf('Corridor build time: %.2f s | skeleton points: %d\n', toc(t_corr), size(corridor.skeleton, 2));

safe_ratio = nnz(corridor.safe_mask) / numel(corridor.safe_mask);
goal_idx = round(corridor.grid.goal_idx);
T_goal = corridor.T(goal_idx(1), goal_idx(2), goal_idx(3));
fprintf('Safe voxel ratio: %.2f%%\n', 100 * safe_ratio);
if ~isfinite(T_goal)
    warning(['FMM could not connect start->goal under strict threshold. ', ...
             'Applying ACCURATE fallback (slightly relaxed corridor threshold)...']);
    corridor_vis_threshold = min(0.70, cfg.visibility_threshold + 0.15);
    corridor = compute_stealth_corridor( ...
        threat.risk_grid, start_world, goal_world, corridor_vis_threshold, cfg.corridor_alpha, corr_params);
    goal_idx = round(corridor.grid.goal_idx);
    T_goal = corridor.T(goal_idx(1), goal_idx(2), goal_idx(3));
    fprintf('Fallback corridor threshold: %.2f | Connected: %d\n', corridor_vis_threshold, isfinite(T_goal));
end

%% 6) Trials: FMM-guided RRT* and FMM-guided BIT*
path_rrt_guided = [];
info_rrt_guided = struct('success', false, 'first_solution_time', NaN, 'iterations', NaN, ...
    'terminated_by_time_limit', false, 'planning_time', NaN);

if cfg.run_fmm_rrt_trial
    fprintf('\n[Trial 1/2] FMM-guided RRT*...\n');
    guided_bounds_world = corridor.bounds_world;
    guided_bounds_ned = [ ...
        max(base_bounds(1), guided_bounds_world(1)), ...
        min(base_bounds(2), guided_bounds_world(2)), ...
        max(base_bounds(3), guided_bounds_world(3)), ...
        min(base_bounds(4), guided_bounds_world(4)), ...
        max(base_bounds(5), -guided_bounds_world(6)), ...
        min(base_bounds(6), -guided_bounds_world(5))];

    if guided_bounds_ned(1) >= guided_bounds_ned(2) || ...
       guided_bounds_ned(3) >= guided_bounds_ned(4) || ...
       guided_bounds_ned(5) >= guided_bounds_ned(6)
        warning('Guided NED bounds invalid. Falling back to base bounds for RRT* trial.');
        guided_bounds_ned = base_bounds;
    end

    planner_params_rrt = planner_params;
    planner_params_rrt.max_iter = min(planner_params.max_iter, rrt_tuning.max_iter);
    planner_params_rrt.max_planning_time = rrt_tuning.max_planning_time;
    planner_params_rrt.bounds = guided_bounds_ned;
    planner_params_rrt.goal_bias = rrt_tuning.goal_bias;
    planner_params_rrt.shadow_bias = rrt_tuning.shadow_bias;
    planner_params_rrt.rewire_mode = rrt_tuning.rewire_mode;
    planner_params_rrt.rewire_k_min = rrt_tuning.rewire_k_min;
    planner_params_rrt.rewire_k_max = rrt_tuning.rewire_k_max;
    planner_params_rrt.edge_sample_spacing = rrt_tuning.edge_sample_spacing;

    fprintf(['RRT* caps: max_iter=%d, max_time=%.1fs, rewire=%s, ', ...
             'k=[%d,%d], edge_ds=%.1f\n'], ...
             planner_params_rrt.max_iter, planner_params_rrt.max_planning_time, ...
             planner_params_rrt.rewire_mode, planner_params_rrt.rewire_k_min, ...
             planner_params_rrt.rewire_k_max, planner_params_rrt.edge_sample_spacing);

    t_rrt_trial = tic;
    try
        [path_rrt_guided, info_rrt_guided] = rrt_star_radar(start_pos, goal_pos, tm, threat, planner_params_rrt);
    catch ME_rrt
        warning('FMM-RRT* trial crashed: %s', ME_rrt.message);
        path_rrt_guided = [];
        info_rrt_guided = struct('success', false, ...
            'first_solution_time', NaN, ...
            'iterations', NaN, ...
            'planning_time', toc(t_rrt_trial), ...
            'terminated_by_time_limit', false);
    end
    fprintf('FMM-RRT* elapsed wall time: %.2f s\n', toc(t_rrt_trial));
end

if cfg.run_fmm_rrt_trial
    fprintf('\n[Trial 2/2] FMM-guided BIT*...\n');
else
    fprintf('\n[Trial 1/1] FMM-guided BIT*...\n');
end
bit_params = struct();
bit_params.max_batches = bit_tuning.max_batches;
bit_params.batch_size = bit_tuning.batch_size;
bit_params.max_nodes = bit_tuning.max_nodes;
bit_params.eps_global = bit_tuning.eps_global;            % Completeness-preserving global branch
if cfg.bit_strict_no_radar_violation
    bit_params.stealth_weight = inf;
    bit_params.strict_zero_risk_tol = cfg.bit_strict_zero_risk_tol;
else
    bit_params.stealth_weight = bit_tuning.stealth_weight;
end
bit_params.edge_check_samples = bit_tuning.edge_check_samples;
bit_params.post_solution_batches = bit_tuning.post_solution_batches;
bit_params.goal_radius = planner_params.base_step_size;
bit_params.min_connection_radius = 0.6 * planner_params.base_step_size;
bit_params.max_connection_radius = 4.0 * planner_params.base_step_size;
bit_params.eta_radius = bit_tuning.eta_radius;
if cfg.bit_force_full_budget
    bit_params.post_solution_batches = inf;
end

[path_bit_world, info_bit] = bit_star_planner(start_world, goal_world, corridor, bit_params);
if info_bit.success
    path_bit_ned = [path_bit_world(1, :); path_bit_world(2, :); -path_bit_world(3, :)];
else
    path_bit_ned = [];
end

%% 7) MinSnap smoothing + flatness check
fprintf('\nRunning MinSnap smoothing and flatness verification...\n');

path_simplified = simplify_path(path_bit_ned, max(25, 0.005 * diag_len));
path_len = 0;
for i = 2:size(path_simplified, 2)
    path_len = path_len + norm(path_simplified(:, i) - path_simplified(:, i-1));
end

smooth_params = struct();
smooth_params.v_max = 25;
smooth_params.a_max = 8;
smooth_params.dt = 0.02;
smooth_params.vel_bc_mode = 'free';
smooth_params.cruise_speed = 15;
smooth_params.max_waypoints = 50;
smooth_params.max_seg_length = 500;
smooth_params.aggressiveness = 5;

traj_ok = false;
traj_poly = struct();
Omega_ref = [];
omega_dot_ref = [];
try
    total_time_guess = max(5, path_len / smooth_params.cruise_speed);
    traj_poly = trajectory_smoother(path_simplified, total_time_guess, smooth_params);
    traj_ok = true;

    % Mission-style chain:
    %   RRT* -> trajectory_smoother (with internal time allocation) ->
    %   yaw_planner('tangent') -> compute_incremental_attitude_cmd -> flatness
    e3 = [0; 0; 1];
    g = 9.81;
    n_samp = size(traj_poly.pos, 2);
    Omega_ref = zeros(3, n_samp);
    omega_dot_ref = zeros(3, n_samp);
    psi_ref = zeros(1, n_samp);
    psi_dot_ref = zeros(1, n_samp);
    psi_ddot_ref = zeros(1, n_samp);

    q_ref = [1; 0; 0; 0];
    yaw_state = struct('last_psi', 0);

    for k = 1:n_samp
        state_ref = struct();
        state_ref.pos = traj_poly.pos(:, k);
        state_ref.v = traj_poly.vel(:, k);
        state_ref.a = traj_poly.acc(:, k);
        state_ref.j = traj_poly.jerk(:, k);

        [psi_ref(k), psi_dot_ref(k), psi_ddot_ref(k)] = ...
            yaw_planner(traj_poly.t(k), 'tangent', state_ref, yaw_state);
        yaw_state.last_psi = psi_ref(k);

        tau_bz_c = state_ref.a - g * e3;
        tau_norm = norm(tau_bz_c);
        if tau_norm < 1e-6
            tau_bz_c = [0; 0; -g];
            tau_norm = norm(tau_bz_c);
        end

        R_curr = quat_to_R(q_ref);
        q_inc = compute_incremental_attitude_cmd(R_curr, tau_bz_c, psi_ref(k));
        q_ref = quat_normalize(quat_mul(q_ref, q_inc));
        R_ideal_next = quat_to_R(q_ref);

        [Omega_ref(:, k), omega_dot_ref(:, k)] = flatness( ...
            R_ideal_next, tau_norm, state_ref.j, traj_poly.snap(:, k), ...
            psi_dot_ref(k), psi_ddot_ref(k));
    end

    fprintf('MinSnap trajectory: %d samples, duration %.2f s (time allocation: ON)\n', n_samp, traj_poly.t(end));
    fprintf('Flatness omega_ref max norm: %.2f rad/s\n', max(vecnorm(Omega_ref, 2, 1)));
    fprintf('Flatness omega_dot_ref max norm: %.2f rad/s^2\n', max(vecnorm(omega_dot_ref, 2, 1)));
catch ME
    warning('MinSnap/flatness stage failed: %s', ME.message);
end

%% 8) Metrics
n_methods = 1 + double(cfg.run_fmm_rrt_trial);
method_names = cell(1, n_methods);
success = false(1, n_methods);
time_first = nan(1, n_methods);
iter_count = nan(1, n_methods);
final_cost = inf(1, n_methods);
paths_ned = cell(1, n_methods);

idx_method = 1;
if cfg.run_fmm_rrt_trial
    method_names{idx_method} = 'FMM-RRT*';
    success(idx_method) = logical(info_rrt_guided.success);
    time_first(idx_method) = info_rrt_guided.first_solution_time;
    iter_count(idx_method) = info_rrt_guided.iterations;
    if info_rrt_guided.success && ~isempty(path_rrt_guided)
        final_cost(idx_method) = compute_stealth_path_cost(path_rrt_guided, threat, cfg.stealth_cost_alpha);
    end
    paths_ned{idx_method} = path_rrt_guided;
    idx_method = idx_method + 1;
end

method_names{idx_method} = 'FMM-BIT*';
success(idx_method) = logical(info_bit.success);
time_first(idx_method) = info_bit.first_solution_time;
iter_count(idx_method) = info_bit.iterations;
if info_bit.success && ~isempty(path_bit_ned)
    final_cost(idx_method) = compute_stealth_path_cost(path_bit_ned, threat, cfg.stealth_cost_alpha);
end
paths_ned{idx_method} = path_bit_ned;

fprintf('\n=== Trial Summary ===\n');
for i = 1:numel(method_names)
    if success(i)
        if strcmp(method_names{i}, 'FMM-BIT*') && isfield(info_bit, 'stop_reason')
            fprintf('%-12s | t_first = %7.3f s | stealth_cost = %10.3f | iter = %6d | stop=%s\n', ...
                method_names{i}, time_first(i), final_cost(i), round(iter_count(i)), info_bit.stop_reason);
        else
            fprintf('%-12s | t_first = %7.3f s | stealth_cost = %10.3f | iter = %6d\n', ...
                method_names{i}, time_first(i), final_cost(i), round(iter_count(i)));
        end
    else
        if strcmp(method_names{i}, 'FMM-RRT*') && isfield(info_rrt_guided, 'terminated_by_time_limit') && info_rrt_guided.terminated_by_time_limit
            fprintf('%-12s | FAILED (time limit)\n', method_names{i});
        else
            fprintf('%-12s | FAILED\n', method_names{i});
        end
    end
end

figure('Name', 'Planner Metrics: FMM-RRT* vs FMM-BIT*', ...
       'Position', [100, 80, 1320, 430], 'Color', 'w');
subplot(1, 3, 1);
bar(time_first, 0.6, 'FaceColor', [0.2 0.55 0.85]);
set(gca, 'XTick', 1:numel(method_names), 'XTickLabel', method_names);
ylabel('Time to First [s]');
title('First Solution Time');
grid on;
subplot(1, 3, 2);
bar(final_cost, 0.6, 'FaceColor', [0.2 0.75 0.45]);
set(gca, 'XTick', 1:numel(method_names), 'XTickLabel', method_names);
ylabel('Stealth Cost');
title('Final Path Cost');
grid on;
subplot(1, 3, 3);
bar(iter_count, 0.6, 'FaceColor', [0.85 0.45 0.2]);
set(gca, 'XTick', 1:numel(method_names), 'XTickLabel', method_names);
ylabel('Iterations');
title('Iteration Count');
grid on;

%% 8) Visualization: raw trajectory comparison (always visible)
figure('Name', 'Raw Trajectory Comparison: FMM-RRT* vs FMM-BIT*', ...
       'Position', [110, 90, 1180, 820], 'Color', 'w');
tm.plot();
hold on;
[sx_cmp, sy_cmp, sz_cmp] = sphere(20);
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    surf(sx_cmp * rr.R_max + rr.position(1), ...
         sy_cmp * rr.R_max + rr.position(2), ...
         sz_cmp * rr.R_max + rr.position(3), ...
         'FaceColor', [1.0, 0.2, 0.2], 'FaceAlpha', 0.02, ...
         'EdgeColor', [1.0, 0.25, 0.25], 'EdgeAlpha', 0.04);
end
plot3(corridor.skeleton(1, :), corridor.skeleton(2, :), corridor.skeleton(3, :), ...
    '-', 'Color', [0.05, 0.55, 0.75], 'LineWidth', 1.2);

if cfg.run_fmm_rrt_trial && ~isempty(path_rrt_guided)
    plot3(path_rrt_guided(1, :), path_rrt_guided(2, :), -path_rrt_guided(3, :), ...
        '-', 'Color', [0.00, 0.65, 1.00], 'LineWidth', 2.1, 'DisplayName', 'FMM-RRT* (raw)');
end
if ~isempty(path_bit_ned)
    plot3(path_bit_ned(1, :), path_bit_ned(2, :), -path_bit_ned(3, :), ...
        '-', 'Color', [1.00, 0.20, 0.90], 'LineWidth', 2.1, 'DisplayName', 'FMM-BIT* (raw)');
end
plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 9);
xlabel('North [m]');
ylabel('East [m]');
zlabel('Altitude [m]');
title('Raw Planner Trajectories');
grid on;
view(38, 30);
axis tight;
legend('Location', 'best');

%% 9) Visualization: best smoothed trajectory (bounded)
best_path_ned = [];
best_label = 'No successful path';
if any(success)
    cost_for_select = final_cost;
    cost_for_select(~success) = inf;
    [~, best_idx] = min(cost_for_select);
    best_path_ned = paths_ned{best_idx};
    best_label = method_names{best_idx};
end

if ~isempty(best_path_ned) && size(best_path_ned, 2) >= 2
    path_simplified = simplify_path(best_path_ned, max(25, 0.005 * diag_len));
    path_smoothed = smooth_path_geometric(path_simplified, tm, 2);
    path_smoothed = enforce_path_bounds_and_clearance(path_smoothed, base_bounds, tm, planner_params.min_clearance);
else
    path_smoothed = best_path_ned;
end

figure('Name', 'Terrain + Stealth Corridor + Final Trajectory', ...
       'Position', [120, 90, 1180, 820], 'Color', 'w');
tm.plot();
hold on;

% Radar coverage hint
[sx, sy, sz] = sphere(20);
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    surf(sx * rr.R_max + rr.position(1), ...
         sy * rr.R_max + rr.position(2), ...
         sz * rr.R_max + rr.position(3), ...
         'FaceColor', [1.0, 0.2, 0.2], 'FaceAlpha', 0.03, ...
         'EdgeColor', [1.0, 0.25, 0.25], 'EdgeAlpha', 0.06);
end

% Transparent strategic corridor (subsampled spheres for speed)
corr_stride = max(1, floor(size(corridor.skeleton, 2) / 36));
[csx, csy, csz] = sphere(12);
for k = 1:corr_stride:size(corridor.skeleton, 2)
    cN = corridor.skeleton(1, k);
    cE = corridor.skeleton(2, k);
    cA = corridor.skeleton(3, k);
    cr = corridor.radii(k);
    surf(csx * cr + cN, csy * cr + cE, csz * cr + cA, ...
        'FaceColor', [0.1, 0.8, 0.95], 'FaceAlpha', 0.06, ...
        'EdgeColor', [0.1, 0.8, 0.95], 'EdgeAlpha', 0.05);
end
plot3(corridor.skeleton(1, :), corridor.skeleton(2, :), corridor.skeleton(3, :), ...
    '-', 'Color', [0.05, 0.55, 0.75], 'LineWidth', 1.8);

if ~isempty(path_smoothed)
    plot3(path_smoothed(1, :), path_smoothed(2, :), -path_smoothed(3, :), ...
        'm-', 'LineWidth', 2.8);
end
plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 9);

xlabel('North [m]');
ylabel('East [m]');
zlabel('Altitude [m]');
title(sprintf('Best Final Trajectory: %s', best_label));
grid on;
view(38, 30);
axis tight;

if cfg.skymap_gui
    if isfield(cfg, 'skymap_mode') && strcmpi(cfg.skymap_mode, 'viewer')
        skymap_ui = skymap_viewer(cfg, tm, threat, path_bit_ned, start_pos, goal_pos, cfg.visibility_threshold); %#ok<NASGU>
    else
        skymap_ui = skymap_app(cfg, tm, threat, path_bit_ned, start_pos, goal_pos, cfg.visibility_threshold); %#ok<NASGU>
    end
end

if traj_ok
    figure('Name', 'Flatness Check', 'Position', [140, 140, 1000, 700]);
    subplot(2, 1, 1);
    plot(traj_poly.t, Omega_ref(1, :), 'r-', traj_poly.t, Omega_ref(2, :), 'g-', traj_poly.t, Omega_ref(3, :), 'b-', 'LineWidth', 1.3);
    grid on;
    ylabel('\Omega_{ref} [rad/s]');
    title('Differential Flatness Angular Rate Feedforward');
    legend('p', 'q', 'r', 'Location', 'best');

    subplot(2, 1, 2);
    plot(traj_poly.t, omega_dot_ref(1, :), 'r--', traj_poly.t, omega_dot_ref(2, :), 'g--', traj_poly.t, omega_dot_ref(3, :), 'b--', 'LineWidth', 1.3);
    grid on;
    ylabel('\dot{\Omega}_{ref} [rad/s^2]');
    xlabel('Time [s]');
    title('Differential Flatness Angular Acceleration Feedforward');
    legend('dp/dt', 'dq/dt', 'dr/dt', 'Location', 'best');
end

%% ---------------- Local helpers ----------------
function radar_ne = select_radar_positions_demo(tm, cfg, n_radars)
n_radars = max(1, round(n_radars));
mode_name = lower(cfg.radar_position_mode);

switch mode_name
    case 'manual_ne'
        radar_ne = cfg.radar_manual_ne;
        if size(radar_ne, 1) ~= 2
            error('cfg.radar_manual_ne must be [2 x M] as [N;E].');
        end

    case 'manual_fraction'
        frac = cfg.radar_manual_fraction;
        if size(frac, 1) ~= 2
            error('cfg.radar_manual_fraction must be [2 x M] in [0,1].');
        end
        frac(1, :) = min(max(frac(1, :), 0), 1);
        frac(2, :) = min(max(frac(2, :), 0), 1);
        N_vals = tm.bounds(1) + frac(1, :) * (tm.bounds(2) - tm.bounds(1));
        E_vals = tm.bounds(3) + frac(2, :) * (tm.bounds(4) - tm.bounds(3));
        radar_ne = [N_vals; E_vals];

    case 'figure_click'
        fig = figure('Name', 'Select Radar Positions', 'Position', [120, 120, 1000, 760], 'Color', 'w');
        imagesc(tm.N_vec, tm.E_vec, tm.Z);
        set(gca, 'YDir', 'normal');
        axis equal tight;
        xlabel('North [m]');
        ylabel('East [m]');
        title(sprintf('Click 1/%d radar position(s), Enter to finish', n_radars));
        colormap(turbo);
        colorbar;
        hold on;
        if isfield(cfg, 'figure_click_use_contour') && cfg.figure_click_use_contour
            [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
            contour(N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
        end

        fprintf('Figure click mode: select %d radar position(s) on the terrain map.\n', n_radars);
        fprintf('  Left click to add points, Enter to finish.\n');
        x = zeros(1, n_radars);
        y = zeros(1, n_radars);
        n_sel = 0;
        while n_sel < n_radars
            [xi, yi, btn] = ginput(1); %#ok<ASGLU>
            if isempty(xi)
                break;
            end
            n_sel = n_sel + 1;
            x(n_sel) = xi;
            y(n_sel) = yi;
            plot(xi, yi, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
            title(sprintf('Click %d/%d radar position(s), Enter to finish', min(n_sel + 1, n_radars), n_radars));
            fprintf('  Selected %d/%d: N=%.1f, E=%.1f\n', n_sel, n_radars, xi, yi);
            drawnow;
        end
        if n_sel == 0
            error('No radar point selected in figure_click mode.');
        end

        radar_ne = [x(1:n_sel); y(1:n_sel)];
        if size(radar_ne, 2) < n_radars
            warning('Selected %d/%d points. Repeating last point for remaining radars.', ...
                size(radar_ne, 2), n_radars);
            radar_ne = [radar_ne, repmat(radar_ne(:, end), 1, n_radars - size(radar_ne, 2))];
        end

        plot(radar_ne(1, :), radar_ne(2, :), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        drawnow;
        if isfield(cfg, 'figure_click_close_after') && cfg.figure_click_close_after
            close(fig);
        end

    case 'terrain_peaks'
        margin_n = cfg.radar_edge_margin_frac * (tm.bounds(2) - tm.bounds(1));
        margin_e = cfg.radar_edge_margin_frac * (tm.bounds(4) - tm.bounds(3));
        N_min = tm.bounds(1) + margin_n;
        N_max = tm.bounds(2) - margin_n;
        E_min = tm.bounds(3) + margin_e;
        E_max = tm.bounds(4) - margin_e;

        [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
        N_flat = N_grid(:);
        E_flat = E_grid(:);
        Z_flat = tm.Z(:);

        valid = (N_flat >= N_min) & (N_flat <= N_max) & (E_flat >= E_min) & (E_flat <= E_max);
        N_flat = N_flat(valid);
        E_flat = E_flat(valid);
        Z_flat = Z_flat(valid);

        [~, order] = sort(Z_flat, 'descend');
        n_candidates = numel(order);
        selected = zeros(2, n_radars);
        selected_count = 0;
        used = false(n_candidates, 1);

        for idx_order = 1:n_candidates
            idx = order(idx_order);
            cand = [N_flat(idx); E_flat(idx)];
            if selected_count == 0
                selected_count = 1;
                selected(:, selected_count) = cand;
                used(idx_order) = true;
            else
                d = vecnorm(selected(:, 1:selected_count) - cand, 2, 1);
                if all(d >= cfg.radar_min_separation)
                    selected_count = selected_count + 1;
                    selected(:, selected_count) = cand;
                    used(idx_order) = true;
                end
            end
            if selected_count >= n_radars
                break;
            end
        end

        if selected_count < n_radars
            for idx_order = 1:n_candidates
                if used(idx_order)
                    continue;
                end
                idx = order(idx_order);
                cand = [N_flat(idx); E_flat(idx)];
                if selected_count == 0
                    selected_count = 1;
                    selected(:, selected_count) = cand;
                else
                    d = vecnorm(selected(:, 1:selected_count) - cand, 2, 1);
                    if all(d > 1e-3)
                        selected_count = selected_count + 1;
                        selected(:, selected_count) = cand;
                    end
                end
                if selected_count >= n_radars
                    break;
                end
            end
        end

        if selected_count == 0
            selected_count = 1;
            selected(:, 1) = [mean(tm.bounds(1:2)); mean(tm.bounds(3:4))];
        end
        radar_ne = selected(:, 1:selected_count);

    otherwise
        error('Unknown cfg.radar_position_mode: %s', cfg.radar_position_mode);
end

if size(radar_ne, 2) < n_radars
    reps = n_radars - size(radar_ne, 2);
    radar_ne = [radar_ne, repmat(radar_ne(:, end), 1, reps)];
elseif size(radar_ne, 2) > n_radars
    radar_ne = radar_ne(:, 1:n_radars);
end

radar_ne(1, :) = min(max(radar_ne(1, :), tm.bounds(1)), tm.bounds(2));
radar_ne(2, :) = min(max(radar_ne(2, :), tm.bounds(3)), tm.bounds(4));
end

function [NE, h, risk] = find_hidden_endpoint_demo(tm, threat, N_fixed, E_center, agl, risk_threshold, E_span)
offset_fracs = [0, 0.08, -0.08, 0.16, -0.16, 0.24, -0.24, 0.32, -0.32, 0.40, -0.40];
E_min = tm.bounds(3);
E_max = tm.bounds(4);
margin = 0.02 * (E_max - E_min);

best_cost = inf;
best_NE = [N_fixed; E_center];
best_h = tm.get_height(N_fixed, E_center);
best_risk = 1.0;

for k = 1:length(offset_fracs)
    E_try = E_center + offset_fracs(k) * E_span;
    E_try = min(max(E_try, E_min + margin), E_max - margin);
    h_try = tm.get_height(N_fixed, E_try);
    if isnan(h_try)
        continue;
    end

    alt_try = h_try + agl;
    risk_try = threat.get_risk(N_fixed, E_try, alt_try);
    cost = 1000 * max(0, risk_try - risk_threshold) + abs(offset_fracs(k));
    if cost < best_cost
        best_cost = cost;
        best_NE = [N_fixed; E_try];
        best_h = h_try;
        best_risk = risk_try;
    end

    if risk_try < risk_threshold
        break;
    end
end

NE = best_NE;
h = best_h;
risk = best_risk;
end

function cost = compute_stealth_path_cost(path_ned, threat, alpha_cost)
if isempty(path_ned) || size(path_ned, 2) < 2
    cost = inf;
    return;
end

segment = diff(path_ned, 1, 2);
seg_len = vecnorm(segment, 2, 1);
n_seg = numel(seg_len);

sample_counts = max(2, ceil(seg_len / 20) + 1);
total_samples = sum(sample_counts);

N = zeros(1, total_samples);
E = zeros(1, total_samples);
A = zeros(1, total_samples);
W = zeros(1, total_samples);

cursor = 1;
for s = 1:n_seg
    ns = sample_counts(s);
    t = linspace(0, 1, ns);
    p0 = path_ned(:, s);
    p1 = path_ned(:, s + 1);
    seg_pts = p0 + (p1 - p0) * t;
    idx = cursor:(cursor + ns - 1);
    N(idx) = seg_pts(1, :);
    E(idx) = seg_pts(2, :);
    A(idx) = -seg_pts(3, :);
    W(idx) = seg_len(s) / max(ns - 1, 1);
    cursor = cursor + ns;
end

risk = threat.get_risk(N, E, A);
risk = reshape(risk, 1, []);
risk(~isfinite(risk)) = 1.0;

cost = sum(exp(alpha_cost .* risk) .* W);
end

function path_out = enforce_path_bounds_and_clearance(path_in, bounds_ned, tm, min_clearance)
% Ensure displayed smoothed paths remain within planner bounds and above terrain.
if isempty(path_in)
    path_out = path_in;
    return;
end

path_out = path_in;
n_pts = size(path_out, 2);
for i = 1:n_pts
    N = min(max(path_out(1, i), bounds_ned(1)), bounds_ned(2));
    E = min(max(path_out(2, i), bounds_ned(3)), bounds_ned(4));
    D = min(max(path_out(3, i), bounds_ned(5)), bounds_ned(6));

    terrain_h = tm.get_height(N, E);
    min_alt = terrain_h + min_clearance;
    alt = -D;
    if alt < min_alt
        alt = min_alt;
        D = -alt;
    end

    path_out(:, i) = [N; E; D];
end
end


function q_inc = compute_incremental_attitude_cmd(R_curr, tau_bz_c, psi_ref)
% Same attitude increment logic used in mission/INDI demos.

t_norm = norm(tau_bz_c);
if t_norm < 1e-6
    t_des_in = [0; 0; -1];
else
    t_des_in = tau_bz_c / t_norm;
end

bz_des_in = -t_des_in;
bz_des_body = R_curr' * bz_des_in;

current_z = [0; 0; 1];
cross_prod = cross(current_z, bz_des_body);
dot_prod = dot(current_z, bz_des_body);

if dot_prod < -0.9999
    q_tilt = [0; 1; 0; 0];
else
    s = sqrt(2 * (1 + dot_prod));
    q_tilt = [0.5 * s; cross_prod / s];
end
q_tilt = quat_normalize(q_tilt);

R_tilt = quat_to_R(q_tilt);
R_new_in = R_curr * R_tilt;
psi_curr = atan2(R_new_in(2, 1), R_new_in(1, 1));

psi_err = psi_ref - psi_curr;
while psi_err > pi,  psi_err = psi_err - 2*pi; end
while psi_err < -pi, psi_err = psi_err + 2*pi; end

q_yaw = [cos(psi_err/2); 0; 0; sin(psi_err/2)];
q_inc = quat_mul(q_tilt, q_yaw);
if q_inc(1) < 0
    q_inc = -q_inc;
end
end