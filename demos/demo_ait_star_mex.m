%% DEMO_AIT_STAR_MEX
% Direct AIT*-style radar-aware planning demo using C++ MEX.

clear; clc; close all;
clear ait_star_mex;
rehash;

this_dir = fileparts(mfilename('fullpath'));
% Add project and ait_star build folder to path
addpath(fullfile(fileparts(this_dir), 'project'));
addpath(fullfile(fileparts(this_dir), 'ait_star', 'build', 'Release'));
setup_project_paths();

fprintf('=== C++ AIT* MEX Demo ===\n');
fprintf('Using ait_star_mex: %s\n', which('ait_star_mex'));

if isempty(which('ait_star_mex'))
    error('ait_star_mex not found. Please ensure it is compiled and on the MATLAB path (ait_star/build/Release).');
end

cfg = struct();
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 30;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 2500;
cfg.visibility_threshold = 0.1;
cfg.start_goal_agl = 80;
cfg.threat_vert_res = 40;
cfg.use_parallel_threat = false;
cfg.use_gpu_threat = true;
cfg.radar_alt_offset = 20;
cfg.radar_range = 2000;
cfg.radar_range_model = 'cylinder';       % 'sphere' | 'cylinder'
cfg.radar_elevation_limits = [-90, 90];   % Wide vertical coverage for cylindrical masking
cfg.rng_seed = 7;
cfg.radar_position_mode = 'manual_fraction'; % 'manual_ne' | 'manual_fraction' | 'figure_click'
cfg.num_radars = 1;
cfg.radar_manual_ne = [];
cfg.radar_manual_fraction = [0.60; 0.20];
cfg.figure_click_use_contour = true;
cfg.figure_click_close_after = true;
cfg.skymap_gui = true;
cfg.skymap_mode = 'app';                  % 'app' | 'viewer'
cfg.skymap_half_window = 1800;            % [m] map window half-size around current waypoint
cfg.skymap_slice_half_length = 4000;      % [m] forward/backward span for vertical slice
cfg.skymap_slice_n_horiz = 181;           % Horizontal samples in vertical slice
cfg.skymap_slice_n_alt = 120;             % Altitude samples in vertical slice
cfg.skymap_play_period = 0.15;            % [s] auto-play update period
cfg.stealth_cost_alpha = 3.0;
cfg.fatrop_enable = true;
cfg.casadi_path = 'D:\casadi-3.7.2-windows64-matlab2018b';
cfg.fatrop_num_samples = 80;
cfg.fatrop_max_iter = 120;
cfg.fatrop_time_step = 0.25;
cfg.fatrop_tracking_weight = 0.04;
cfg.fatrop_smooth_weight = 1e-3;
cfg.fatrop_control_weight = 1e-4;
cfg.fatrop_radar_weight = 1.0;
cfg.fatrop_v_max = 25;
cfg.fatrop_a_max = 15;
cfg.fatrop_min_clearance = 10;
cfg.fatrop_max_clearance = 250;

rng(cfg.rng_seed, 'twister');

terrain_data = [];
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
    end
catch ME
    warning('%s', sprintf('DEM load failed (%s). Using synthetic terrain fallback.', ME.message));
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
    terrain_data.type = 'synthetic_ait_direct';
end

tm = terrain_map(terrain_data);
los = los_checker(tm);

N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));

threat_params = struct();
threat_params.resolution = [max(40, diag_len / 120), cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];

threat = threat_map(tm, los, 0.1, threat_params);
radar_ne = select_radar_positions_demo(tm, cfg, cfg.num_radars);

for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + cfg.radar_alt_offset;
    radar_obj = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar_obj.R_max = cfg.radar_range;
    radar_obj.range_model = cfg.radar_range_model;
    radar_obj.elevation_limits = cfg.radar_elevation_limits;
    radar_obj.P_t = 150e3;
    threat.add_radar(radar_obj);
end

fprintf('Computing binary threat map...\n');
t_threat = tic;
threat.compute_map('binary', struct( ...
    'use_gpu', cfg.use_gpu_threat, ...
    'use_parallel', cfg.use_parallel_threat, ...
    'show_progress', true));
t_threat_s = toc(t_threat);
fprintf('Threat map compute time: %.2f s\n', t_threat_s);

[start_NE, start_h, start_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(1) + 0.08 * N_span; tm.bounds(3) + 0.10 * E_span], cfg.start_goal_agl, cfg.visibility_threshold);
[goal_NE, goal_h, goal_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(2) - 0.08 * N_span; tm.bounds(4) - 0.10 * E_span], cfg.start_goal_agl, cfg.visibility_threshold);

start_world = [start_NE; start_h + cfg.start_goal_agl];
goal_world = [goal_NE; goal_h + cfg.start_goal_agl];
start_pos = [start_world(1:2); -start_world(3)];
goal_pos = [goal_world(1:2); -goal_world(3)];
fprintf('Start risk: %.2f | Goal risk: %.2f\n', start_risk, goal_risk);

planner_params = struct();
planner_params.max_batches = 200;
planner_params.batch_size = 300;
planner_params.max_nodes = 15000;
planner_params.stealth_weight = 30000.0;
planner_params.risk_limit = cfg.visibility_threshold;
planner_params.gamma = 20.0;
planner_params.preferred_agl = cfg.start_goal_agl;
planner_params.edge_check_samples = 9;
planner_params.post_solution_batches = 50;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.max_time = 120;
planner_params.min_clearance = 10;
planner_params.max_clearance = 250;

% Initialize MEX State
fprintf('\nInitializing MEX C++ Planner...\n');
ait_star_mex('init_map', double(tm.Z), double(tm.resolution), double([tm.bounds(1), tm.bounds(3)]), double(planner_params.min_clearance));
ait_star_mex('init_radar3d', single(threat.risk_grid), double(threat.N_vec), double(threat.E_vec), double(threat.alt_vec), 'asl', double(planner_params.stealth_weight));

fprintf('\nRunning C++ AIT* Planner...\n');
fprintf('---------------------------------\n');
fprintf('DEBUGGING MEX MEMORY:\n');
[t_elev, r_risk] = ait_star_mex('eval', double([start_world(1), start_world(2), start_world(3)]));
fprintf('Start C++ Elev: %.2f m, Risk: %.3f\n', t_elev, r_risk);
[t_elev, r_risk] = ait_star_mex('eval', double([goal_world(1), goal_world(2), goal_world(3)]));
fprintf('Goal C++ Elev: %.2f m, Risk: %.3f\n', t_elev, r_risk);
valid_count = 0;
for i = 1:100000
    px = -3000 + 6000 * rand();
    py = -3000 + 6000 * rand();
    pz = 2000 + 2800 * rand();
    [t, r] = ait_star_mex('eval', double([px, py, pz]));
    if pz >= t + 50 && r < 1.1
        valid_count = valid_count + 1;
    end
end
fprintf('Random Sample Valid Rate: %d / 100000\n', valid_count);
fprintf('---------------------------------\n');

% Run planner
t_total = tic;
res = ait_star_mex('plan', double(start_world), double(goal_world), planner_params);
fprintf('Total C++ planner:     %.2f s\n', res.planning_time_s);

path_world = res.path;
final_cost = res.cost;

if ~isempty(path_world)
    path_ned = [path_world(1, :); path_world(2, :); -path_world(3, :)];
    [path_max_risk, path_avg_risk, path_int_risk] = threat.evaluate_path(path_ned, 1);
else
    path_ned = [];
    path_max_risk = NaN;
    path_avg_risk = NaN;
    path_int_risk = NaN;
end

fprintf('\nRunning MinSnap smoothing and flatness verification...\n');
traj_ok = false;
traj_poly = struct();
Omega_ref = [];
omega_dot_ref = [];
path_smoothed = path_ned;
smoothed_cost = NaN;
if ~isempty(path_ned) && size(path_ned, 2) >= 2
    try
        path_simplified = simplify_path(path_ned, max(25, 0.005 * diag_len));
        path_len = sum(vecnorm(diff(path_simplified, 1, 2), 2, 1));

        smooth_params = struct();
        smooth_params.v_max = 25;
        smooth_params.a_max = 15;
        smooth_params.dt = 0.02;
        smooth_params.vel_bc_mode = 'free';
        smooth_params.cruise_speed = 15;
        smooth_params.max_waypoints = 50;
        smooth_params.max_seg_length = 500;
        smooth_params.aggressiveness = 4;

        total_time_guess = max(5, path_len / smooth_params.cruise_speed);
        traj_poly = trajectory_smoother(path_simplified, total_time_guess, smooth_params);
        traj_ok = true;

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
            q_inc = incremental_attitude_cmd_local(R_curr, tau_bz_c, psi_ref(k));
            q_ref = quat_normalize(quat_mul(q_ref, q_inc));
            R_ideal_next = quat_to_R(q_ref);

            [Omega_ref(:, k), omega_dot_ref(:, k)] = flatness( ...
                R_ideal_next, tau_norm, state_ref.j, traj_poly.snap(:, k), ...
                psi_dot_ref(k), psi_ddot_ref(k));
        end

        path_smoothed = traj_poly.pos;
        smoothed_cost = compute_stealth_path_cost_local(path_smoothed, threat, cfg.stealth_cost_alpha);
        fprintf('MinSnap trajectory: %d samples, duration %.2f s\n', n_samp, traj_poly.t(end));
        fprintf('Flatness omega_ref max norm: %.2f rad/s\n', max(vecnorm(Omega_ref, 2, 1)));
        fprintf('Flatness omega_dot_ref max norm: %.2f rad/s^2\n', max(vecnorm(omega_dot_ref, 2, 1)));
    catch ME
        warning('MinSnap/flatness stage failed: %s', ME.message);
    end
end

fprintf('\nRunning FATROP local trajectory optimization...\n');
fatrop_ok = false;
fatrop_info = struct('message', '', 'solve_time_s', NaN, 'setup_time_s', NaN, ...
    'solver_time_s', NaN, 'post_time_s', NaN, 'solver', 'fatrop');
path_fatrop_ned = [];
fatrop_cost = NaN;
if cfg.fatrop_enable && ~isempty(path_ned) && size(path_ned, 2) >= 2
    try
        t_fatrop = tic;
        [path_fatrop_ned, fatrop_info] = optimize_path_fatrop_local(path_ned, tm, threat, cfg);
        fatrop_info.solve_time_s = toc(t_fatrop);
        fatrop_cost = compute_stealth_path_cost_local(path_fatrop_ned, threat, cfg.stealth_cost_alpha);
        fatrop_ok = true;
        fprintf('FATROP trajectory: %d samples, total %.2f s | setup %.2f s | solver %.2f s | post %.2f s\n', ...
            size(path_fatrop_ned, 2), fatrop_info.solve_time_s, fatrop_info.setup_time_s, ...
            fatrop_info.solver_time_s, fatrop_info.post_time_s);
    catch ME
        fatrop_info.message = ME.message;
        warning('FATROP stage failed: %s', ME.message);
    end
end

fig_path = figure(1);
clf(fig_path);
set(fig_path, 'Name', 'MEX AIT* Path on Terrain', 'Position', [100, 80, 1180, 820], ...
    'Color', 'w', 'Visible', 'on', 'WindowStyle', 'normal');
tm.plot(fig_path);
hold on;
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    draw_radar_range(rr);
    plot3(rr.position(1), rr.position(2), rr.position(3), 'kd', ...
        'MarkerFaceColor', [1.0, 0.55, 0.1], 'MarkerSize', 9);
end
if ~isempty(path_world)
    plot3(path_world(1, :), path_world(2, :), path_world(3, :), '-', ...
        'Color', [0.92, 0.15, 0.18], 'LineWidth', 2.6);
end
if traj_ok
    plot3(path_smoothed(1, :), path_smoothed(2, :), -path_smoothed(3, :), '-', ...
        'Color', [0.08, 0.72, 0.92], 'LineWidth', 2.3);
end
if fatrop_ok
    plot3(path_fatrop_ned(1, :), path_fatrop_ned(2, :), -path_fatrop_ned(3, :), '-', ...
        'Color', [0.58, 0.16, 0.88], 'LineWidth', 2.3);
end
plot3(start_world(1), start_world(2), start_world(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot3(goal_world(1), goal_world(2), goal_world(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
title('Terrain Mesh With MEX C++ AIT* Path and Radar Coverage');
view(38, 30);
axis tight;
grid on;
legend({'Terrain', 'Radar range', 'Radar', 'Planned path', 'MinSnap trajectory', 'FATROP trajectory', 'Start', 'Goal'}, ...
    'Location', 'northeastoutside');

fig_diag = figure(2);
clf(fig_diag);
set(fig_diag, 'Name', 'MEX AIT* Diagnostics', 'Position', [130, 100, 600, 520], ...
    'Color', 'w', 'Visible', 'on', 'WindowStyle', 'normal');

imagesc(threat.N_vec, threat.E_vec, max(threat.risk_grid, [], 3));
set(gca, 'YDir', 'normal');
axis equal tight;
hold on;
if ~isempty(path_world)
    plot(path_world(1, :), path_world(2, :), '-', 'Color', [0.92, 0.15, 0.18], 'LineWidth', 2.2);
end
if traj_ok
    plot(path_smoothed(1, :), path_smoothed(2, :), '-', 'Color', [0.08, 0.72, 0.92], 'LineWidth', 2.0);
end
if fatrop_ok
    plot(path_fatrop_ned(1, :), path_fatrop_ned(2, :), '-', 'Color', [0.58, 0.16, 0.88], 'LineWidth', 2.0);
end
plot(start_world(1), start_world(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(goal_world(1), goal_world(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
title('Top-Down Threat Map and Final Path (C++ MEX)');
xlabel('North [m]');
ylabel('East [m]');
colormap(gca, turbo);
colorbar;

if cfg.skymap_gui && ~isempty(path_ned)
    try
        if isfield(cfg, 'skymap_mode') && strcmpi(cfg.skymap_mode, 'viewer')
            skymap_ui = skymap_viewer(cfg, tm, threat, path_ned, start_pos, goal_pos, cfg.visibility_threshold); %#ok<NASGU>
        else
            skymap_ui = skymap_app(cfg, tm, threat, path_ned, start_pos, goal_pos, cfg.visibility_threshold); %#ok<NASGU>
        end
    catch ME
        warning('SkyMap launch failed: %s', ME.message);
    end
end

if traj_ok || fatrop_ok
    fig_compare = figure(4);
    clf(fig_compare);
    set(fig_compare, 'Name', 'MinSnap vs FATROP Comparison', ...
        'Position', [180, 160, 820, 420], 'Color', 'w', ...
        'Visible', 'on', 'WindowStyle', 'normal');
    cost_vals = [final_cost, smoothed_cost, fatrop_cost];
    cost_vals(~isfinite(cost_vals)) = NaN;
    bar(cost_vals, 0.62, 'FaceColor', [0.22, 0.48, 0.74]);
    set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'AIT* raw', 'MinSnap', 'FATROP'});
    ylabel('Stealth-weighted path cost');
    title('Post-Processing Cost Comparison');
    grid on;
end

fprintf('\n=================================\n');
fprintf('        MEX PLANNER STATS        \n');
fprintf('=================================\n');
fprintf('Threat map prep:       %.2f s\n', t_threat_s);
fprintf('Total C++ planner:     %.2f s\n', res.planning_time_s);
fprintf('AIT* batches run:      %d\n', res.batches_run);
fprintf('Forward tree size:     %d\n', res.nodes_in_tree);
fprintf('Total sampled nodes:   %d\n', res.total_nodes);
fprintf('Final cost:            %.2f\n', final_cost);
fprintf('Forward edge checks:   %d\n', res.edge_checks);
fprintf('First solution batch:  %d\n', res.first_solution_batch);
fprintf('First solution time:   %.2f s\n', res.first_solution_time_s);
fprintf('Stop reason:           %s\n', res.stop_reason);
fprintf('MinSnap cost:          %.2f\n', smoothed_cost);
fprintf('FATROP cost:           %.2f\n', fatrop_cost);
fprintf('FATROP solve time:     %.2f s\n', fatrop_info.solve_time_s);
fprintf('FATROP setup time:     %.2f s\n', fatrop_info.setup_time_s);
fprintf('FATROP solver time:    %.2f s\n', fatrop_info.solver_time_s);
fprintf('FATROP post time:      %.2f s\n', fatrop_info.post_time_s);
fprintf('Path max risk:         %.3f\n', path_max_risk);
fprintf('Path avg risk:         %.3f\n', path_avg_risk);
fprintf('Path integrated risk:  %.3f\n', path_int_risk);
fprintf('RNG seed:              %d\n', cfg.rng_seed);
fprintf('=================================\n\n');

% Release persistent state memory in MEX
ait_star_mex('release');

function radar_ne = select_radar_positions_demo(tm, cfg, n_radars)
n_radars = max(1, round(n_radars));
mode_name = lower(string(cfg.radar_position_mode));

switch mode_name
    case "manual_ne"
        radar_ne = cfg.radar_manual_ne;
        if isempty(radar_ne)
            error('cfg.radar_manual_ne is empty. Provide [2 x M] [N;E] radar coordinates.');
        end
        if size(radar_ne, 1) ~= 2
            error('cfg.radar_manual_ne must be [2 x M] as [N;E].');
        end

    case "manual_fraction"
        frac = cfg.radar_manual_fraction;
        if isempty(frac)
            error('cfg.radar_manual_fraction is empty. Provide [2 x M] fractions in [0,1].');
        end
        if size(frac, 1) ~= 2
            error('cfg.radar_manual_fraction must be [2 x M] as [N_frac;E_frac].');
        end
        frac = min(max(double(frac), 0), 1);
        radar_ne = [tm.bounds(1) + frac(1, :) * (tm.bounds(2) - tm.bounds(1)); ...
                    tm.bounds(3) + frac(2, :) * (tm.bounds(4) - tm.bounds(3))];

    case "figure_click"
        fig = figure('Name', 'Select Radar Positions', 'Position', [120, 120, 1000, 760], 'Color', 'w');
        imagesc(tm.N_vec, tm.E_vec, tm.Z);
        set(gca, 'YDir', 'normal');
        axis equal tight;
        xlabel('North [m]');
        ylabel('East [m]');
        title(sprintf('Click %d radar position(s), then press Enter', n_radars));
        colormap(turbo);
        colorbar;
        hold on;
        if isfield(cfg, 'figure_click_use_contour') && cfg.figure_click_use_contour
            [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
            contour(N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
        end

        fprintf('Figure click mode: select %d radar position(s) on the terrain map.\n');
        fprintf('  Left click to add points, Enter to finish.\n');
        [x, y] = ginput(n_radars);
        if isempty(x)
            error('No radar point selected in figure_click mode.');
        end

        radar_ne = [x(:)'; y(:)'];
        if size(radar_ne, 2) < n_radars
            warning('Selected %d/%d points. Repeating last point for remaining radars.', ...
                size(radar_ne, 2), n_radars);
            radar_ne = [radar_ne, repmat(radar_ne(:, end), 1, n_radars - size(radar_ne, 2))];
        end
        radar_ne(1, :) = min(max(radar_ne(1, :), tm.bounds(1)), tm.bounds(2));
        radar_ne(2, :) = min(max(radar_ne(2, :), tm.bounds(3)), tm.bounds(4));

        plot(radar_ne(1, :), radar_ne(2, :), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        drawnow;
        if isfield(cfg, 'figure_click_close_after') && cfg.figure_click_close_after
            close(fig);
        end

    otherwise
        error('Unknown cfg.radar_position_mode: %s', cfg.radar_position_mode);
end

if size(radar_ne, 2) < n_radars
    radar_ne = [radar_ne, repmat(radar_ne(:, end), 1, n_radars - size(radar_ne, 2))];
elseif size(radar_ne, 2) > n_radars
    radar_ne = radar_ne(:, 1:n_radars);
end
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

function draw_radar_range(radar_obj)
[sx, sy, sz] = sphere(28);
surf( ...
    sx * radar_obj.R_max + radar_obj.position(1), ...
    sy * radar_obj.R_max + radar_obj.position(2), ...
    sz * radar_obj.R_max + radar_obj.position(3), ...
    'FaceColor', [1.0, 0.35, 0.15], ...
    'FaceAlpha', 0.05, ...
    'EdgeColor', [1.0, 0.45, 0.18], ...
    'EdgeAlpha', 0.10, ...
    'HandleVisibility', 'off');
end

function q_inc = incremental_attitude_cmd_local(R_curr, tau_bz_c, psi_ref)
b3_des = tau_bz_c(:);
b3_norm = norm(b3_des);
if b3_norm < 1e-9
    b3_des = R_curr(:, 3);
else
    b3_des = b3_des / b3_norm;
end

b1_yaw = [cos(psi_ref); sin(psi_ref); 0];
b2_des = cross(b3_des, b1_yaw);
if norm(b2_des) < 1e-9
    b1_yaw = [cos(psi_ref + pi / 2); sin(psi_ref + pi / 2); 0];
    b2_des = cross(b3_des, b1_yaw);
end
b2_des = b2_des / max(norm(b2_des), 1e-9);
b1_des = cross(b2_des, b3_des);
b1_des = b1_des / max(norm(b1_des), 1e-9);
R_des = [b1_des, b2_des, b3_des];

q_curr = R_to_quat(R_curr);
q_des = R_to_quat(R_des);
q_inc = quat_mul(quat_conj(q_curr), q_des);
q_inc = quat_normalize(q_inc);
end

function [path_opt_ned, info] = optimize_path_fatrop_local(path_ned, tm, threat, cfg)
info = struct('solver', 'fatrop', 'solve_time_s', NaN, 'setup_time_s', NaN, ...
    'solver_time_s', NaN, 'post_time_s', NaN, 'message', '');
setup_timer = tic;
if ~isempty(cfg.casadi_path) && exist(cfg.casadi_path, 'dir') == 7
    addpath(cfg.casadi_path);
end

try
    import casadi.*
catch ME
    error('CasADi is not available. Set cfg.casadi_path to a valid CasADi MATLAB folder. (%s)', ME.message);
end
if ~(exist('casadi.Opti', 'class') == 8 || exist('casadi.MX', 'class') == 8)
    error('CasADi classes were not found on the MATLAB path.');
end

N = max(6, round(cfg.fatrop_num_samples));
warm_ned = resample_path_arclength_local(path_ned, N);
warm_asl = [warm_ned(1, :); warm_ned(2, :); -warm_ned(3, :)];
terrain_warm = tm.get_height(warm_asl(1, :), warm_asl(2, :));
terrain_warm = reshape(terrain_warm, 1, []);
terrain_warm(~isfinite(terrain_warm)) = min(tm.Z(:));

path_len = sum(vecnorm(diff(warm_asl, 1, 2), 2, 1));
dt_min_speed = path_len / max(1e-6, 0.8 * cfg.fatrop_v_max * max(1, N - 1));
dt = max([1e-3, cfg.fatrop_time_step, dt_min_speed]);
opti = Opti();
P = opti.variable(3, N);
V = opti.variable(3, N);
A = opti.variable(3, N - 1);

opti.subject_to(P(:, 1) == warm_asl(:, 1));
opti.subject_to(P(:, end) == warm_asl(:, end));
for k = 1:(N - 1)
    opti.subject_to(P(:, k + 1) == P(:, k) + dt * V(:, k) + 0.5 * dt^2 * A(:, k));
    opti.subject_to(V(:, k + 1) == V(:, k) + dt * A(:, k));
    opti.subject_to(sumsqr(V(:, k)) <= cfg.fatrop_v_max^2);
    opti.subject_to(sumsqr(A(:, k)) <= cfg.fatrop_a_max^2);
end
opti.subject_to(sumsqr(V(:, N)) <= cfg.fatrop_v_max^2);
opti.subject_to(P(1, :) >= tm.bounds(1));
opti.subject_to(P(1, :) <= tm.bounds(2));
opti.subject_to(P(2, :) >= tm.bounds(3));
opti.subject_to(P(2, :) <= tm.bounds(4));
opti.subject_to(P(3, :) >= terrain_warm + cfg.fatrop_min_clearance);
opti.subject_to(P(3, :) <= terrain_warm + cfg.fatrop_max_clearance);

V0 = finite_difference_velocity_local(warm_asl, dt);
A0 = finite_difference_acceleration_local(V0, dt);
opti.set_initial(P, warm_asl);
opti.set_initial(V, V0);
opti.set_initial(A, A0);

J = cfg.fatrop_tracking_weight * sumsqr(P - warm_asl) + cfg.fatrop_control_weight * sumsqr(A);
if N > 2
    J = J + cfg.fatrop_smooth_weight * sumsqr(A(:, 2:end) - A(:, 1:end-1));
end
for k = 1:N
    risk_k = threat.get_risk(warm_asl(1, k), warm_asl(2, k), warm_asl(3, k));
    if isfinite(risk_k)
        for r = 1:length(threat.radars)
            radar_pos = threat.radars{r}.position(:);
            range2 = sumsqr(P(:, k) - radar_pos) + 1.0;
            J = J + cfg.fatrop_radar_weight * exp(cfg.stealth_cost_alpha * risk_k) / range2;
        end
    end
end
opti.minimize(J);

opts = struct();
opts.print_time = 0;
opts.fatrop.max_iter = cfg.fatrop_max_iter;
info.setup_time_s = toc(setup_timer);
try
    solver_timer = tic;
    evalc('opti.solver(''fatrop'', opts); sol = opti.solve();');
    info.solver_time_s = toc(solver_timer);
catch ME
    error('FATROP solver failed or is unavailable in this CasADi install: %s', ME.message);
end

post_timer = tic;
P_sol = full(sol.value(P));
P_sol(1, :) = min(max(P_sol(1, :), tm.bounds(1)), tm.bounds(2));
P_sol(2, :) = min(max(P_sol(2, :), tm.bounds(3)), tm.bounds(4));
terrain_sol = tm.get_height(P_sol(1, :), P_sol(2, :));
P_sol(3, :) = max(P_sol(3, :), terrain_sol + cfg.fatrop_min_clearance);
P_sol(3, :) = min(P_sol(3, :), terrain_sol + cfg.fatrop_max_clearance);
path_opt_ned = [P_sol(1, :); P_sol(2, :); -P_sol(3, :)];
info.post_time_s = toc(post_timer);
end

function cost_val = compute_stealth_path_cost_local(path_ned, threat, alpha)
cost_val = Inf;
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

function path_out = resample_path_arclength_local(path_in, n_samples)
if size(path_in, 2) < 2
    path_out = repmat(path_in(:, 1), 1, n_samples);
    return;
end
d = [0, cumsum(vecnorm(diff(path_in, 1, 2), 2, 1))];
if d(end) < 1e-9
    path_out = repmat(path_in(:, 1), 1, n_samples);
    return;
end
s = linspace(0, d(end), n_samples);
path_out = zeros(size(path_in, 1), n_samples);
for i = 1:size(path_in, 1)
    path_out(i, :) = interp1(d, path_in(i, :), s, 'linear');
end
end

function V = finite_difference_velocity_local(P, dt)
V = zeros(size(P));
if size(P, 2) == 1
    return;
end
V(:, 1:end-1) = diff(P, 1, 2) ./ dt;
V(:, end) = V(:, end-1);
end

function A = finite_difference_acceleration_local(V, dt)
A = zeros(size(V, 1), max(0, size(V, 2) - 1));
if isempty(A)
    return;
end
A(:, :) = diff(V, 1, 2) ./ dt;
end
