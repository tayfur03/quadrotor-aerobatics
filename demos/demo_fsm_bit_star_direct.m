%% DEMO_FSM_BIT_STAR_DIRECT
% Direct FSM-guided BIT* demo without corridor computation.

clear; clc; close all;
clear fsm_bit_star_planner_direct compute_fsm_heuristic;
rehash;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
setup_project_paths();

fprintf('=== Direct FSM-BIT* Demo ===\n');
fprintf('Using fsm_bit_star_planner_direct: %s\n', which('fsm_bit_star_planner_direct'));
fprintf('Using compute_fsm_heuristic: %s\n', which('compute_fsm_heuristic'));

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
cfg.skymap_gui = true;
cfg.skymap_mode = 'app';
cfg.stealth_cost_alpha = 2.0;
cfg.rng_seed = 7;
cfg.strict_zero_risk = true;
cfg.fsm_max_sweeps = 96;
cfg.fsm_tol = 1e-4;

rng(cfg.rng_seed, 'twister');

%% 1) Terrain
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
    terrain_data.type = 'synthetic_fsm_bit_direct';
end

tm = terrain_map(terrain_data);
los = los_checker(tm);

%% 2) Radar / threat volume
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));

threat_params = struct();
threat_params.resolution = [max(40, diag_len / 120), cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];

threat = threat_map(tm, los, 0.1, threat_params);
radar_ne = [tm.bounds(1) + 0.22 * N_span, tm.bounds(1) + 0.78 * N_span; ...
            tm.bounds(3) + 0.25 * E_span, tm.bounds(3) + 0.72 * E_span];
for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + cfg.radar_alt_offset;
    radar_obj = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar_obj.R_max = cfg.radar_range;
    radar_obj.P_t = 150e3;
    threat.add_radar(radar_obj);
end

fprintf('Computing binary threat map...\n');
t_threat = tic;
threat.compute_map('binary', struct('use_parallel', cfg.use_parallel_threat, 'show_progress', false));
t_threat_s = toc(t_threat);
fprintf('Threat map compute time: %.2f s\n', t_threat_s);

radar_los_map = threat.risk_grid >= cfg.visibility_threshold;

%% 3) Start / goal
[start_NE, start_h, start_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(1) + 0.08 * N_span; tm.bounds(3) + 0.10 * E_span], cfg.start_goal_agl, cfg.visibility_threshold);
[goal_NE, goal_h, goal_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(2) - 0.08 * N_span; tm.bounds(4) - 0.10 * E_span], cfg.start_goal_agl, cfg.visibility_threshold);

start_world = [start_NE; start_h + cfg.start_goal_agl];
goal_world = [goal_NE; goal_h + cfg.start_goal_agl];
fprintf('Start risk: %.2f | Goal risk: %.2f\n', start_risk, goal_risk);

%% 4) Direct env and planner
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
if cfg.strict_zero_risk
    planner_params.stealth_weight = inf;
else
    planner_params.stealth_weight = 3.0;
end
planner_params.strict_zero_risk_tol = 1e-6;
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

t_total = tic;
[path_world, final_cost, info] = fsm_bit_star_planner_direct(start_world, goal_world, env, planner_params);
t_total_s = toc(t_total);
if ~isempty(path_world)
    path_ned = [path_world(1, :); path_world(2, :); -path_world(3, :)];
else
    path_ned = [];
end
start_pos = [start_world(1); start_world(2); -start_world(3)];
goal_pos = [goal_world(1); goal_world(2); -goal_world(3)];
raw_path_max_risk = NaN;
raw_path_avg_risk = NaN;
raw_path_integrated_risk = NaN;
if ~isempty(path_ned)
    [raw_path_max_risk, raw_path_avg_risk, raw_path_integrated_risk] = ...
        threat.evaluate_path(path_ned, 1);
end

%% 5) MinSnap smoothing + flatness check
fprintf('\nRunning MinSnap smoothing and flatness verification...\n');
traj_ok = false;
traj_poly = struct();
Omega_ref = [];
omega_dot_ref = [];
path_simplified = path_ned;
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

%% 6) Figures
fig_path = figure(1);
clf(fig_path);
set(fig_path, 'Name', 'Direct FSM-BIT* Path on Terrain', ...
    'Position', [100, 80, 1180, 820], 'Color', 'w', 'Visible', 'on', ...
    'WindowStyle', 'normal');
tm.plot(fig_path);
hold on;
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    draw_radar_range(rr);
    plot3(rr.position(1), rr.position(2), rr.position(3), 'kd', ...
        'MarkerFaceColor', [1.0, 0.55, 0.1], 'MarkerSize', 9);
end
if ~isempty(path_world)
    plot3(path_world(1, :), path_world(2, :), path_world(3, :), ...
        '-', 'Color', [0.92, 0.15, 0.18], 'LineWidth', 2.6);
end
if traj_ok
    plot3(path_smoothed(1, :), path_smoothed(2, :), -path_smoothed(3, :), ...
        '-', 'Color', [0.08, 0.72, 0.92], 'LineWidth', 2.3);
end
plot3(start_world(1), start_world(2), start_world(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot3(goal_world(1), goal_world(2), goal_world(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
title('Terrain Mesh With Direct FSM-BIT* Path and Radar Coverage');
view(38, 30);
axis tight;
grid on;
legend({'Terrain', 'Radar range', 'Radar', 'Planned path', 'MinSnap trajectory', 'Start', 'Goal'}, ...
    'Location', 'northeastoutside');
drawnow;
figure(fig_path);

fig_diag = figure(2);
clf(fig_diag);
set(fig_diag, 'Name', 'Direct FSM-BIT* Diagnostics', ...
    'Position', [130, 100, 1240, 500], 'Color', 'w', 'Visible', 'on', ...
    'WindowStyle', 'normal');
subplot(1, 2, 1);
if ~isempty(path_ned)
    imagesc(threat.N_vec, threat.E_vec, max(threat.risk_grid, [], 3));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    hold on;
    plot(path_world(1, :), path_world(2, :), '-', 'Color', [0.92, 0.15, 0.18], 'LineWidth', 2.2);
    if traj_ok
        plot(path_smoothed(1, :), path_smoothed(2, :), '-', 'Color', [0.08, 0.72, 0.92], 'LineWidth', 2.0);
    end
    for r = 1:length(threat.radars)
        rr = threat.radars{r};
        th = linspace(0, 2 * pi, 180);
        plot(rr.position(1) + rr.R_max * cos(th), rr.position(2) + rr.R_max * sin(th), ...
            '--', 'Color', [1.0, 0.45, 0.18], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        plot(rr.position(1), rr.position(2), 'kd', 'MarkerFaceColor', [1.0, 0.55, 0.1], 'MarkerSize', 8, ...
            'HandleVisibility', 'off');
    end
    plot(start_world(1), start_world(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(goal_world(1), goal_world(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    title('Top-Down Threat Map and Post-Processed Trajectory');
    xlabel('North [m]');
    ylabel('East [m]');
    colormap(gca, turbo);
    colorbar;
else
    text(0.1, 0.5, 'No path available for post-planner visualization.', 'FontSize', 12);
    axis off;
end

subplot(1, 2, 2);
warm_cost_plot = finite_or_nan(info.warmstart_cost);
final_cost_plot = finite_or_nan(final_cost);
smoothed_cost_plot = finite_or_nan(smoothed_cost);
bar([warm_cost_plot, final_cost_plot, smoothed_cost_plot], 0.6, 'FaceColor', [0.18, 0.55, 0.82]);
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Warm-start', 'Final', 'MinSnap'});
ylabel('Cost');
title('Warm-start, Final, and Post-Processed Cost');
grid on;
yl = ylim;
text(1.05, yl(2) * 0.85, sprintf('Volume ratio \\rho = %.2f', info.volume_ratio), 'FontSize', 11);
text(1.05, yl(2) * 0.72, sprintf('Raw path max risk = %.3f', raw_path_max_risk), 'FontSize', 11);

if info.fsm_used && isfinite(final_cost) && ~isempty(info.fsm_T_fwd) && ~isempty(info.fsm_T_bwd)
    fig_heur = figure(3);
    clf(fig_heur);
    set(fig_heur, 'Name', 'Direct FSM-BIT* Heuristic Region', ...
        'Position', [160, 120, 1240, 520], 'Color', 'w', 'Visible', 'on', ...
        'WindowStyle', 'normal');
    if ~isempty(path_world)
        slice_alt = median(path_world(3, :));
    else
        slice_alt = median([start_world(3), goal_world(3)]);
    end
    [~, k_slice] = min(abs(info.fsm_meta.grid_z - slice_alt));
    T_sum_slice = double(info.fsm_T_fwd(:, :, k_slice)) + double(info.fsm_T_bwd(:, :, k_slice));
    informed_slice = isfinite(T_sum_slice) & (T_sum_slice < final_cost);

    subplot(1, 2, 1);
    imagesc(info.fsm_meta.grid_x, info.fsm_meta.grid_y, T_sum_slice.');
    set(gca, 'YDir', 'normal');
    axis equal tight;
    hold on;
    if ~isempty(path_world)
        plot(path_world(1, :), path_world(2, :), 'w-', 'LineWidth', 2.2);
    end
    plot(start_world(1), start_world(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(goal_world(1), goal_world(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    title(sprintf('T_{fwd} + T_{bwd} Slice at %.1f m', info.fsm_meta.grid_z(k_slice)));
    xlabel('North [m]');
    ylabel('East [m]');
    colormap(gca, turbo);
    colorbar;

    subplot(1, 2, 2);
    imagesc(info.fsm_meta.grid_x, info.fsm_meta.grid_y, informed_slice.');
    set(gca, 'YDir', 'normal');
    axis equal tight;
    hold on;
    if ~isempty(path_world)
        plot(path_world(1, :), path_world(2, :), '-', 'Color', [0.92, 0.15, 0.18], 'LineWidth', 2.2);
    end
    plot(start_world(1), start_world(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(goal_world(1), goal_world(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    title(sprintf('Heuristic Region: T_{sum} < %.1f', final_cost));
    xlabel('North [m]');
    ylabel('East [m]');
    colormap(gca, parula(2));
    colorbar('Ticks', [0, 1], 'TickLabels', {'Outside', 'Inside'});
end

if traj_ok
    fig_flat = figure(4);
    clf(fig_flat);
    set(fig_flat, 'Name', 'Flatness Check', 'Position', [140, 140, 1000, 700], ...
        'Color', 'w', 'Visible', 'on', 'WindowStyle', 'normal');
    subplot(2, 1, 1);
    plot(traj_poly.t, Omega_ref(1, :), 'r-', traj_poly.t, Omega_ref(2, :), 'g-', traj_poly.t, Omega_ref(3, :), 'b-', 'LineWidth', 1.3);
    grid on;
    ylabel('Omega_ref [rad/s]', 'Interpreter', 'none');
    title('Differential Flatness Angular Rate Feedforward');
    legend('p', 'q', 'r', 'Location', 'best');

    subplot(2, 1, 2);
    plot(traj_poly.t, omega_dot_ref(1, :), 'r--', traj_poly.t, omega_dot_ref(2, :), 'g--', traj_poly.t, omega_dot_ref(3, :), 'b--', 'LineWidth', 1.3);
    grid on;
    ylabel('dOmega_ref/dt [rad/s^2]', 'Interpreter', 'none');
    xlabel('Time [s]');
    title('Differential Flatness Angular Acceleration Feedforward');
    legend('dp/dt', 'dq/dt', 'dr/dt', 'Location', 'best');
end

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

%% 7) Console summary
fprintf('\nThreat map:        %.2f s\n', t_threat_s);
fprintf('FSM preprocessing: %.2f s\n', info.fsm_time_s);
fprintf('BIT* search:       %.2f s\n', info.planning_time);
fprintf('Total planner:     %.2f s\n', t_total_s);
fprintf('BIT* iterations:   %d\n', info.iterations);
fprintf('Warm-start cost:   %.2f\n', info.warmstart_cost);
fprintf('Final cost:        %.2f\n', final_cost);
fprintf('MinSnap cost:      %.2f\n', smoothed_cost);
fprintf('Raw path max risk: %.3f\n', raw_path_max_risk);
fprintf('Raw path avg risk: %.3f\n', raw_path_avg_risk);
fprintf('Raw path int risk: %.3f\n', raw_path_integrated_risk);
fprintf('Volume ratio rho:  %.2f\n', info.volume_ratio);
fprintf('FSM used:          %s\n', string(info.fsm_used));
if isfield(info, 'fsm_meta') && ~isempty(info.fsm_meta)
    fprintf('FSM sweeps fwd/bwd:%d / %d\n', read_meta_field(info.fsm_meta, 'sweeps_fwd'), read_meta_field(info.fsm_meta, 'sweeps_bwd'));
    fprintf('FSM delta fwd/bwd: %.3e / %.3e\n', read_meta_field(info.fsm_meta, 'final_delta_fwd'), read_meta_field(info.fsm_meta, 'final_delta_bwd'));
end
fprintf('Strict zero-risk:  %s\n', string(cfg.strict_zero_risk));
fprintf('RNG seed:          %d\n', cfg.rng_seed);

%% ---------------- Local helpers ----------------
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

function val = finite_or_nan(val)
if ~isfinite(val)
    val = NaN;
end
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

function val = read_meta_field(meta, field_name)
val = NaN;
if isstruct(meta) && isfield(meta, field_name) && ~isempty(meta.(field_name))
    val = meta.(field_name);
end
end
