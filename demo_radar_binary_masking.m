%% DEMO_RADAR_BINARY_MASKING - Deterministic radar visibility with terrain masking
%
% Goal:
%   Use a hard visibility model instead of probabilistic detection:
%   - Risk = 1 if target is visible to any radar (in range + LOS)
%   - Risk = 0 if terrain masked or out of coverage
%
% RRT* is configured with a hard radar constraint, so visible edges are
% rejected. This allows sampling/planning in radar shadow regions.

clear; clc; close all;

addpath('terrain');
addpath('radar');
addpath('motion_planner');

fprintf('=== Binary Radar Masking Demo ===\n');

%% 0) Configuration
cfg = struct();
cfg.dem_file = 'DEM/artvin.tif';          % Any local DEM/*.tif
cfg.dem_target_resolution = 50;           % [m] downsample for speed
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 4000;            % Crop around DEM center [m]
cfg.use_mesh_los_threat = false;          % Faster for dense threat-map computation
cfg.use_mesh_los_verify = true;           % Keep precise LOS for final path checks
cfg.mesh_los_eps = 0.75;
cfg.threat_horiz_res = [];                % [] => auto
cfg.threat_vert_res = 40;                 % Larger = fewer altitude layers (faster)
cfg.use_parallel_threat = true;           % Use parfor in threat.compute_map if available
cfg.radar_position_mode = 'figure_click';% 'terrain_peaks' | 'manual_ne' | 'manual_fraction' | 'figure_click'
cfg.num_radars = 1;
cfg.radar_min_separation = 800;           % [m] min separation between selected radar sites
cfg.radar_edge_margin_frac = 0.08;        % Ignore edge band when auto-selecting peaks
cfg.radar_alt_offset = [];                % [] => auto offset above terrain
cfg.radar_manual_ne = [0; 0];             % [2 x M] [N;E] if mode='manual_ne'
cfg.radar_manual_fraction = [0.50; 0.50]; % [2 x M] fraction in [0,1] if mode='manual_fraction'
cfg.figure_click_use_contour = true;      % Overlay contours in click figure
cfg.figure_click_close_after = true;      % Auto-close click figure after selection
cfg.rrt_live_animation = false;            % Live RRT* animation on 3D terrain
cfg.rrt_live_plot_interval = 500;         % Iteration interval for animation refresh
cfg.skymap_gui = true;                    % Sliding SkyMap-style GUI
cfg.skymap_mode = 'app';                  % 'app' | 'viewer'
cfg.skymap_half_window = 1800;            % [m] map window half-size around current waypoint
cfg.skymap_slice_half_length = 2500;      % [m] forward/backward span for vertical slice
cfg.skymap_slice_n_horiz = 181;           % Horizontal samples in vertical slice
cfg.skymap_slice_n_alt = 120;             % Altitude samples in vertical slice
cfg.skymap_play_period = 0.15;            % [s] auto-play update period

%% 1) Terrain (DEM or synthetic)
use_dem = true;
if ~exist(cfg.dem_file, 'file')
    warning('DEM file not found: %s. Falling back to synthetic terrain.', cfg.dem_file);
    use_dem = false;
end
terrain_label = 'SYNTHETIC';
if use_dem
    fprintf('Terrain source: DEM (%s)\n', cfg.dem_file);
    terrain_label = sprintf('DEM: %s', cfg.dem_file);
    dem_params = struct();
    dem_params.target_resolution = cfg.dem_target_resolution;
    dem_params.fill_nodata = cfg.dem_fill_nodata;
    if ~isempty(cfg.dem_crop_half_size) && cfg.dem_crop_half_size > 0
        h = cfg.dem_crop_half_size;
        dem_params.crop_bounds = [-h, h, -h, h];
    end
    terrain_data = dem_loader(cfg.dem_file, dem_params);
else
    fprintf('Terrain source: synthetic mountain\n');
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
los_verify = los_checker(tm, struct('use_mesh', cfg.use_mesh_los_verify, 'mesh', mesh, 'los_eps', cfg.mesh_los_eps));

%% 2) Radar and binary threat map
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));
center_N = mean(tm.bounds(1:2));
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

nN = floor((tm.bounds(2) - tm.bounds(1)) / threat_params.resolution(1)) + 1;
nE = floor((tm.bounds(4) - tm.bounds(3)) / threat_params.resolution(1)) + 1;
nA = floor((threat_params.alt_range(2) - threat_params.alt_range(1)) / threat_params.resolution(2)) + 1;
fprintf('Estimated threat voxels: %d x %d x %d = %.2fM\n', nN, nE, nA, (nN*nE*nA)/1e6);

threat = threat_map(tm, los_threat, 0.1, threat_params);
radar_ne = select_radar_positions(tm, cfg, cfg.num_radars);
radar_range = max(500, 0.35 * diag_len);
for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + radar_alt_offset;
    radar = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar.R_max = radar_range;
    radar.P_t = 250e3;
    threat.add_radar(radar);
end
fprintf('Radar placement mode: %s | count: %d | range: %.0f m\n', ...
    cfg.radar_position_mode, length(threat.radars), radar_range);
for r = 1:length(threat.radars)
    rp = threat.radars{r}.position;
    fprintf('  Radar-%d @ [N=%.1f, E=%.1f, Alt=%.1f]\n', r, rp(1), rp(2), rp(3));
end

fprintf('Computing deterministic visibility map...\n');
fprintf('Threat grid resolution: horiz=%.1f m, vert=%.1f m\n', threat_params.resolution(1), threat_params.resolution(2));
tic;
threat.compute_map('binary', struct('use_parallel', cfg.use_parallel_threat, 'show_progress', true));
t_threat = toc;
fprintf('Threat map compute time: %.2f s\n', t_threat);

%% 3) Mission and hard-constrained RRT*
visibility_threshold = 0.5;
start_N_fixed = tm.bounds(1) + 0.08 * N_span;
goal_N_fixed  = tm.bounds(2) - 0.08 * N_span;
start_agl = 40;
goal_agl = 40;

[start_NE, start_h, start_risk] = find_hidden_endpoint(tm, threat, start_N_fixed, center_E, start_agl, visibility_threshold, E_span);
[goal_NE, goal_h, goal_risk] = find_hidden_endpoint(tm, threat, goal_N_fixed, center_E, goal_agl, visibility_threshold, E_span);
fprintf('Start risk (binary): %.2f | Goal risk (binary): %.2f\n', start_risk, goal_risk);

start_pos = [start_NE; -(start_h + start_agl)];
goal_pos = [goal_NE; -(goal_h + goal_agl)];

planner_params = struct();
planner_params.max_iter = cfg.dem_crop_half_size * 6;  % Scale max iterations with terrain size
planner_params.base_step_size = max(30, min(75, diag_len / 85));
planner_params.max_step_size = planner_params.base_step_size * 3.5;
planner_params.rewire_radius = planner_params.base_step_size * 2.8;
planner_params.rewire_mode = 'fixed';  % 'fixed' | 'dynamic' | 'knearest'
planner_params.rewire_k_const = 3.2 * exp(1) * (1 + 1/3);
planner_params.rewire_k_min = 28;
planner_params.rewire_k_max = 320;

planner_params.goal_bias = 0.15;
planner_params.goal_bias_after_goal = 0.04;
planner_params.min_clearance = 20;
planner_params.max_flight_alt = 500;
planner_params.alpha = 1.0;
planner_params.beta = 0.0; % no probabilistic radar penalty in binary mode
planner_params.gamma = 0.2;  % preferred AGL penalty weight
planner_params.preferred_agl = 40; % target altitude above terrain [m]
planner_params.shadow_bias = 0.85;
planner_params.use_parallel_rewire = false;
planner_params.animate = cfg.rrt_live_animation;
planner_params.plot_interval = cfg.rrt_live_plot_interval;

% Key settings for deterministic safety:
planner_params.radar_hard_constraint = true;
planner_params.radar_visibility_threshold = visibility_threshold;

% Explicit 3D planning bounds compatible with DEM elevation scale (NED D)
min_alt = terrain_min + planner_params.min_clearance;
max_alt = terrain_max + planner_params.max_flight_alt;
planner_params.bounds = [tm.bounds(1), tm.bounds(2), tm.bounds(3), tm.bounds(4), -max_alt, -min_alt];

if cfg.rrt_live_animation
    fig_rrt_live = figure('Name', 'RRT* Live Planning', 'Position', [70, 70, 1150, 820], 'Color', 'w');
    ax_rrt_live = axes('Parent', fig_rrt_live);
    axes(ax_rrt_live);
    tm.plot();
    hold(ax_rrt_live, 'on');
    [sx_live, sy_live, sz_live] = sphere(24);
    for r = 1:length(threat.radars)
        rr = threat.radars{r};
        plot3(ax_rrt_live, rr.position(1), rr.position(2), rr.position(3), ...
            'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        surf(ax_rrt_live, sx_live * rr.R_max + rr.position(1), ...
            sy_live * rr.R_max + rr.position(2), ...
            sz_live * rr.R_max + rr.position(3), ...
            'FaceColor', [1 0 0], 'FaceAlpha', 0.02, ...
            'EdgeColor', [1 0 0], 'EdgeAlpha', 0.05);
    end
    plot3(ax_rrt_live, start_pos(1), start_pos(2), -start_pos(3), ...
        'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot3(ax_rrt_live, goal_pos(1), goal_pos(2), -goal_pos(3), ...
        'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    xlabel(ax_rrt_live, 'North [m]');
    ylabel(ax_rrt_live, 'East [m]');
    zlabel(ax_rrt_live, 'Altitude [m]');
    title(ax_rrt_live, 'RRT* Live Planning');
    grid(ax_rrt_live, 'on');
    view(ax_rrt_live, 35, 32);
    planner_params.anim_axes = ax_rrt_live;
end

[path_rrt, info] = rrt_star_radar(start_pos, goal_pos, tm, threat, planner_params);
if ~info.success
    error(['RRT* could not find a fully masked path. Try: ', ...
           'increase max_iter, reduce radar range, increase max_flight_alt, ', ...
           'or relax hard constraint for this terrain.']);
end

%% 4) MinSnap smoothing + flatness check
fprintf('\nRunning MinSnap smoothing and flatness verification...\n');

path_simplified = simplify_path(path_rrt, max(25, 0.005 * diag_len));
path_len = 0;
for i = 2:size(path_simplified, 2)
    path_len = path_len + norm(path_simplified(:, i) - path_simplified(:, i-1));
end

smooth_params = struct();
smooth_params.v_max = 14;
smooth_params.a_max = 6;
smooth_params.dt = 0.02;
smooth_params.vel_bc_mode = 'free';
smooth_params.cruise_speed = 9;
smooth_params.max_waypoints = 80;
smooth_params.max_seg_length = 150;
smooth_params.aggressiveness = 2.6;

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

%% 5) LOS verification on final path
n_wp = size(path_rrt, 2);
vis_wp = false(1, n_wp);
in_range_wp = false(1, n_wp);
for i = 1:n_wp
    p_alt = -path_rrt(3, i);
    p = [path_rrt(1, i); path_rrt(2, i); p_alt];
    visible_any = false;
    in_range_any = false;
    for r = 1:length(threat.radars)
        rr = threat.radars{r};
        d = norm(p - rr.position);
        if d <= rr.R_max
            in_range_any = true;
            if los_verify.has_los(rr.position, p)
                visible_any = true;
                break;
            end
        end
    end
    in_range_wp(i) = in_range_any;
    vis_wp(i) = visible_any;
end

n_in_range = sum(in_range_wp);
n_masked_in_range = sum(in_range_wp & ~vis_wp);
n_visible = sum(vis_wp);

fprintf('Path waypoints: %d\n', n_wp);
fprintf('Inside radar sphere: %d\n', n_in_range);
fprintf('Inside sphere but terrain-masked (safe): %d\n', n_masked_in_range);
fprintf('Visible waypoints (should be 0 in hard mode): %d\n', n_visible);

%% 6) Visualization
figure('Name', 'Binary Radar Masking - 3D', 'Position', [80, 80, 1200, 820]);
tm.plot(); hold on;

% Radar markers + coverage spheres
[sx, sy, sz] = sphere(28);
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    plot3(rr.position(1), rr.position(2), rr.position(3), ...
        'r^', 'MarkerSize', 11, 'MarkerFaceColor', 'r');
    surf(sx*rr.R_max + rr.position(1), ...
         sy*rr.R_max + rr.position(2), ...
         sz*rr.R_max + rr.position(3), ...
         'FaceColor', [1 0 0], 'FaceAlpha', 0.03, ...
         'EdgeColor', [1 0 0], 'EdgeAlpha', 0.08);
end

% Path segments: green = masked/safe, red = visible
for i = 1:n_wp-1
    seg_color = [0 0.7 0];
    if vis_wp(i) || vis_wp(i+1)
        seg_color = [0.9 0.1 0.1];
    end
    plot3(path_rrt(1, i:i+1), path_rrt(2, i:i+1), -path_rrt(3, i:i+1), ...
        '-', 'Color', seg_color, 'LineWidth', 3);
end

% Smoothed MinSnap trajectory overlay
if traj_ok
    plot3(traj_poly.pos(1, :), traj_poly.pos(2, :), -traj_poly.pos(3, :), ...
        'c-', 'LineWidth', 2.2);
end

plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

view(35, 32);
grid on;
xlabel('North [m]');
ylabel('East [m]');
zlabel('Altitude [m]');
title(sprintf('Hard Radar Constraint + Terrain Masking (%s)', terrain_label));

% 2D visibility slice at median path altitude
figure('Name', 'Binary Threat Slice', 'Position', [120, 120, 980, 760]);
alt_slice = median(-path_rrt(3, :));
[class_slice, alt_slice_used] = threat.get_binary_horizontal_slice(alt_slice, visibility_threshold, tm);

imagesc(threat.N_vec, threat.E_vec, class_slice);
set(gca, 'YDir', 'normal');
axis equal tight;
colormap([0.35 0.35 0.35; 0.1 0.6 0.1; 0.8 0.1 0.1]);
caxis([1 3]);
colorbar('Ticks', [1, 2, 3], 'TickLabels', {'Below Terrain', 'Hidden/Safe', 'Visible'});
hold on;
plot(path_rrt(1, :), path_rrt(2, :), 'w-', 'LineWidth', 2.5);
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    plot(rr.position(1), rr.position(2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end
plot(start_pos(1), start_pos(2), 'go', 'MarkerSize', 9, 'MarkerFaceColor', 'g');
plot(goal_pos(1), goal_pos(2), 'bs', 'MarkerSize', 9, 'MarkerFaceColor', 'b');

title(sprintf('Binary Visibility Slice at %.1f m Altitude', alt_slice_used));
xlabel('North [m]');
ylabel('East [m]');

if cfg.skymap_gui
    if isfield(cfg, 'skymap_mode') && strcmpi(cfg.skymap_mode, 'viewer')
        skymap_ui = skymap_viewer(cfg, tm, threat, path_rrt, start_pos, goal_pos, visibility_threshold); %#ok<NASGU>
    else
        skymap_ui = skymap_app(cfg, tm, threat, path_rrt, start_pos, goal_pos, visibility_threshold); %#ok<NASGU>
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

function radar_ne = select_radar_positions(tm, cfg, n_radars)
% Select radar positions from terrain data.
% Output: [2 x n_radars] matrix [N; E].

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
        title(sprintf('Click %d radar position(s), then press Enter', n_radars));
        colormap(turbo);
        colorbar;
        hold on;
        if isfield(cfg, 'figure_click_use_contour') && cfg.figure_click_use_contour
            [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
            contour(N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
        end

        fprintf('Figure click mode: select %d radar position(s) on the terrain map.\n', n_radars);
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

        valid = N_flat >= N_min & N_flat <= N_max & E_flat >= E_min & E_flat <= E_max;
        N_flat = N_flat(valid);
        E_flat = E_flat(valid);
        Z_flat = Z_flat(valid);

        [~, order] = sort(Z_flat, 'descend');
        selected = zeros(2, 0);

        for idx = order'
            cand = [N_flat(idx); E_flat(idx)];
            if isempty(selected)
                selected(:, end+1) = cand; %#ok<AGROW>
            else
                d = vecnorm(selected - cand, 2, 1);
                if all(d >= cfg.radar_min_separation)
                    selected(:, end+1) = cand; %#ok<AGROW>
                end
            end
            if size(selected, 2) >= n_radars
                break;
            end
        end

        if size(selected, 2) < n_radars
            for idx = order'
                cand = [N_flat(idx); E_flat(idx)];
                if isempty(selected)
                    selected(:, end+1) = cand; %#ok<AGROW>
                else
                    d = vecnorm(selected - cand, 2, 1);
                    if all(d > 1e-3)
                        selected(:, end+1) = cand; %#ok<AGROW>
                    end
                end
                if size(selected, 2) >= n_radars
                    break;
                end
            end
        end

        if isempty(selected)
            selected = [mean(tm.bounds(1:2)); mean(tm.bounds(3:4))];
        end
        radar_ne = selected;

    otherwise
        error('Unknown cfg.radar_position_mode: %s', cfg.radar_position_mode);
end

if size(radar_ne, 2) < n_radars
    radar_ne = repmat(radar_ne(:, end), 1, n_radars);
elseif size(radar_ne, 2) > n_radars
    radar_ne = radar_ne(:, 1:n_radars);
end

% Clamp to terrain bounds.
radar_ne(1, :) = min(max(radar_ne(1, :), tm.bounds(1)), tm.bounds(2));
radar_ne(2, :) = min(max(radar_ne(2, :), tm.bounds(3)), tm.bounds(4));
end

function [NE, h, risk] = find_hidden_endpoint(tm, threat, N_fixed, E_center, agl, risk_threshold, E_span)
% Find a nearby endpoint with low binary risk by sweeping lateral offsets.

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

    % Prefer hidden points first, then smaller lateral offset.
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
