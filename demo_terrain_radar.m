%% DEMO_TERRAIN_RADAR - Terrain-Following with Radar Avoidance
%
% Demonstrates the complete pipeline for terrain-aware, radar-minimizing
% trajectory planning and control for a quadrotor UAV.
%
% Pipeline:
%   1. Generate synthetic terrain
%   2. Create radar threat environment
%   3. Plan path using RRT* with radar cost
%   4. Smooth trajectory with minsnappolytraj
%   5. Simulate flight with CBF safety filter
%
% Author: Quadrotor Terrain Following Project
% Based on: Tal & Karaman differential flatness control

clear; clc; close all;

%% Add paths
addpath('terrain');
addpath('radar');
addpath('safety');
addpath('motion_planner');

%% ==================== CONFIGURATION ====================

% Scenario selection
scenario = 'ridge_penetration';  % 'ridge_penetration', 'valley_flight', 'multi_radar'

fprintf('=== Terrain-Following Radar Avoidance Demo ===\n');
fprintf('Scenario: %s\n\n', scenario);

%% ==================== TERRAIN SETUP ====================

fprintf('Step 1: Generating terrain...\n');

switch scenario
    case 'ridge_penetration'
        terrain_params.bounds = [0, 800, -300, 300];
        terrain_params.resolution = 10;
        terrain_params.amplitude = 100;
        terrain_params.wavelength = 300;
        terrain_params.type = 'ridge';
        terrain_data = terrain_generator(terrain_params);

    case 'valley_flight'
        terrain_params.bounds = [0, 1000, -400, 400];
        terrain_params.resolution = 10;
        terrain_params.amplitude = 100;
        terrain_data = terrain_generator('valley', terrain_params);

    case 'multi_radar'
        terrain_params.bounds = [0, 1000, -500, 500];
        terrain_params.resolution = 10;
        terrain_params.amplitude = 60;
        terrain_data = terrain_generator('hills', terrain_params);

    otherwise
        error('Unknown scenario: %s', scenario);
end

tm = terrain_map(terrain_data);
los = los_checker(tm);

fprintf('  Terrain: %s, bounds: [%.0f, %.0f] x [%.0f, %.0f]\n', ...
    terrain_data.type, terrain_params.bounds);

%% ==================== RADAR SETUP ====================

fprintf('\nStep 2: Setting up radar threat environment...\n');

% UAV radar cross section (small quadrotor)
rcs = 0.01;  % 0.01 m^2 typical for small UAV

% Create threat map
threat_params.alt_range = [0, 300];
threat_params.resolution = [20, 20];
threat = threat_map(tm, los, rcs, threat_params);

% Add radar(s) based on scenario
switch scenario
    case 'ridge_penetration'
        % Single radar behind the ridge
        radar_pos = [600; 0; tm.get_height(600, 0) + 30];  % On terrain
        radar1 = radar_site(radar_pos, 'SAM-1', 'tracking');
        radar1.R_max = 200;  % Reduced for demo scale
        threat.add_radar(radar1);

    case 'valley_flight'
        % Radars on valley walls
        radar_pos1 = [500; 200; tm.get_height(500, 200) + 20];
        radar_pos2 = [500; -200; tm.get_height(500, -200) + 20];
        radar1 = radar_site(radar_pos1, 'AA-East', 'tracking');
        radar2 = radar_site(radar_pos2, 'AA-West', 'tracking');
        radar1.R_max = 400;
        radar2.R_max = 400;
        threat.add_radar(radar1);
        threat.add_radar(radar2);

    case 'multi_radar'
        % Multiple overlapping radars
        positions = [300, 100; 500, -150; 700, 50];
        for i = 1:size(positions, 1)
            pos = [positions(i,1); positions(i,2); ...
                tm.get_height(positions(i,1), positions(i,2)) + 25];
            radar = radar_site(pos, sprintf('Radar-%d', i), 'tracking');
            radar.R_max = 350;
            threat.add_radar(radar);
        end
end

fprintf('  Added %d radar site(s)\n', length(threat.radars));
for i = 1:length(threat.radars)
    fprintf('    %s\n', threat.radars{i}.get_info());
end

% Compute threat map
fprintf('\nStep 3: Computing threat map...\n');
threat.compute_map('max');

%% ==================== PATH PLANNING ====================

fprintf('\nStep 4: Planning radar-aware path...\n');

% Define mission
switch scenario
    case 'ridge_penetration'
        start_pos = [50; -250; -(tm.get_height(50, -250) + 40)];   % 40m AGL, E=-250 offset (outside 200m radar range)
        goal_pos = [750; 250; -(tm.get_height(750, -250) + 40)];   % E=-250 offset

    case 'valley_flight'
        start_pos = [50; 0; -(tm.get_height(50, 0) + 30)];
        goal_pos = [950; 0; -(tm.get_height(950, 0) + 30)];

    case 'multi_radar'
        start_pos = [50; -200; -(tm.get_height(50, -200) + 50)];
        goal_pos = [950; 200; -(tm.get_height(950, 200) + 50)];
end

% RRT* parameters - TUNED for reliable path finding
planner_params.max_iter = 5000;     % Increased for better convergence
planner_params.step_size = 40;      % Slightly smaller for better resolution
planner_params.goal_bias = 0.20;    % Increased from 0.15 for faster goal connection
planner_params.rewire_radius = 80;  % Proportional to step size
planner_params.alpha = 1.0;         % Distance weight
planner_params.beta = 30;           % REDUCED from 100 - was causing path failure
planner_params.gamma = 0.3;         % Altitude cost weight (moderate)
planner_params.min_clearance = 15;  % 15m terrain clearance
planner_params.preferred_alt = 30;  % Prefer low-moderate altitude for terrain masking
planner_params.shadow_bias = 0.4;   % Reduced sampling bias for visible points
planner_params.max_flight_alt = 400; % Maximum AGL altitude for sampling

% Plan path
[path_rrt, plan_info] = rrt_star_radar(start_pos, goal_pos, tm, threat, planner_params);

fprintf('  Path found: %d waypoints, %.1f m length, %.2f integrated risk\n', ...
    size(path_rrt, 2), plan_info.path_length, plan_info.path_risk);

% Visualize planned path
figure('Name', 'Path Planning Result', 'Position', [100, 100, 1000, 800]);

% 3D view
subplot(2, 2, [1, 3]);
tm.plot_with_path(path_rrt);
hold on;

% Plot radar positions
for i = 1:length(threat.radars)
    radar = threat.radars{i};
    plot3(radar.position(1), radar.position(2), radar.position(3), ...
        'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
end

title('RRT* Path with Radar Avoidance');
view(30, 30);

% Risk along path
subplot(2, 2, 2);
[risks, cumulative] = threat.get_risk_along_path(path_rrt);
plot(risks, 'r-', 'LineWidth', 2);
xlabel('Waypoint Index');
ylabel('Detection Probability');
title('Risk Along Path');
ylim([0, 1]);
grid on;

% Altitude profile
subplot(2, 2, 4);
path_dist = [0, cumsum(vecnorm(diff(path_rrt, 1, 2), 2, 1))];
path_alt = -path_rrt(3, :);
terrain_under_path = zeros(1, size(path_rrt, 2));
for i = 1:size(path_rrt, 2)
    terrain_under_path(i) = tm.get_height(path_rrt(1,i), path_rrt(2,i));
end

plot(path_dist, path_alt, 'b-', 'LineWidth', 2);
hold on;
plot(path_dist, terrain_under_path, 'k-', 'LineWidth', 2);
fill([path_dist, fliplr(path_dist)], ...
    [terrain_under_path, zeros(size(terrain_under_path))], ...
    [0.5 0.4 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xlabel('Distance Along Path [m]');
ylabel('Altitude [m]');
title('Altitude Profile');
legend('Flight Path', 'Terrain', 'Location', 'best');
grid on;

drawnow;

%% ==================== TRAJECTORY GENERATION ====================

fprintf('\nStep 5: Generating smooth trajectory...\n');

% Smoothen path geometrically first
waypoints_smooth = smooth_path_geometric(path_rrt, tm, 2);

% Global Minimum Snap Trajectory Generation
V_cruise = 8;

% Filter out duplicate waypoints to prevent zero-duration segments
wps_unique = waypoints_smooth(:, 1);
for i = 2:size(waypoints_smooth, 2)
    if norm(waypoints_smooth(:, i) - wps_unique(:, end)) > 1e-3
        wps_unique = [wps_unique, waypoints_smooth(:, i)];
    end
end
waypoints_smooth = wps_unique;

dists = sqrt(sum(diff(waypoints_smooth, 1, 2).^2, 1));
% Ensure minimum time duration per segment
dt_segments = max(dists / V_cruise, 0.1);
times = [0, cumsum(dt_segments)];
total_duration = times(end);

% Boundary conditions (Rest-to-Rest)
num_wps = size(waypoints_smooth, 2);
vel_bc = nan(3, num_wps);
vel_bc(:, 1) = [0;0;0];
vel_bc(:, end) = [0;0;0];

acc_bc = nan(3, num_wps);
acc_bc(:, 1) = [0;0;0];
acc_bc(:, end) = [0;0;0];

% Generate PP
numSamples = num_wps * 10;
% Robustly capture outputs (Output count varies by MATLAB version/toolbox)
[~, ~, ~, ~, ~, pp] = minsnappolytraj(waypoints_smooth, times, numSamples, ...
    'MinSegmentTime', 1e-3, ...
    'VelocityBoundaryCondition', vel_bc, ...
    'AccelerationBoundaryCondition', acc_bc);



fprintf('  Trajectory: %.1f s duration (Global Minimum Snap)\n', total_duration);

%% ==================== SIMULATION PRE-PROCESSING ====================

% Resample trajectory for high-frequency control (500 Hz)
% This is CRITICAL for INDI filter stability (tau_f ~ 0.004s requires dt < 0.002s)
dt_sim = 0.002;
t_sim = 0:dt_sim:total_duration;
N_sim = length(t_sim);

fprintf('\nStep 5b: Resampling for control (%.0f Hz)...\n', 1/dt_sim);

% pp is already defined above

% Evaluate derivatives using spline utilities
xref = ppval(pp, t_sim);
vref = ppval(fnder(pp, 1), t_sim);
aref = ppval(fnder(pp, 2), t_sim);
jref = ppval(fnder(pp, 3), t_sim);
sref = ppval(fnder(pp, 4), t_sim);

N = N_sim;
dt = dt_sim;

%% ==================== SIMULATION WITH CBF ====================

fprintf('\nStep 6: Simulating with CBF safety filter...\n');

% Initialize simulation
params = quad_params_indi();

state.x = xref(:, 1);
state.v = [0; 0; 0];
state.q = [1; 0; 0; 0];
state.omega = [0; 0; 0];

motors.omega = params.omega_hover * ones(4, 1);
filters.a_f = zeros(3, 1);
filters.tau_bz_f = zeros(3, 1);

% CBF safety filter (disabled for now - can be enabled later)
USE_CBF = true;

if USE_CBF
    cbf_params.h_min = 15;
    cbf_params.P_max = 0.5;
    cbf_params.v_max = 15;
    cbf_params.alpha_terrain = 1.5;
    cbf_params.enable_radar = false;  % Disable radar CBF for stability
    cbf = cbf_safety_filter(tm, threat, cbf_params);
end

% Initialize logs
log.t = zeros(1, N);
log.x = zeros(3, N);
log.x_ref = zeros(3, N);
log.v = zeros(3, N);
log.v_ref = zeros(3, N);
log.q = zeros(4, N);
log.alt = zeros(1, N);
log.terrain_h = zeros(1, N);
log.risk = zeros(1, N);
log.clearance = zeros(1, N);

% Yaw planner parameters
yaw_params = params;
yaw_params.last_psi = 0;

% Simulation loop
fprintf('  Running simulation...\n');

for k = 1:N
    t = (k-1) * dt;

    % Get reference
    ref.x = xref(:, k);
    ref.v = vref(:, k);
    ref.a = aref(:, k);
    ref.j = jref(:, k);
    ref.s = sref(:, k);

    % Yaw planning (tangent to velocity)
    state_ref.v = ref.v;
    state_ref.a = ref.a;
    state_ref.j = ref.j; % Pass jerk for accurate derivatives
    state_ref.pos = ref.x;
    [ref.psi, ref.psi_dot, ref.psi_ddot] = yaw_planner(t, 'tangent', state_ref, yaw_params);
    yaw_params.last_psi = ref.psi;

    % Standard simulation step (same as demo_aerobatics)
    [state, motors, filters, out] = sim_step_terrain(state, motors, filters, ref, params, dt);

    % Log data
    log.t(k) = t;
    log.x(:, k) = state.x;
    log.x_ref(:, k) = ref.x;
    log.v(:, k) = state.v;
    log.v_ref(:, k) = ref.v;
    log.q(:, k) = state.q;
    log.alt(k) = -state.x(3);
    log.terrain_h(k) = tm.get_height(state.x(1), state.x(2));
    log.clearance(k) = log.alt(k) - log.terrain_h(k);

    if threat.computed
        log.risk(k) = threat.get_risk(state.x(1), state.x(2), -state.x(3));
    end

    % Progress
    if mod(k, round(N/10)) == 0
        fprintf('    Progress: %d%%\n', round(100*k/N));
    end
end

fprintf('  Simulation complete.\n');

%% ==================== RESULTS ====================

fprintf('\n=== RESULTS ===\n');

% Final position error
final_error = norm(state.x - waypoints_smooth(:, end));
fprintf('Final position error: %.2f m\n', final_error);

% Tracking error
pos_error = vecnorm(log.x - log.x_ref, 2, 1);
fprintf('Max position error: %.2f m\n', max(pos_error));
fprintf('RMS position error: %.2f m\n', sqrt(mean(pos_error.^2)));

% Clearance and risk summaries
min_clearance = min(log.clearance);
mean_clearance = mean(log.clearance);
max_risk = max(log.risk);
mean_risk = mean(log.risk);

fprintf('Min clearance: %.2f m\n', min_clearance);
fprintf('Mean clearance: %.2f m\n', mean_clearance);
fprintf('Max detection probability: %.3f\n', max_risk);
fprintf('Mean detection probability: %.3f\n', mean_risk);

% Collect summary for report tables
summary.scenario = scenario;
summary.num_waypoints = size(path_rrt, 2);
summary.path_length_m = plan_info.path_length;
summary.path_risk = plan_info.path_risk;
summary.traj_duration_s = total_duration;
summary.final_error_m = final_error;
summary.max_error_m = max(pos_error);
summary.rms_error_m = sqrt(mean(pos_error.^2));
summary.min_clearance_m = min_clearance;
summary.mean_clearance_m = mean_clearance;
summary.max_risk = max_risk;
summary.mean_risk = mean_risk;

fprintf('\n=== SUMMARY (copy to report table) ===\n');
fprintf(['Scenario: %s\n', ...
    'Waypoints: %d\n', ...
    'Path length: %.1f m\n', ...
    'Integrated risk: %.2f\n', ...
    'Trajectory duration: %.1f s\n', ...
    'Final error: %.2f m\n', ...
    'Max error: %.2f m\n', ...
    'RMS error: %.2f m\n', ...
    'Min clearance: %.2f m\n', ...
    'Mean clearance: %.2f m\n', ...
    'Max detection prob.: %.3f\n', ...
    'Mean detection prob.: %.3f\n'], ...
    summary.scenario, summary.num_waypoints, summary.path_length_m, ...
    summary.path_risk, summary.traj_duration_s, summary.final_error_m, ...
    summary.max_error_m, summary.rms_error_m, summary.min_clearance_m, ...
    summary.mean_clearance_m, summary.max_risk, summary.mean_risk);

% Annotate plots with key metrics (kept short for readability)
figure(findobj('Name', 'Path Planning Result'));
subplot(2, 2, 2);
text(0.02, 0.90, sprintf('Path risk: %.2f', plan_info.path_risk), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'w');

subplot(2, 2, 4);
text(0.02, 0.90, sprintf('Min clearance: %.1f m', min_clearance), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'w');

% Safety metrics
fprintf('Min terrain clearance: %.1f m\n', min(log.clearance));
fprintf('Max detection probability: %.2f\n', max(log.risk));
fprintf('Integrated risk exposure: %.2f\n', sum(log.risk) * dt);

%% ==================== VISUALIZATION ====================

% 3D trajectory with threat overlay
figure('Name', 'Mission Results', 'Position', [100, 100, 1200, 800]);

subplot(2, 2, [1, 3]);
tm.plot_with_path(log.x);
hold on;

% Flight path is already shown by the terrain plot

% Radar positions and range spheres
for i = 1:length(threat.radars)
    radar = threat.radars{i};
    plot3(radar.position(1), radar.position(2), radar.position(3), ...
        'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

    % Add radar range sphere
    [Xs, Ys, Zs] = sphere(30);
    Xs = radar.position(1) + radar.R_max * Xs;
    Ys = radar.position(2) + radar.R_max * Ys;
    Zs = radar.position(3) + radar.R_max * Zs;
    surf(Xs, Ys, Zs, 'FaceColor', 'red', 'FaceAlpha', 0.15, ...
        'EdgeColor', 'none');
end

title('Flight Path with Radar Range Visualization');
view(30, 30);

% Tracking performance
subplot(2, 2, 2);
plot(log.t, pos_error * 100, 'b-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Position Error [cm]');
title('Tracking Error');
grid on;

% Risk exposure over time (always show, regardless of USE_CBF)
subplot(2, 2, 4);
plot(log.t, log.risk, 'r-', 'LineWidth', 1.5);
hold on;
plot(log.t, 0.5 * ones(size(log.t)), 'k--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Detection Probability');
title('Radar Exposure Over Time');
legend('P_{det}', 'Threshold (0.5)');
ylim([0, 1]);
grid on;

drawnow;

%% ==================== DEDICATED 3D MAP ====================
% Full-screen 3D visualization with terrain, path, and radar
figure('Name', '3D Mission Visualization', 'Position', [50, 50, 1400, 900]);

% Plot terrain surface (colored by elevation)
tm.plot(gcf);
hold on;

% Plot flight path (colored by altitude)
N_path = log.x(1, :);
E_path = log.x(2, :);
alt_path = -log.x(3, :);  % Convert NED to altitude
plot3(N_path, E_path, alt_path, 'g-', 'LineWidth', 3);

% Plot RRT* waypoints (planning nodes)
plot3(path_rrt(1, :), path_rrt(2, :), -path_rrt(3, :), ...
    'ko--', 'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', 'k');

% Start and end markers
plot3(N_path(1), E_path(1), alt_path(1), 'go', 'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot3(N_path(end), E_path(end), alt_path(end), 'rs', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'LineWidth', 2);

% Goal marker (target position) - larger and with label
goal_alt = -waypoints_smooth(3, end); % Assuming waypoints(:, end) is the goal position
plot3(waypoints_smooth(1, end), waypoints_smooth(2, end), goal_alt, 'mp', 'MarkerSize', 25, 'MarkerFaceColor', 'm', 'LineWidth', 3);
text(waypoints_smooth(1, end), waypoints_smooth(2, end), goal_alt + 20, sprintf('GOAL \n[%.0f, %.0f]', waypoints_smooth(1, end), waypoints_smooth(2, end)), ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'magenta', 'HorizontalAlignment', 'center');

% Start marker with label
text(waypoints_smooth(1, 1), waypoints_smooth(2, 1), -waypoints_smooth(3, 1) + 20, 'START', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'green', 'HorizontalAlignment', 'center');

% Radar positions and range spheres
for i = 1:length(threat.radars)
    radar = threat.radars{i};

    % Radar marker
    plot3(radar.position(1), radar.position(2), radar.position(3), ...
        'r^', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'LineWidth', 2);

    % Radar range sphere (semi-transparent)
    [Xs, Ys, Zs] = sphere(40);
    Xs = radar.position(1) + radar.R_max * Xs;
    Ys = radar.position(2) + radar.R_max * Ys;
    Zs = radar.position(3) + radar.R_max * Zs;
    % Only show upper hemisphere (above terrain)
    Zs(Zs < 0) = NaN;
    surf(Xs, Ys, Zs, 'FaceColor', 'red', 'FaceAlpha', 0.2, ...
        'EdgeColor', 'none');

    % Add label
    text(radar.position(1), radar.position(2), radar.position(3) + 50, ...
        sprintf('%s \nR_{max}=%dm', radar.name, round(radar.R_max)), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red', ...
        'HorizontalAlignment', 'center');
end

% Formatting
xlabel('North [m]', 'FontSize', 12);
ylabel('East [m]', 'FontSize', 12);
zlabel('Altitude [m]', 'FontSize', 12);
title('3D Mission Map: Terrain + Flight Path + Radar Coverage', 'FontSize', 14);
axis equal;
grid on;
view(45, 30);
rotate3d on;  % Enable mouse rotation

% Add legend
legend('Terrain', 'Flight Path', 'RRT* Waypoints', 'Start', 'End', 'Radar', 'Radar Range', ...
    'Location', 'northeast');

fprintf('\n=== 3D Visualization Ready ===\n');
fprintf('Use mouse to rotate the 3D view!\n');

drawnow;

% 2D map view with trajectory
figure('Name', '2D Mission Map', 'Position', [150, 150, 900, 700]);
tm.plot_2d();
hold on;

% Plot trajectory
plot(log.x(1, :), log.x(2, :), 'r-', 'LineWidth', 2);
plot(log.x(1, 1), log.x(2, 1), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g');
plot(log.x(1, end), log.x(2, end), 'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

% Plot radars
for i = 1:length(threat.radars)
    radar = threat.radars{i};
    plot(radar.position(1), radar.position(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    % Draw range circle
    theta = linspace(0, 2*pi, 100);
    plot(radar.position(1) + radar.R_max * cos(theta), ...
        radar.position(2) + radar.R_max * sin(theta), 'r--', 'LineWidth', 2);
end

title('Mission Map (Top View)');
legend('Terrain Contours', 'Flight Path', 'Start', 'End', 'Radar', 'Radar Range', 'Location', 'best');

fprintf('\nDemo complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function [state, motors, filters, out] = sim_step_terrain(state, motors, filters, ref, params, dt)
% Simulation step for terrain following demo
% Same as demo_aerobatics sim_step
g = params.g;
m = params.m;
J = params.J;
e3 = [0; 0; 1];

R = quat_to_R(state.q);
bz = R(:, 3);

u_vec = params.G1 * (motors.omega.^2);
mu = u_vec(1:3);
T = u_vec(4);

a_true = g*e3 - (T/m)*bz;

% Integrate state (dynamics first, then control)
x_dot = state.v;
v_dot = a_true;
q_dot = 0.5 * quat_mul(state.q, [0; state.omega]);
omega_dot = J \ (mu - cross(state.omega, J * state.omega));

state.x = state.x + dt * x_dot;
state.v = state.v + dt * v_dot;
state.q = quat_normalize(state.q + dt * q_dot);
state.omega = state.omega + dt * omega_dot;

% Filters
a_b = R' * (a_true - g*e3);
a_meas = R * a_b + g*e3;
filters.a_f = filters.a_f + (dt/params.tau_f) * (a_meas - filters.a_f);
tau_bz = a_true - g*e3;
filters.tau_bz_f = filters.tau_bz_f + (dt/params.tau_f) * (tau_bz - filters.tau_bz_f);

% Position control
a_c = params.Kx*(ref.x - state.x) + params.Kv*(ref.v - state.v) + ...
    params.Ka*(ref.a - filters.a_f) + ref.a;
tau_bz_c = filters.tau_bz_f + a_c - filters.a_f;

tau_norm = norm(tau_bz_c);
if tau_norm < 1e-6
    T_cmd = m * 9.81;
else
    T_cmd = m * tau_norm;
end

% Attitude control
R_curr = quat_to_R(state.q);
q_inc = compute_incremental_attitude_cmd_terrain(R_curr, tau_bz_c, ref.psi);

q_cmd_fake = q_inc;
q_curr_fake = [1; 0; 0; 0];

q_ideal_next = quat_mul(state.q, q_inc);
R_ideal_next = quat_to_R(q_ideal_next);

[omega_ref, omega_dot_ref] = flatness(R_ideal_next, tau_norm, ref.j, ref.s, ref.psi_dot, ref.psi_ddot);

[alpha_cmd, ~] = attitude_pd(q_cmd_fake, q_curr_fake, state.omega, omega_ref, omega_dot_ref, params);
mu_cmd = J * alpha_cmd + cross(state.omega, J * state.omega);

omega_mot_cmd = motor_inversion(mu_cmd, T_cmd, params);
motors.omega = motors.omega + (dt/params.tau_m) * (omega_mot_cmd - motors.omega);

out.a_true = a_true;
out.T_cmd = T_cmd;
end

function q_inc = compute_incremental_attitude_cmd_terrain(R_curr, tau_bz_c, psi_ref)
t_norm = norm(tau_bz_c);
if t_norm < 1e-6
    t_des_in = [0;0;-1];
else
    t_des_in = tau_bz_c / t_norm;
end

bz_des_in = -t_des_in;
bz_des_body = R_curr' * bz_des_in;

current_z = [0;0;1];
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
psi_curr = atan2(R_new_in(2,1), R_new_in(1,1));

psi_err = psi_ref - psi_curr;
while psi_err > pi,  psi_err = psi_err - 2*pi; end
while psi_err < -pi, psi_err = psi_err + 2*pi; end

q_yaw = [cos(psi_err/2); 0; 0; sin(psi_err/2)];
q_inc = quat_mul(q_tilt, q_yaw);

if q_inc(1) < 0
    q_inc = -q_inc;
end
end
