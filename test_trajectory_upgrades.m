function test_trajectory_upgrades()
%TEST_TRAJECTORY_UPGRADES Verification script for radar path planning upgrades
%
% Verifies:
% 1. Dynamic RCS and Log-Barrier cost in RRT*
% 2. Informed RRT* sampling with Shadow Bias
% 3. Sliding Window Trajectory Generation
%
% Author: Quadrotor Terrain Following Project

clear; clc; close all;

% Add paths
addpath('terrain');
addpath('radar');
addpath('motion_planner');
addpath('safety');

fprintf('=== Verifying Radar Path Planning Upgrades ===\n\n');

%% 1. Setup Environment
fprintf('1. Setting up Environment (Ridge Penetration)...\n');
terrain_params.bounds = [0, 800, -300, 300];
terrain_params.resolution = 10;
terrain_params.amplitude = 80;
terrain_params.wavelength = 300;
terrain_data = terrain_generator('ridge', terrain_params);
tm = terrain_map(terrain_data);
los = los_checker(tm);

% Radar Setup
rcs = 0.01;
threat_params.alt_range = [0, 300];
threat_params.resolution = [20, 20];
threat = threat_map(tm, los, rcs, threat_params);

radar_pos = [600; 0; tm.get_height(600, 0) + 30];
radar1 = radar_site(radar_pos, 'SAM-1', 'tracking');
radar1.R_max = 250;
threat.add_radar(radar1);
threat.compute_map('max');

%% 2. Run RRT* with Upgrades
fprintf('\n2. Running Upgraded RRT*...\n');
start_pos = [50; -200; -(tm.get_height(50, -200) + 40)];
goal_pos = [750; -200; -(tm.get_height(750, -200) + 40)]; % Goal is near radar but offset

params.max_iter = 2000;
params.step_size = 40;
params.goal_bias = 0.15;
params.rewire_radius = 80;
params.alpha = 1.0;
params.beta = 2000; % High weight for radar to test avoidance
params.gamma = 0.5;
params.min_clearance = 15;
params.preferred_alt = 30;

[path_rrt, info] = rrt_star_radar(start_pos, goal_pos, tm, threat, params);

if ~info.success
    warning('RRT* failed to reach goal!');
end

%% 3. Generate Trajectory (Sliding Window)
fprintf('\n3. Generating Trajectory (Sliding Window N=4)...\n');
V_cruise = 10;
window_size = 4;

% Smoothen path geometrically first
path_smooth = smooth_path_geometric(path_rrt, tm, 2);

% Downsample path slightly if needed, but sliding window handles many points well
[traj, debug] = generate_trajectory_sliding_window(path_smooth, V_cruise, window_size);

fprintf('  Generated trajectory: %.2f s duration, %d polynomials\n', ...
    traj.t(end), debug.n_poly);

%% 4. Visualization & Validation
fprintf('\n4. Validating Results...\n');

figure('Name', 'System Verification', 'Position', [100, 100, 1200, 800]);

% 3D View
subplot(2, 2, [1, 3]);
tm.plot_with_path(traj.x);
hold on;
plot3(radar_pos(1), radar_pos(2), radar_pos(3), 'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
[Xs, Ys, Zs] = sphere(20);
surf(radar_pos(1)+radar1.R_max*Xs, radar_pos(2)+radar1.R_max*Ys, radar_pos(3)+radar1.R_max*Zs, ...
    'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('3D Trajectory Verification');
view(30, 30);

% Velocity Profile (Continuity Check)
subplot(2, 2, 2);
vel_norm = vecnorm(traj.v, 2, 1);
plot(traj.t, vel_norm, 'b-', 'LineWidth', 1.5);
hold on;
yline(V_cruise, 'k--');
xlabel('Time [s]');
ylabel('Speed [m/s]');
title('Velocity Profile (Continuity Check)');
grid on;

% Risk Profile
subplot(2, 2, 4);
% Re-evaluate risk using threat map
risks = zeros(1, length(traj.t));
for i = 1:length(traj.t)
    risks(i) = threat.get_risk(traj.x(1,i), traj.x(2,i), -traj.x(3,i));
end
plot(traj.t, risks, 'r-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Detection Prob');
title('Risk Profile (Log-Barrier Validation)');
ylim([0, 1]);
grid on;

%% 5. Continuity Validation
% Check for jumps in velocity or acceleration
max_vel_jump = max(vecnorm(diff(traj.v, 1, 2), 2, 1));
max_acc_jump = max(vecnorm(diff(traj.a, 1, 2), 2, 1));

fprintf('\nContinuity Check:\n');
fprintf('  Max Velocity Delta per step: %.4f m/s\n', max_vel_jump);
fprintf('  Max Acceleration Delta per step: %.4f m/s^2\n', max_acc_jump);

if max_vel_jump < 0.5
    fprintf('  Velocity Continuity: PASS\n');
else
    fprintf('  Velocity Continuity: WARN (Possible discontinuities)\n');
end

fprintf('\nVerification Complete.\n');

end
