% Test script to diagnose Mission 2 trajectory issues
% Run this to see exactly what's happening with RRT waypoints

clear; clc; close all;
addpath('motion_planner');

%% Replicate Mission 2 Setup
start_pos = [0; 0; -2];
goal_pos = [12; 6; -2];

% Define obstacles
obstacles = {};
obstacles{1} = struct('type', 'sphere', 'center', [5; 2; -2], 'radius', 1.5, 'active', true);
obstacles{2} = struct('type', 'sphere', 'center', [4; 4; -2.5], 'radius', 1.2, 'active', true);
obstacles{3} = struct('type', 'sphere', 'center', [10; 5; -1.5], 'radius', 0.8, 'active', true);

bounds = struct('min', [-2; -2; -6], 'max', [15; 10; 0]);

%% Setup Waypoint Manager
wm = waypoint_manager();
wm.rrt_params.max_iter = 8000;
wm.rrt_params.step_size = 0.4;
wm.rrt_params.goal_bias = 0.15;
wm.rrt_params.safety_margin = 0.5;

wm.set_mission(start_pos, goal_pos, obstacles, bounds);

% Plan path
fprintf('Planning RRT path...\n');
[success, waypoints] = wm.plan_path();

if ~success
    error('RRT failed');
end

fprintf('RRT found %d waypoints:\n', size(waypoints, 2));
disp(waypoints');

%% Calculate mission time
path_length = 0;
for i = 1:size(waypoints, 2)-1
    path_length = path_length + norm(waypoints(:,i+1) - waypoints(:,i));
end
mission_speed = 2.0;
Tf = max(8.0, path_length / mission_speed + 4.0);
fprintf('Path length: %.2f m, Duration: %.1f s\n\n', path_length, Tf);

%% Generate trajectory
fprintf('Generating trajectory...\n');
traj = wm.generate_trajectory(Tf);

%% Analyze trajectory
fprintf('\n=== Trajectory Analysis ===\n');

% Check for issues
has_nan = any(isnan(traj.pos(:))) || any(isnan(traj.vel(:))) || any(isnan(traj.acc(:)));
has_inf = any(isinf(traj.pos(:))) || any(isinf(traj.vel(:))) || any(isinf(traj.acc(:)));
fprintf('Has NaN: %d, Has Inf: %d\n', has_nan, has_inf);

% Check ranges
fprintf('\nPosition range:\n');
fprintf('  N: [%.2f, %.2f]\n', min(traj.pos(1,:)), max(traj.pos(1,:)));
fprintf('  E: [%.2f, %.2f]\n', min(traj.pos(2,:)), max(traj.pos(2,:)));
fprintf('  D: [%.2f, %.2f]\n', min(traj.pos(3,:)), max(traj.pos(3,:)));

fprintf('\nVelocity range:\n');
fprintf('  vN: [%.2f, %.2f]\n', min(traj.vel(1,:)), max(traj.vel(1,:)));
fprintf('  vE: [%.2f, %.2f]\n', min(traj.vel(2,:)), max(traj.vel(2,:)));
fprintf('  vD: [%.2f, %.2f]\n', min(traj.vel(3,:)), max(traj.vel(3,:)));

fprintf('\nAcceleration range:\n');
fprintf('  aN: [%.2f, %.2f]\n', min(traj.acc(1,:)), max(traj.acc(1,:)));
fprintf('  aE: [%.2f, %.2f]\n', min(traj.acc(2,:)), max(traj.acc(2,:)));
fprintf('  aD: [%.2f, %.2f]\n', min(traj.acc(3,:)), max(traj.acc(3,:)));

fprintf('\nJerk range:\n');
fprintf('  jN: [%.2f, %.2f]\n', min(traj.jerk(1,:)), max(traj.jerk(1,:)));
fprintf('  jE: [%.2f, %.2f]\n', min(traj.jerk(2,:)), max(traj.jerk(2,:)));
fprintf('  jD: [%.2f, %.2f]\n', min(traj.jerk(3,:)), max(traj.jerk(3,:)));

fprintf('\nSnap range:\n');
fprintf('  sN: [%.2f, %.2f]\n', min(traj.snap(1,:)), max(traj.snap(1,:)));
fprintf('  sE: [%.2f, %.2f]\n', min(traj.snap(2,:)), max(traj.snap(2,:)));
fprintf('  sD: [%.2f, %.2f]\n', min(traj.snap(3,:)), max(traj.snap(3,:)));

%% Check segment times
fprintf('\n=== Segment Analysis ===\n');
fprintf('Segment times: ');
fprintf('%.2f ', traj.segment_times);
fprintf('\n');

n_segments = length(traj.segment_times);
t_segments = [0, cumsum(traj.segment_times)];
fprintf('Cumulative times: ');
fprintf('%.2f ', t_segments);
fprintf('\n');

%% Check speed at segment boundaries (waypoints)
fprintf('\n=== Speed at Waypoints ===\n');
for i = 1:n_segments+1
    if i == 1
        t_wp = 0;
    elseif i == n_segments+1
        t_wp = traj.t(end);
    else
        t_wp = t_segments(i);
    end

    % Find closest index
    [~, idx] = min(abs(traj.t - t_wp));
    speed = norm(traj.vel(:, idx));
    fprintf('Waypoint %d (t=%.2f): pos=[%.2f, %.2f, %.2f], speed=%.3f m/s\n', ...
        i, t_wp, traj.pos(:,idx), speed);
end

%% Plot comprehensive analysis
figure('Name', 'Mission 2 Trajectory Analysis', 'Position', [50, 50, 1400, 900]);

% 3D trajectory
subplot(2,3,1);
plot3(traj.pos(1,:), traj.pos(2,:), traj.pos(3,:), 'b-', 'LineWidth', 1.5);
hold on;
plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% Draw obstacles
for i = 1:length(obstacles)
    obs = obstacles{i};
    [X, Y, Z] = sphere(20);
    X = obs.radius * X + obs.center(1);
    Y = obs.radius * Y + obs.center(2);
    Z = obs.radius * Z + obs.center(3);
    surf(X, Y, Z, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
xlabel('N'); ylabel('E'); zlabel('D');
title('3D Trajectory');
set(gca, 'ZDir', 'reverse');
grid on; axis equal;
view(30, 25);

% Position vs time
subplot(2,3,2);
plot(traj.t, traj.pos(1,:), 'r-', traj.t, traj.pos(2,:), 'g-', traj.t, traj.pos(3,:), 'b-');
hold on;
% Mark waypoint times
for i = 2:n_segments
    xline(t_segments(i), 'k--', 'Alpha', 0.5);
end
xlabel('Time [s]'); ylabel('Position [m]');
title('Position vs Time');
legend('N', 'E', 'D');
grid on;

% Velocity vs time
subplot(2,3,3);
plot(traj.t, traj.vel(1,:), 'r-', traj.t, traj.vel(2,:), 'g-', traj.t, traj.vel(3,:), 'b-');
hold on;
for i = 2:n_segments
    xline(t_segments(i), 'k--', 'Alpha', 0.5);
end
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Velocity vs Time');
legend('vN', 'vE', 'vD');
grid on;

% Speed vs time
subplot(2,3,4);
speed = vecnorm(traj.vel, 2, 1);
plot(traj.t, speed, 'k-', 'LineWidth', 1.5);
hold on;
for i = 2:n_segments
    xline(t_segments(i), 'r--', 'Alpha', 0.7);
end
xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Speed vs Time (red lines = waypoint times)');
grid on;
ylim([0, max(speed)*1.2]);

% Acceleration vs time
subplot(2,3,5);
plot(traj.t, traj.acc(1,:), 'r-', traj.t, traj.acc(2,:), 'g-', traj.t, traj.acc(3,:), 'b-');
hold on;
for i = 2:n_segments
    xline(t_segments(i), 'k--', 'Alpha', 0.5);
end
xlabel('Time [s]'); ylabel('Acceleration [m/s^2]');
title('Acceleration vs Time');
legend('aN', 'aE', 'aD');
grid on;

% Jerk vs time
subplot(2,3,6);
plot(traj.t, traj.jerk(1,:), 'r-', traj.t, traj.jerk(2,:), 'g-', traj.t, traj.jerk(3,:), 'b-');
hold on;
for i = 2:n_segments
    xline(t_segments(i), 'k--', 'Alpha', 0.5);
end
xlabel('Time [s]'); ylabel('Jerk [m/s^3]');
title('Jerk vs Time');
legend('jN', 'jE', 'jD');
grid on;

%% Check if large values are the issue
fprintf('\n=== Checking for Large Values ===\n');
max_acc = max(abs(traj.acc(:)));
max_jerk = max(abs(traj.jerk(:)));
max_snap = max(abs(traj.snap(:)));
fprintf('Max |acceleration|: %.2f m/s^2\n', max_acc);
fprintf('Max |jerk|: %.2f m/s^3\n', max_jerk);
fprintf('Max |snap|: %.2f m/s^4\n', max_snap);

if max_acc > 50
    fprintf('WARNING: Acceleration seems too high!\n');
end
if max_jerk > 200
    fprintf('WARNING: Jerk seems too high!\n');
end
if max_snap > 1000
    fprintf('WARNING: Snap seems too high!\n');
end

fprintf('\n=== Test Complete ===\n');
