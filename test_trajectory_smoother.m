% Test script for trajectory_smoother
% Run this to diagnose issues with the minimum-snap trajectory generation

addpath('motion_planner');

%% Simple 3-waypoint test
fprintf('=== Testing Trajectory Smoother ===\n\n');

% Simple waypoints: start -> middle -> end
waypoints = [0, 5, 10;    % N
             0, 2, 0;     % E
             0, -1, -2];  % D

total_time = 10.0;

fprintf('Waypoints:\n');
disp(waypoints);

% Generate trajectory
try
    traj = trajectory_smoother(waypoints, total_time);
    fprintf('Trajectory generated successfully!\n');
    fprintf('  Time samples: %d\n', length(traj.t));
    fprintf('  Duration: %.2f s\n', traj.t(end));
catch ME
    fprintf('ERROR generating trajectory:\n');
    fprintf('  %s\n', ME.message);
    fprintf('  Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    return;
end

%% Check trajectory values
fprintf('\n=== Checking Trajectory Values ===\n');

% Check for NaN or Inf
has_nan_pos = any(isnan(traj.pos(:)));
has_nan_vel = any(isnan(traj.vel(:)));
has_nan_acc = any(isnan(traj.acc(:)));
has_inf_pos = any(isinf(traj.pos(:)));
has_inf_vel = any(isinf(traj.vel(:)));
has_inf_acc = any(isinf(traj.acc(:)));

fprintf('NaN check - pos: %d, vel: %d, acc: %d\n', has_nan_pos, has_nan_vel, has_nan_acc);
fprintf('Inf check - pos: %d, vel: %d, acc: %d\n', has_inf_pos, has_inf_vel, has_inf_acc);

% Check ranges
fprintf('\nPosition range:\n');
fprintf('  N: [%.2f, %.2f]\n', min(traj.pos(1,:)), max(traj.pos(1,:)));
fprintf('  E: [%.2f, %.2f]\n', min(traj.pos(2,:)), max(traj.pos(2,:)));
fprintf('  D: [%.2f, %.2f]\n', min(traj.pos(3,:)), max(traj.pos(3,:)));

fprintf('\nVelocity range:\n');
fprintf('  N: [%.2f, %.2f]\n', min(traj.vel(1,:)), max(traj.vel(1,:)));
fprintf('  E: [%.2f, %.2f]\n', min(traj.vel(2,:)), max(traj.vel(2,:)));
fprintf('  D: [%.2f, %.2f]\n', min(traj.vel(3,:)), max(traj.vel(3,:)));

fprintf('\nAcceleration range:\n');
fprintf('  N: [%.2f, %.2f]\n', min(traj.acc(1,:)), max(traj.acc(1,:)));
fprintf('  E: [%.2f, %.2f]\n', min(traj.acc(2,:)), max(traj.acc(2,:)));
fprintf('  D: [%.2f, %.2f]\n', min(traj.acc(3,:)), max(traj.acc(3,:)));

%% Check boundary conditions
fprintf('\n=== Checking Boundary Conditions ===\n');

% Start position
fprintf('Start position: expected [%.1f, %.1f, %.1f], got [%.2f, %.2f, %.2f]\n', ...
    waypoints(:,1), traj.pos(:,1));

% End position
fprintf('End position: expected [%.1f, %.1f, %.1f], got [%.2f, %.2f, %.2f]\n', ...
    waypoints(:,end), traj.pos(:,end));

% Start velocity (should be zero)
fprintf('Start velocity: expected [0, 0, 0], got [%.2f, %.2f, %.2f]\n', traj.vel(:,1));

% End velocity (should be zero)
fprintf('End velocity: expected [0, 0, 0], got [%.2f, %.2f, %.2f]\n', traj.vel(:,end));

%% Check middle waypoint (should pass through but NOT stop)
fprintf('\n=== Checking Middle Waypoint ===\n');

% Find time when we should be at middle waypoint (approximately t=5)
mid_time = total_time / 2;
mid_idx = find(traj.t >= mid_time, 1, 'first');

fprintf('At t=%.1f (mid-trajectory):\n', traj.t(mid_idx));
fprintf('  Position: [%.2f, %.2f, %.2f] (expected near [5, 2, -1])\n', traj.pos(:,mid_idx));
fprintf('  Velocity: [%.2f, %.2f, %.2f] (should NOT be zero)\n', traj.vel(:,mid_idx));
fprintf('  Speed: %.2f m/s\n', norm(traj.vel(:,mid_idx)));

%% Plot results
figure('Name', 'Trajectory Smoother Test', 'Position', [100, 100, 1200, 800]);

% 3D trajectory
subplot(2,3,1);
plot3(traj.pos(1,:), traj.pos(2,:), traj.pos(3,:), 'b-', 'LineWidth', 1.5);
hold on;
plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('N'); ylabel('E'); zlabel('D');
title('3D Trajectory');
grid on;
set(gca, 'ZDir', 'reverse');
legend('Trajectory', 'Waypoints');

% Position vs time
subplot(2,3,2);
plot(traj.t, traj.pos(1,:), 'r-', traj.t, traj.pos(2,:), 'g-', traj.t, traj.pos(3,:), 'b-');
xlabel('Time [s]'); ylabel('Position [m]');
title('Position vs Time');
legend('N', 'E', 'D');
grid on;

% Velocity vs time
subplot(2,3,3);
plot(traj.t, traj.vel(1,:), 'r-', traj.t, traj.vel(2,:), 'g-', traj.t, traj.vel(3,:), 'b-');
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Velocity vs Time');
legend('vN', 'vE', 'vD');
grid on;

% Speed vs time
subplot(2,3,4);
speed = vecnorm(traj.vel, 2, 1);
plot(traj.t, speed, 'k-', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Speed vs Time (should NOT drop to zero at waypoints)');
grid on;

% Acceleration vs time
subplot(2,3,5);
plot(traj.t, traj.acc(1,:), 'r-', traj.t, traj.acc(2,:), 'g-', traj.t, traj.acc(3,:), 'b-');
xlabel('Time [s]'); ylabel('Acceleration [m/s^2]');
title('Acceleration vs Time');
legend('aN', 'aE', 'aD');
grid on;

% Jerk vs time
subplot(2,3,6);
plot(traj.t, traj.jerk(1,:), 'r-', traj.t, traj.jerk(2,:), 'g-', traj.t, traj.jerk(3,:), 'b-');
xlabel('Time [s]'); ylabel('Jerk [m/s^3]');
title('Jerk vs Time');
legend('jN', 'jE', 'jD');
grid on;

fprintf('\n=== Test Complete ===\n');
fprintf('Check the figure for visual verification.\n');
fprintf('Key: Speed plot should NOT drop to zero at t=5s (middle waypoint)\n');
