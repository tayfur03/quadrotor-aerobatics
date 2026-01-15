%% DEMO_AEROBATICS - Aerobatic Maneuvers Demonstration
%
% Demonstrates various aerobatic maneuvers using the quadrotor simulation.
% Based on maneuvers from Tal et al. (2023) "Aerobatic Trajectory Generation
% for a VTOL Fixed-Wing Aircraft Using Differential Flatness"
%
% Available maneuvers:
%   1. Straight Line    - For testing different yaw modes
%   2. Vertical Loop    - Classic vertical loop in N-D plane (Fig. 3 of paper)
%   3. Barrel Roll      - Helical roll around velocity vector
%   4. Immelmann Turn   - Half loop up (direction reversal, gains altitude)
%   5. Split-S          - Half loop down (direction reversal, loses altitude)
%   6. Climbing Turn    - 270-degree banked turn with altitude gain
%   7. Figure Eight     - Horizontal figure-8 pattern
%
% Author: Quadrotor Simulation Framework
% Based on: Tal & Karaman differential flatness control

clear; clc; close all;

%% Select Maneuver
fprintf('=== Aerobatic Maneuvers Demo ===\n\n');
fprintf('Available maneuvers:\n');
fprintf('  1. Straight Line (yaw mode testing)\n');
fprintf('  2. Vertical Loop\n');
fprintf('  3. Barrel Roll\n');
fprintf('  4. Immelmann Turn\n');
fprintf('  5. Split-S\n');
fprintf('  6. Climbing Turn (270 deg)\n');
fprintf('  7. Figure Eight\n');
fprintf('  8. Racetrack Circuit\n');
fprintf('  9. Run All (Sequential)\n\n');

% Default selection (can be changed)
% selection = 2;  % Run All
selection = input('Select maneuver (1-9): ');

% Default maneuver
%maneuver_type_str = 'split_s';
%maneuver_type_str = 'racetrack';
%maneuver_type_str = 'immelmann';

maneuver_names = {'straight_line', 'vertical_loop', 'barrel_roll', 'immelmann', ...
    'split_s', 'climbing_turn', 'figure_eight', 'racetrack'};

% Determine selection based on maneuver_type_str
% selection = find(strcmpi(maneuver_names, maneuver_type_str), 1);
if isempty(selection)
    warning('Specified maneuver_type_str "%s" not found. Defaulting to Vertical Loop.', maneuver_type_str);
    selection = 2; % Default to Vertical Loop if not found
end

% If you want to run all, uncomment the line below and comment out maneuver_type_str
% selection = 9; % Run All

if selection == 9
    % Run all maneuvers sequentially
    run_all_maneuvers(maneuver_names);
else
    % Run single maneuver
    run_single_maneuver(maneuver_names{selection});
end

%% ==================== MAIN FUNCTIONS ====================

function run_single_maneuver(maneuver_type)
%% Initialize Parameters
params = quad_params_indi();
dt  = 0.002;

%% Generate Maneuver Trajectory
maneuver_params = struct();

% Default parameters (matching paper Fig. 3 for vertical loop)
maneuver_params.radius = 1.0;         % 1m radius as in paper
maneuver_params.V = 3.5;              % ~3.5 m/s nominal speed
maneuver_params.start_pos = [0; 0; -3.5];  % Start at 2.5m altitude
maneuver_params.num_rolls = 1;
maneuver_params.climb_height = 1.0;
maneuver_params.turn_angle = 270;

% Special parameters for straight line
if strcmpi(maneuver_type, 'straight_line')
    maneuver_params.end_pos = [10; 0; -2.5];  % 10m forward
    maneuver_params.V = 3.0;
end

if strcmpi(maneuver_type, 'rectangular')
    maneuver_params.length = 12.0;
    maneuver_params.width = 6.0;
    maneuver_params.radius = 2.0;
    maneuver_params.V = 4.0;
end

fprintf('Generating %s trajectory...\n', maneuver_type);
[waypoints, timePoints, velBC, accBC, info] = ...
    aerobatic_maneuvers(maneuver_type, maneuver_params);

fprintf('Maneuver: %s\n', info.name);
fprintf('Description: %s\n', info.description);
fprintf('Duration: %.2f s\n', info.total_time);
fprintf('Yaw mode: %s\n', info.yaw_mode);
fprintf('Waypoints: %d\n\n', size(waypoints, 2));

%% Plot Waypoints with Velocity Vectors (Before Simulation)
plot_waypoints_with_velocity(waypoints, velBC, info);

%% Generate Smooth Trajectory with minsnappolytraj
Tf = info.total_time;
N = round(Tf / dt);
numSamples = N;

fprintf('Generating minimum-snap trajectory...\n');
[xref, vref, aref, jref, sref] = minsnappolytraj(...
    waypoints, timePoints, numSamples, ...
    'VelocityBoundaryCondition', velBC, ...
    'AccelerationBoundaryCondition', accBC);

%% Initialize Simulation State
[state, motors, filters] = init_sim_state(params, waypoints(:,1));

% Initialize logs
log = init_log(N);

% Yaw parameters
yaw_params = params;
if isfield(info, 'yaw_params')
    fn = fieldnames(info.yaw_params);
    for i = 1:length(fn)
        yaw_params.(fn{i}) = info.yaw_params.(fn{i});
    end
end

%% Main Simulation Loop
fprintf('Running simulation...\n');

for k = 1:N
    t = (k-1) * dt;

    % Get reference
    ref.x = xref(:, k);
    ref.v = vref(:, k);
    ref.a = aref(:, k);
    ref.j = jref(:, k);
    ref.s = sref(:, k);

    % Yaw planning based on recommended mode
    state_ref.v = ref.v;
    state_ref.a = ref.a;
    state_ref.pos = ref.x;

    [ref.psi, ref.psi_dot, ref.psi_ddot] = yaw_planner(t, info.yaw_mode, state_ref, yaw_params);
    yaw_params.last_psi = ref.psi; % Store for next iteration unwrap

    % Update yaw_params with last_psi for continuity
    yaw_params.last_psi = ref.psi;

    % Simulation step
    [state, motors, filters, out] = sim_step(state, motors, filters, ref, params, dt);

    % Log data
    log = log_step(log, k, t, state, ref, out);
end

%% Results
fprintf('\n=== Results ===\n');
final_pos = state.x;
final_error = norm(final_pos - xref(:, end));
fprintf('Final position: [%.3f, %.3f, %.3f] m\n', final_pos);
fprintf('Final error: %.3f m\n', final_error);

position_error = vecnorm(log.x - log.x_ref, 2, 1);
fprintf('Max position error: %.3f m\n', max(position_error));
fprintf('RMS position error: %.3f m\n', sqrt(mean(position_error.^2)));

%% Plots
plot_maneuver_results(log, params, waypoints, velBC, info);

%% Animation
animate_maneuver(log, waypoints, velBC, info, true);
end

function run_all_maneuvers(maneuver_names)
fprintf('Running all maneuvers sequentially...\n\n');

for i = 1:length(maneuver_names)
    fprintf('\n========== Maneuver %d/%d: %s ==========\n\n', ...
        i, length(maneuver_names), upper(maneuver_names{i}));
    run_single_maneuver(maneuver_names{i});
    fprintf('\nPress any key to continue to next maneuver...\n');
    pause;
    close all;
end

fprintf('\nAll maneuvers completed!\n');
end

%% ==================== PLOTTING FUNCTIONS ====================

function plot_waypoints_with_velocity(waypoints, velBC, info)
%PLOT_WAYPOINTS_WITH_VELOCITY Plot waypoints with velocity vectors
%
% Shows the waypoints and their associated velocity boundary conditions
% as arrows, similar to Fig. 3 of the Tal et al. paper.

figure('Name', [info.name ' - Waypoints & Velocity BCs'], ...
    'Position', [50, 50, 800, 700]);

n_wp = size(waypoints, 2);

% Plot waypoints
plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ko', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Waypoints');
hold on;

% Plot velocity vectors at waypoints
vel_scale = 0.3;  % Scale factor for velocity arrows
for i = 1:n_wp
    pos = waypoints(:, i);
    vel = velBC(:, i);

    if ~any(isnan(vel))
        % Draw velocity arrow
        quiver3(pos(1), pos(2), pos(3), ...
            vel(1)*vel_scale, vel(2)*vel_scale, vel(3)*vel_scale, ...
            0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    end
end

% Add legend entry for velocity (only once)
quiver3(nan, nan, nan, nan, nan, nan, 0, 'r', 'LineWidth', 2, ...
    'DisplayName', sprintf('Velocity (x%.1f)', vel_scale));

% Connect waypoints with dashed line
plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'k--', ...
    'LineWidth', 1, 'DisplayName', 'Path (approx)');

% Label waypoints
for i = 1:n_wp
    text(waypoints(1,i), waypoints(2,i), waypoints(3,i) - 0.15, ...
        sprintf('WP%d', i), 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

% Formatting
xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');
title(sprintf('%s - Waypoints with Velocity BCs', info.name));
legend('Location', 'best');
set(gca, 'ZDir', 'reverse');
grid on; axis equal;
view(30, 25);

% Add velocity magnitude annotations
annotation_str = sprintf('Waypoint Velocities:\n');
for i = 1:n_wp
    vel = velBC(:, i);
    if ~any(isnan(vel))
        v_mag = norm(vel);
        annotation_str = sprintf('%sWP%d: [%.1f, %.1f, %.1f] m/s (|v|=%.1f)\n', ...
            annotation_str, i, vel(1), vel(2), vel(3), v_mag);
    else
        annotation_str = sprintf('%sWP%d: [free]\n', annotation_str, i);
    end
end

% Position text box in figure
dim = [0.02 0.02 0.3 0.25];
annotation('textbox', dim, 'String', annotation_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'FontSize', 8, 'FontName', 'FixedWidth');

drawnow;
end

function plot_maneuver_results(log, params, waypoints, velBC, info)
t = log.t;

% 3D Trajectory with velocity vectors
figure('Name', info.name, 'Position', [100, 100, 1200, 800]);

% Main 3D plot
subplot(2,2,[1,3]);

% Plot trajectory colored by speed
speed = vecnorm(log.v, 2, 1);
scatter3(log.x(1,:), log.x(2,:), log.x(3,:), 10, speed, 'filled');
hold on;
colormap(jet);
cb = colorbar;
cb.Label.String = 'Speed [m/s]';

% Plot reference trajectory
plot3(log.x_ref(1,:), log.x_ref(2,:), log.x_ref(3,:), 'k--', ...
    'LineWidth', 1, 'DisplayName', 'Reference');

% For vertical loop: overlay ideal circle for comparison
if isfield(info, 'center') && isfield(info, 'radius')
    theta_ideal = linspace(0, 2*pi, 100);
    center = info.center;
    R = info.radius;
    N_ideal = center(1) + R * sin(theta_ideal);
    E_ideal = center(2) * ones(size(theta_ideal));
    D_ideal = center(3) + R * cos(theta_ideal);
    plot3(N_ideal, E_ideal, D_ideal, 'm-', 'LineWidth', 2, 'DisplayName', 'Ideal Circle');
end

% Plot waypoints
plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ko', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'y', 'DisplayName', 'Waypoints');

% Plot velocity vectors at waypoints
vel_scale = 0.2;
n_wp = size(waypoints, 2);
for i = 1:n_wp
    pos = waypoints(:, i);
    vel = velBC(:, i);
    if ~any(isnan(vel))
        quiver3(pos(1), pos(2), pos(3), ...
            vel(1)*vel_scale, vel(2)*vel_scale, vel(3)*vel_scale, ...
            0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.8);
    end
end

% Start and end markers
plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'gs', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot3(log.x(1,end), log.x(2,end), log.x(3,end), 'r^', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'End');

xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');
title(sprintf('%s - 3D Trajectory', info.name));
legend('Location', 'best');
set(gca, 'ZDir', 'reverse');
grid on; axis equal;
view(30, 25);

% Position error
subplot(2,2,2);
position_error = vecnorm(log.x - log.x_ref, 2, 1);
plot(t, position_error * 100, 'LineWidth', 2);  % Convert to cm
xlabel('Time [s]'); ylabel('Position Error [cm]');
title('Position Tracking Error');
grid on;

% Yaw reference
subplot(2,2,4);
plot(t, rad2deg(log.psi_ref), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Yaw [deg]');
title(sprintf('Yaw Reference (Mode: %s)', info.yaw_mode));
grid on;

% For vertical loop: add 2D N-D plane view
if isfield(info, 'center') && isfield(info, 'radius')
    figure('Name', [info.name ' - N-D Plane (Side View)'], 'Position', [150, 150, 800, 700]);

    % Plot ideal circle
    theta_ideal = linspace(0, 2*pi, 100);
    center = info.center;
    R = info.radius;
    N_ideal = center(1) + R * sin(theta_ideal);
    D_ideal = center(3) + R * cos(theta_ideal);
    plot(N_ideal, D_ideal, 'm-', 'LineWidth', 3, 'DisplayName', 'Ideal Circle');
    hold on;

    % Plot reference trajectory
    plot(log.x_ref(1,:), log.x_ref(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'minsnappolytraj');

    % Plot actual trajectory
    plot(log.x(1,:), log.x(3,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Actual (simulated)');

    % Plot waypoints
    plot(waypoints(1,:), waypoints(3,:), 'ko', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'y', 'DisplayName', 'Waypoints');

    % Plot center of circle
    plot(center(1), center(3), 'rx', 'MarkerSize', 15, 'LineWidth', 3, 'DisplayName', 'Center');

    % Velocity vectors at waypoints
    vel_scale = 0.15;
    for i = 1:size(waypoints, 2)
        pos = waypoints(:, i);
        vel = velBC(:, i);
        if ~any(isnan(vel))
            quiver(pos(1), pos(3), vel(1)*vel_scale, vel(3)*vel_scale, ...
                0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        end
    end

    xlabel('North [m]'); ylabel('Down [m]');
    title(sprintf('%s - Side View (N-D Plane)', info.name));
    legend('Location', 'best');
    set(gca, 'YDir', 'reverse');  % Down is positive
    axis equal; grid on;

    % Add circle deviation statistics
    % Compute distance from reference trajectory to ideal circle
    ref_N = log.x_ref(1,:);
    ref_D = log.x_ref(3,:);
    dist_to_center = sqrt((ref_N - center(1)).^2 + (ref_D - center(3)).^2);
    deviation = abs(dist_to_center - R);

    % Only consider loop section (not approach/exit)
    loop_mask = dist_to_center < R * 1.5;  % Points near the loop
    if any(loop_mask)
        max_dev = max(deviation(loop_mask)) * 100;  % cm
        mean_dev = mean(deviation(loop_mask)) * 100;  % cm
        text_str = sprintf('Circle Deviation:\n  Max: %.1f cm\n  Mean: %.1f cm', max_dev, mean_dev);
        text(0.02, 0.98, text_str, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 10, 'FontName', 'FixedWidth', 'BackgroundColor', 'white');
    end
end

% Additional figure: velocity components
figure('Name', [info.name ' - Velocity'], 'Position', [150, 150, 800, 600]);

subplot(3,1,1);
plot(t, log.v(1,:), 'b', t, log.v_ref(1,:), 'r--', 'LineWidth', 1.5);
ylabel('v_N [m/s]'); legend('Actual', 'Reference');
title('North Velocity');
grid on;

subplot(3,1,2);
plot(t, log.v(2,:), 'b', t, log.v_ref(2,:), 'r--', 'LineWidth', 1.5);
ylabel('v_E [m/s]'); legend('Actual', 'Reference');
title('East Velocity');
grid on;

subplot(3,1,3);
plot(t, log.v(3,:), 'b', t, log.v_ref(3,:), 'r--', 'LineWidth', 1.5);
ylabel('v_D [m/s]'); legend('Actual', 'Reference');
xlabel('Time [s]');
title('Down Velocity');
grid on;

% --- NEW PLOTS ---

% Body Rates (p, q, r)
figure('Name', [info.name ' - Body Rates'], 'Position', [150, 150, 800, 600]);
subplot(3,1,1);
plot(t, log.omega(1,:), 'b', t, log.omega_ref(1,:), 'r--', 'LineWidth', 1.5);
ylabel('p [rad/s]'); legend('p', 'p_{ref}'); title('Roll Rate'); grid on;

subplot(3,1,2);
plot(t, log.omega(2,:), 'b', t, log.omega_ref(2,:), 'r--', 'LineWidth', 1.5);
ylabel('q [rad/s]'); legend('q', 'q_{ref}'); title('Pitch Rate'); grid on;

subplot(3,1,3);
plot(t, log.omega(3,:), 'b', t, log.omega_ref(3,:), 'r--', 'LineWidth', 1.5);
ylabel('r [rad/s]'); legend('r', 'r_{ref}'); xlabel('Time [s]'); title('Yaw Rate'); grid on;

% Attitude Error & Alpha Command
figure('Name', [info.name ' - Attitude Error & Control'], 'Position', [150, 150, 800, 600]);
subplot(2,1,1);
if isfield(log, 'xi_e')
    plot(t, log.xi_e, 'LineWidth', 1.5);
    ylabel('\xi_e [rad]'); legend('x','y','z');
    title('Attitude Error (Vector Part)'); grid on;
end

subplot(2,1,2);
if isfield(log, 'alpha_cmd')
    plot(t, log.alpha_cmd, 'LineWidth', 1.5);
    ylabel('\alpha_{cmd} [rad/s^2]'); xlabel('Time [s]');
    legend('x','y','z');
    title('Angular Acceleration Command'); grid on;
end

% Motor Speeds & Thrust
figure('Name', [info.name ' - Motors & Thrust'], 'Position', [150, 150, 800, 600]);
subplot(2,1,1);
if isfield(log, 'omega_mot')
    plot(t, log.omega_mot, 'LineWidth', 1.5);
    ylabel('\omega_{mot} [rad/s]');
    title('Motor Speeds'); legend('m1','m2','m3','m4'); grid on;

    % Check for saturation (approximate)
    omega_max = 1200;
    max_spd = max(log.omega_mot, [], 'all');
    if max_spd > omega_max * 0.95
        text(t(end)*0.05, max_spd*0.9, sprintf('SATURATION: %.1f', max_spd), 'Color', 'r', 'FontWeight', 'bold');
    end
end

subplot(2,1,2);
if isfield(log, 'T_cmd')
    plot(t, log.T_cmd, 'LineWidth', 1.5);
    ylabel('Thrust [N]'); xlabel('Time [s]');
    title('Total Thrust Command'); grid on;
end
end

function animate_maneuver(log, waypoints, velBC, info, save_video, video_filename)
if nargin < 5
    save_video = false;
end
if nargin < 6
    video_filename = sprintf('%s_animation.mp4', lower(strrep(info.name, ' ', '_')));
end

n_steps = numel(log.t);

fig = figure('Name', [info.name ' Animation'], 'Color', 'w', ...
    'Position', [100, 100, 1000, 800]);
axis_len = 0.3;
lw = 2;

% Draw reference trajectory (faded)
plot3(log.x_ref(1,:), log.x_ref(2,:), log.x_ref(3,:), ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 1);
hold on;

% Draw waypoints with velocity vectors
plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ko', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'y');

vel_scale = 0.2;
n_wp = size(waypoints, 2);
for i = 1:n_wp
    pos = waypoints(:, i);
    vel = velBC(:, i);
    if ~any(isnan(vel))
        quiver3(pos(1), pos(2), pos(3), ...
            vel(1)*vel_scale, vel(2)*vel_scale, vel(3)*vel_scale, ...
            0, 'm', 'LineWidth', 1.5, 'MaxHeadSize', 0.8);
    end
end

grid on; axis equal;
set(gca, 'ZDir', 'reverse');
xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');

% Drone body axes
h_body_x = plot3([0,0], [0,0], [0,0], 'r-', 'LineWidth', lw);
h_body_y = plot3([0,0], [0,0], [0,0], 'g-', 'LineWidth', lw);
h_body_z = plot3([0,0], [0,0], [0,0], 'b-', 'LineWidth', lw);

% Actual trail
h_trail = plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'c-', 'LineWidth', 2);

% Velocity vector on drone
h_vel = quiver3(0, 0, 0, 0, 0, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);

h_title = title(sprintf('%s - Time: 0.00 s', info.name));
view(45, 25);

% Setup video writer
if save_video
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 95;
    open(v);
end

steps_skip = 10; % Daha akici ve yavas animasyon icin 50'den 10'a dusuruldu
trail_x = []; trail_y = []; trail_z = [];

for k = 1:steps_skip:n_steps
    pos = log.x(:,k);
    vel = log.v(:,k);
    q_k = log.q(:,k);
    R = quat_to_R(q_k);

    % Update trail
    trail_x = [trail_x, pos(1)];
    trail_y = [trail_y, pos(2)];
    trail_z = [trail_z, pos(3)];
    set(h_trail, 'XData', trail_x, 'YData', trail_y, 'ZData', trail_z);

    % Body axes
    tip_x = pos + R(:,1) * axis_len;
    tip_y = pos + R(:,2) * axis_len;
    tip_z = pos + R(:,3) * axis_len;

    set(h_body_x, 'XData', [pos(1), tip_x(1)], 'YData', [pos(2), tip_x(2)], 'ZData', [pos(3), tip_x(3)]);
    set(h_body_y, 'XData', [pos(1), tip_y(1)], 'YData', [pos(2), tip_y(2)], 'ZData', [pos(3), tip_y(3)]);
    set(h_body_z, 'XData', [pos(1), tip_z(1)], 'YData', [pos(2), tip_z(2)], 'ZData', [pos(3), tip_z(3)]);

    % Velocity vector
    v_scale = 0.15;
    set(h_vel, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), ...
        'UData', vel(1)*v_scale, 'VData', vel(2)*v_scale, 'WData', vel(3)*v_scale);

    % Compute current yaw from rotation matrix
    yaw_deg = rad2deg(atan2(R(2,1), R(1,1)));
    speed = norm(vel);
    set(h_title, 'String', sprintf('%s - Time: %.2f s | Yaw: %.1f deg | Speed: %.1f m/s', ...
        info.name, log.t(k), yaw_deg, speed));

    drawnow;

    if save_video
        frame = getframe(fig);
        writeVideo(v, frame);
    else
        pause(0.05); % Yavaslatmak icin bekleme suresi artirildi
    end
end

if save_video
    close(v);
    fprintf('Video saved: %s\n', video_filename);
end
end

%% ==================== SIMULATION FUNCTIONS ====================

function [state, motors, filters] = init_sim_state(params, start_pos)
state.x = start_pos;
state.v = [0;0;0];
state.q = [1;0;0;0];
state.omega = [0;0;0];
motors.omega = params.omega_hover * ones(4,1);
filters.a_f = zeros(3,1);
filters.tau_bz_f = zeros(3,1);
end

function log = init_log(N)
log.t               = zeros(1,N);
log.x               = zeros(3,N);
log.x_ref           = zeros(3,N);
log.v               = zeros(3,N);
log.v_ref           = zeros(3,N);
log.a_true          = zeros(3,N);
log.a_ref           = zeros(3,N);
log.q               = zeros(4,N);
log.omega           = zeros(3,N);
log.omega_ref       = zeros(3,N);
log.psi_ref         = zeros(1,N);
log.T_cmd           = zeros(1,N);
log.omega_mot       = zeros(4,N);
% Additional log fields for detailed analysis
log.xi_e            = zeros(3,N);
log.alpha_cmd       = zeros(3,N);
end

function [state, motors, filters, out] = sim_step(state, motors, filters, ref, params, dt)
g = params.g;
m = params.m;
J = params.J;
e3 = [0;0;1];

R  = quat_to_R(state.q);
bz = R(:,3);

u_vec = params.G1 * (motors.omega.^2);
mu    = u_vec(1:3);
T     = u_vec(4);

a_true = g*e3 - (T/m)*bz;

x_dot     = state.v;
v_dot     = a_true;
q_dot     = 0.5 * quat_mul(state.q, [0; state.omega]);
omega_dot = J \ (mu - cross(state.omega, J*state.omega));

state.x     = state.x + dt*x_dot;
state.v     = state.v + dt*v_dot;
state.q     = quat_normalize(state.q + dt*q_dot);
state.omega = state.omega + dt*omega_dot;

a_b = R' * (a_true - g*e3);
a_meas = R*a_b + g*e3;
filters.a_f = filters.a_f + (dt/params.tau_f) * (a_meas - filters.a_f);

tau_bz = a_true - g*e3;
filters.tau_bz_f = filters.tau_bz_f + (dt/params.tau_f) * (tau_bz - filters.tau_bz_f);

a_c = params.Kx*(ref.x - state.x) + params.Kv*(ref.v - state.v) + ...
    params.Ka*(ref.a - filters.a_f) + ref.a;

tau_bz_c = filters.tau_bz_f + a_c - filters.a_f;

tau_norm = norm(tau_bz_c);
if tau_norm < 1e-6
    T_cmd = m * 9.81;
else
    T_cmd = m * tau_norm;
end

R_curr = quat_to_R(state.q);
q_inc = compute_incremental_attitude_cmd(R_curr, tau_bz_c, ref.psi);

q_cmd_fake  = q_inc;
q_curr_fake = [1;0;0;0];

q_ideal_next = quat_mul(state.q, q_inc);
R_ideal_next = quat_to_R(q_ideal_next);

[omega_ref, omega_dot_ref] = flatness(R_ideal_next, tau_norm, ref.j, ref.s, ref.psi_dot, ref.psi_ddot);

[alpha_cmd, xi_e] = attitude_pd(q_cmd_fake, q_curr_fake, state.omega, omega_ref, omega_dot_ref, params);
mu_cmd = J*alpha_cmd + cross(state.omega, J*state.omega);

omega_mot_cmd = motor_inversion(mu_cmd, T_cmd, params);
motors.omega = motors.omega + (dt/params.tau_m) * (omega_mot_cmd - motors.omega);

out.a_true        = a_true;
out.a_ref         = ref.a;
out.omega_ref     = omega_ref;
out.T_cmd         = T_cmd;
out.omega_mot     = motors.omega;
out.psi_ref       = ref.psi;
out.xi_e          = xi_e;
out.alpha_cmd     = alpha_cmd;
end

function log = log_step(log, k, t, state, ref, out)
log.t(k)          = t;
log.x(:,k)        = state.x;
log.x_ref(:,k)    = ref.x;
log.v(:,k)        = state.v;
log.v_ref(:,k)    = ref.v;
log.a_true(:,k)   = out.a_true;
log.a_ref(:,k)    = out.a_ref;
log.q(:,k)        = state.q;
log.omega(:,k)    = state.omega;
log.omega_ref(:,k) = out.omega_ref;
log.psi_ref(k)    = out.psi_ref;
log.T_cmd(k)      = out.T_cmd;
log.omega_mot(:,k) = out.omega_mot;
log.xi_e(:,k)      = out.xi_e;
log.alpha_cmd(:,k) = out.alpha_cmd;
end

function q_inc = compute_incremental_attitude_cmd(R_curr, tau_bz_c, psi_ref)
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
dot_prod   = dot(current_z, bz_des_body);

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


