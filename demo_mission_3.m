%% DEMO_MISSION_3 - Pop-up Threat Avoidance
%
% Task 3: Navigate to goal with a pop-up threat appearing mid-flight.
% This demonstrates dynamic replanning:
%   1. Start mission with no obstacles (clear path)
%   2. Pop-up threat appears at specified time
%   3. Detect collision risk with remaining path
%   4. Replan trajectory in real-time to avoid threat
%   5. Continue to goal
%
% Author: Motion Planner Module
% Based on: demo_indi_6dof_ned.m

clear; clc; close all;

% Add motion_planner folder to path
addpath('motion_planner');

%% Initialize Parameters
params = quad_params_indi();
g   = params.g;
m   = params.m;
J   = params.J;
e3  = [0;0;1];

dt = 0.002;       % [s] Simulation timestep

%% Define Mission with Pop-up Threat
% Start position
start_pos = [0; 0; -2];  % Starting at 2m altitude

% Goal position (direct path will pass through threat location)
goal_pos = [15; 0; -2];

% Pop-up threat configuration
% Threat placed further ahead to give drone time to react and replan
popup_time = 4.0;  % Threat appears at t=3.0 seconds (earlier detection)
popup_pos = [10; 0; -2];  % Appears further along path (was [8; 0; -2])
popup_radius = 2.0;  % Large enough to require significant deviation

% Workspace bounds for RRT replanning
bounds = struct(...
    'min', [-2; -8; -6], ...
    'max', [20; 8; 0]);

% Mission duration (will extend if needed for replanning)
Tf_initial = 12.0;  % Initial mission time

fprintf('=== Mission 3: Pop-up Threat Avoidance ===\n');
fprintf('Start:       [%.1f, %.1f, %.1f] m\n', start_pos);
fprintf('Goal:        [%.1f, %.1f, %.1f] m\n', goal_pos);
fprintf('Threat pos:  [%.1f, %.1f, %.1f] m\n', popup_pos);
fprintf('Threat radius: %.1f m\n', popup_radius);
fprintf('Threat appears at: t=%.1f s\n\n', popup_time);

%% Setup Obstacle Manager for Pop-up Threats
obs_mgr = obstacle_manager();

% Add pop-up threat (inactive until activation time)
threat_params = struct('radius', popup_radius);
obs_mgr.add_popup_threat('sphere', popup_pos, threat_params, popup_time);

fprintf('Obstacle manager initialized with 1 pop-up threat\n');

%% Setup Initial Waypoint Manager (No Obstacles Initially)
wm = waypoint_manager();

% Configure RRT parameters for quick replanning (use standard RRT for speed)
wm.rrt_params.algorithm = 'rrt';  % 'rrt' for fast replanning (vs 'rrt_star' for quality)
wm.rrt_params.max_iter = 5000;
wm.rrt_params.step_size = 0.5;
wm.rrt_params.goal_bias = 0.2;
wm.rrt_params.safety_margin = 0.8;  % Increased margin for better avoidance

% Check RRT algorithm type and set wiring radius if using RRT*
if strcmp(wm.rrt_params.algorithm, 'rrt_star')
    wm.rrt_params.wiring_radius = 1.5;  % Set wiring radius for RRT*
    fprintf('Wiring radius set to %.1f for RRT* algorithm\n', wm.rrt_params.wiring_radius);
end

% Plan initial path (no obstacles)
wm.set_mission(start_pos, goal_pos, {}, bounds);
[success, waypoints] = wm.plan_path();
if ~success
    error('Initial path planning failed!');
end

% Generate initial trajectory
traj = wm.generate_trajectory(Tf_initial);
fprintf('Initial trajectory planned (direct path to goal)\n\n');

%% Initialize Simulation State
Tf = Tf_initial + 5.0;  % Extra time for potential replanning
N = round(Tf/dt);
[state, motors, filters] = init_sim_state(params, start_pos);

% Initialize logs
log = init_log(N);

% Replanning state
has_replanned = false;
replan_time = -1;
pre_replan_waypoints = waypoints;

% Saving video
save_video = false;

%% Main Simulation Loop with Dynamic Replanning
fprintf('Running simulation with real-time replanning...\n');

actual_steps = 0;
for k = 1:N
    t = (k-1)*dt;
    actual_steps = k;

    % Update obstacle manager (check for pop-up activation)
    [new_threat, threat_indices] = obs_mgr.update(t);

    % If new threat appeared, check if replanning is needed
    if new_threat && ~has_replanned
        fprintf('\n[%.2fs] NEW THREAT DETECTED!\n', t);

        % Get currently active obstacles
        active_obs = obs_mgr.get_active_obstacles();

        % Check if remaining path collides with new threat
        [needs_replan, reason] = wm.check_replan_needed(state.x, active_obs);

        if needs_replan
            fprintf('[%.2fs] %s - Initiating replanning...\n', t, reason);

            % Replan from current position with current velocity for smooth transition
            % Pass current velocity and filtered acceleration for C1/C2 continuity
            [replan_success, new_traj] = wm.replan(state.x, t, active_obs, state.v, filters.a_f);

            if replan_success
                has_replanned = true;
                replan_time = t;
                fprintf('[%.2fs] Replanning complete - new path generated\n', t);
            else
                fprintf('[%.2fs] WARNING: Replanning failed! Continuing on original path.\n', t);
            end
        else
            fprintf('[%.2fs] Current path is still collision-free\n', t);
        end
    end

    % Get reference from waypoint manager
    try
        ref = wm.get_reference(t);
    catch
        % If trajectory ended, hold at final position
        ref.pos = goal_pos;
        ref.vel = [0;0;0];
        ref.acc = [0;0;0];
        ref.jerk = [0;0;0];
        ref.snap = [0;0;0];
        ref.psi = 0;
        ref.psi_dot = 0;
        ref.psi_ddot = 0;
    end

    % Create state_ref struct for yaw_planner (expects .v and .a fields)
    state_ref.v = ref.vel;
    state_ref.a = ref.acc;

    % Add yaw planning
    [ref.psi, ref.psi_dot, ref.psi_ddot] = yaw_planner(t, 'tangent', state_ref, params);

    % Run simulation step
    [state, motors, filters, out] = sim_step(state, motors, filters, ref, params, dt);

    % Real-time collision monitoring
    active_obs = obs_mgr.get_active_obstacles();
    [is_collision, min_dist, ~] = collision_checker(state.x, active_obs, 0.1);
    if is_collision
        fprintf('[%.2fs] COLLISION WARNING! Min distance: %.3f m\n', t, min_dist);
    end

    % Log data
    log = log_step(log, k, t, state, ref, out);

    % Check if goal reached
    if wm.goal_reached(state.x, 0.3) && t > popup_time + 2.0
        fprintf('[%.2fs] Goal reached! Ending simulation.\n', t);
        break;
    end
end

% Trim logs to actual steps
log = trim_log(log, actual_steps);

%% Results
fprintf('\n=== Results ===\n');
final_pos = state.x;
final_error = norm(final_pos - goal_pos);
fprintf('Final position: [%.3f, %.3f, %.3f] m\n', final_pos);
fprintf('Final error:    %.3f m\n', final_error);

if has_replanned
    fprintf('Replanning occurred at t=%.2f s\n', replan_time);
else
    fprintf('No replanning was needed\n');
end

% Check minimum clearance during entire flight
active_obs = obs_mgr.get_active_obstacles();
min_clearance = inf;
min_clearance_time = 0;
for k = 1:size(log.x, 2)
    if log.t(k) >= popup_time  % Only check after threat appeared
        [~, dist, ~] = collision_checker(log.x(:,k), active_obs, 0);
        if dist < min_clearance
            min_clearance = dist;
            min_clearance_time = log.t(k);
        end
    end
end
fprintf('Minimum threat clearance: %.3f m (at t=%.2f s)\n', min_clearance, min_clearance_time);

if final_error < 0.5 && min_clearance > 0
    fprintf('Mission SUCCESS: Goal reached, threat avoided!\n');
elseif min_clearance <= 0
    fprintf('Mission FAILED: Collision with threat!\n');
else
    fprintf('Mission INCOMPLETE: Did not reach goal within tolerance.\n');
end

%% Plots
active_obs = obs_mgr.get_active_obstacles();
plot_mission_results(log, params, goal_pos, active_obs, wm.waypoints, ...
    pre_replan_waypoints, popup_time, replan_time, 'Mission 3: Pop-up Threat Avoidance');

%% Animation
animate_mission_popup(log, goal_pos, popup_pos, popup_radius, popup_time, ...
    pre_replan_waypoints, wm.waypoints, replan_time, save_video);

%% ==================== LOCAL FUNCTIONS ====================

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
log.a_c             = zeros(3,N);
log.a_f             = zeros(3,N);
log.a_ref           = zeros(3,N);
log.err_eq19        = zeros(3,N);
log.omega           = zeros(3,N);
log.omega_ref       = zeros(3,N);
log.omega_dot_ref   = zeros(3,N);
log.q               = zeros(4,N);
log.q_cmd           = zeros(4,N);
log.xi_e            = zeros(3,N);
log.theta_e         = zeros(1,N);
log.alpha_cmd       = zeros(3,N);
log.tau_norm        = zeros(1,N);
log.T_cmd           = zeros(1,N);
log.omega_mot       = zeros(4,N);
log.omega_mot_cmd   = zeros(4,N);
end

function log = trim_log(log, N)
fields = fieldnames(log);
for i = 1:length(fields)
    f = fields{i};
    if size(log.(f), 2) > N
        log.(f) = log.(f)(:, 1:N);
    end
end
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

a_c = params.Kx*(ref.pos - state.x) + params.Kv*(ref.vel - state.v) + ...
    params.Ka*(ref.acc - filters.a_f) + ref.acc;

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

[omega_ref, omega_dot_ref] = flatness(R_ideal_next, tau_norm, ref.jerk, ref.snap, ref.psi_dot, ref.psi_ddot);

[alpha_cmd, xi_e] = attitude_pd(q_cmd_fake, q_curr_fake, state.omega, omega_ref, omega_dot_ref, params);
mu_cmd = J*alpha_cmd + cross(state.omega, J*state.omega);

omega_mot_cmd = motor_inversion(mu_cmd, T_cmd, params);
motors.omega = motors.omega + (dt/params.tau_m) * (omega_mot_cmd - motors.omega);

err_eq19 = a_true - (tau_bz - filters.tau_bz_f + filters.a_f);

out.a_true        = a_true;
out.a_c           = a_c;
out.a_f           = filters.a_f;
out.a_ref         = ref.acc;
out.err_eq19      = err_eq19;
out.omega_ref     = omega_ref;
out.omega_dot_ref = omega_dot_ref;
out.alpha_cmd     = alpha_cmd;
out.xi_e          = xi_e;
out.q_cmd         = q_inc;
out.tau_norm      = tau_norm;
out.T_cmd         = T_cmd;
out.omega_mot     = motors.omega;
out.omega_mot_cmd = omega_mot_cmd;
end

function log = log_step(log, k, t, state, ref, out)
log.t(k)          = t;
log.x(:,k)        = state.x;
log.x_ref(:,k)    = ref.pos;
log.v(:,k)        = state.v;
log.v_ref(:,k)    = ref.vel;
log.a_true(:,k)   = out.a_true;
log.a_c(:,k)      = out.a_c;
log.a_f(:,k)      = out.a_f;
log.a_ref(:,k)    = out.a_ref;
log.err_eq19(:,k) = out.err_eq19;
log.q(:,k)        = state.q;
log.q_cmd(:,k)    = out.q_cmd;
log.omega(:,k)    = state.omega;
log.omega_ref(:,k) = out.omega_ref;
log.omega_dot_ref(:,k) = out.omega_dot_ref;
log.alpha_cmd(:,k) = out.alpha_cmd;
log.xi_e(:,k)      = out.xi_e;
log.tau_norm(k)    = out.tau_norm;
log.T_cmd(k)       = out.T_cmd;
log.omega_mot(:,k) = out.omega_mot;
log.omega_mot_cmd(:,k) = out.omega_mot_cmd;
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

function plot_mission_results(log, params, goal, obstacles, new_waypoints, ...
    old_waypoints, popup_time, replan_time, title_str)
t = log.t;

% Find replan index
if replan_time > 0
    replan_idx = find(t >= replan_time, 1);
else
    replan_idx = length(t);
end

% 3D Trajectory with threat
figure('Name', title_str, 'Position', [100, 100, 1000, 700]);

% Path before replanning
plot3(log.x(1,1:replan_idx), log.x(2,1:replan_idx), log.x(3,1:replan_idx), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'Path (Before Threat)');
hold on;

% Path after replanning
if replan_idx < length(t)
    plot3(log.x(1,replan_idx:end), log.x(2,replan_idx:end), log.x(3,replan_idx:end), ...
        'c-', 'LineWidth', 2, 'DisplayName', 'Path (After Replanning)');
end

% Original planned path (dashed)
plot3(old_waypoints(1,:), old_waypoints(2,:), old_waypoints(3,:), ...
    'r--', 'LineWidth', 1.5, 'DisplayName', 'Original Plan');

% New planned path
if ~isequal(old_waypoints, new_waypoints)
    plot3(new_waypoints(1,:), new_waypoints(2,:), new_waypoints(3,:), ...
        'g-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Replanned Path');
end

% Start and goal
plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'gs', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot3(goal(1), goal(2), goal(3), 'r^', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Goal');

% Mark replan point
if replan_time > 0 && replan_idx <= length(t)
    plot3(log.x(1,replan_idx), log.x(2,replan_idx), log.x(3,replan_idx), ...
        'm*', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'Replan Point');
end

% Draw pop-up threat
for i = 1:length(obstacles)
    obs = obstacles{i};
    if strcmp(obs.type, 'sphere')
        [X, Y, Z] = sphere(30);
        X = obs.radius * X + obs.center(1);
        Y = obs.radius * Y + obs.center(2);
        Z = obs.radius * Z + obs.center(3);
        surf(X, Y, Z, 'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0.6, ...
            'EdgeColor', 'none');
    end
end

xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');
title(sprintf('%s\nThreat appeared at t=%.1fs, Replanned at t=%.1fs', ...
    title_str, popup_time, replan_time));
legend('Location', 'best');
set(gca, 'ZDir', 'reverse');
grid on; axis equal;
view(30, 25);

% Top-down view
figure('Name', [title_str ' - Top View']);

% Path segments
plot(log.x(1,1:replan_idx), log.x(2,1:replan_idx), 'b-', 'LineWidth', 2);
hold on;
if replan_idx < length(t)
    plot(log.x(1,replan_idx:end), log.x(2,replan_idx:end), 'c-', 'LineWidth', 2);
end

% Original and new paths
plot(old_waypoints(1,:), old_waypoints(2,:), 'r--', 'LineWidth', 1.5);
if ~isequal(old_waypoints, new_waypoints)
    plot(new_waypoints(1,:), new_waypoints(2,:), 'g-o', 'LineWidth', 1.5, 'MarkerSize', 6);
end

% Markers
plot(log.x(1,1), log.x(2,1), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot(goal(1), goal(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

if replan_time > 0 && replan_idx <= length(t)
    plot(log.x(1,replan_idx), log.x(2,replan_idx), 'm*', 'MarkerSize', 15, 'LineWidth', 2);
end

% Draw threat circle
theta_circle = linspace(0, 2*pi, 50);
for i = 1:length(obstacles)
    obs = obstacles{i};
    if strcmp(obs.type, 'sphere')
        x_circle = obs.center(1) + obs.radius * cos(theta_circle);
        y_circle = obs.center(2) + obs.radius * sin(theta_circle);
        fill(x_circle, y_circle, [1 0.3 0.3], 'FaceAlpha', 0.5, ...
            'EdgeColor', 'r', 'LineWidth', 2);
    end
end

xlabel('North [m]'); ylabel('East [m]');
title([title_str ' - Top View']);
legend('Before Threat', 'After Replan', 'Original Plan', 'New Plan', ...
    'Start', 'Goal', 'Replan Point', 'Location', 'best');
grid on; axis equal;

% Timeline plot
figure('Name', [title_str ' - Timeline']);

% Distance to goal
dist_to_goal = zeros(1, length(t));
for k = 1:length(t)
    dist_to_goal(k) = norm(log.x(:,k) - goal);
end

subplot(2,1,1);
plot(t, dist_to_goal, 'b-', 'LineWidth', 2);
hold on;
xline(popup_time, 'r--', 'LineWidth', 2, 'Label', 'Threat Appears');
if replan_time > 0
    xline(replan_time, 'm--', 'LineWidth', 2, 'Label', 'Replan');
end
xlabel('Time [s]'); ylabel('Distance to Goal [m]');
title('Distance to Goal Over Time');
grid on;

% Distance to threat (after it appears)
dist_to_threat = inf(1, length(t));
for k = 1:length(t)
    if t(k) >= popup_time
        [~, dist, ~] = collision_checker(log.x(:,k), obstacles, 0);
        dist_to_threat(k) = dist;
    end
end

subplot(2,1,2);
valid_idx = t >= popup_time;
plot(t(valid_idx), dist_to_threat(valid_idx), 'r-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 1);
if replan_time > 0
    xline(replan_time, 'm--', 'LineWidth', 2, 'Label', 'Replan');
end
xlabel('Time [s]'); ylabel('Distance to Threat [m]');
title('Distance to Threat Over Time (After Appearance)');
grid on;
ylim([-1, max(dist_to_threat(valid_idx)) + 1]);
end

function animate_mission_popup(log, goal, popup_pos, popup_radius, popup_time, ...
    old_waypoints, new_waypoints, replan_time, save_video, video_filename)
% ANIMATE_MISSION_POPUP Animate drone flight with pop-up threat and optional MP4 saving
%
% Inputs:
%   save_video    - (optional) true to save as MP4 (default: true)
%   video_filename - (optional) Output filename (default: 'mission_3_animation.mp4')

if nargin < 9
    save_video = true;
end
if nargin < 10
    video_filename = 'mission_3_animation.mp4';
end

n_steps = numel(log.t);

fig = figure('Name', 'Mission 3 Animation - Pop-up Threat', 'Color', 'w', ...
    'Position', [100, 100, 1000, 700]);
axis_len = 0.5;
lw = 2;

% Draw complete path
plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'b:', 'LineWidth', 0.5); hold on;

% Original plan
h_old_plan = plot3(old_waypoints(1,:), old_waypoints(2,:), old_waypoints(3,:), ...
    'r--', 'LineWidth', 1.5, 'DisplayName', 'Original Plan');

% New plan (hidden initially)
h_new_plan = plot3(nan, nan, nan, 'g-', 'LineWidth', 2, 'DisplayName', 'New Plan');

% Goal
plot3(goal(1), goal(2), goal(3), 'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

% Threat (hidden initially)
[X, Y, Z] = sphere(25);
X_threat = popup_radius * X + popup_pos(1);
Y_threat = popup_radius * Y + popup_pos(2);
Z_threat = popup_radius * Z + popup_pos(3);
h_threat = surf(X_threat, Y_threat, Z_threat, ...
    'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0, 'EdgeColor', 'none', 'Visible', 'off');

grid on; axis equal;

h_body_x = plot3([0,0], [0,0], [0,0], 'r-', 'LineWidth', lw);
h_body_y = plot3([0,0], [0,0], [0,0], 'g-', 'LineWidth', lw);
h_body_z = plot3([0,0], [0,0], [0,0], 'b-', 'LineWidth', lw);
h_trail = plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'c-', 'LineWidth', 2);

set(gca, 'ZDir', 'reverse');
xlabel('North'); ylabel('East'); zlabel('Down');
h_title = title('Time: 0.00 s - No Threats');
view(40, 20);

xlim([-2, 18]); ylim([-6, 6]); zlim([-5, 1]);

% Setup video writer if saving
if save_video
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 95;
    open(v);
    fprintf('Recording animation to %s...\n', video_filename);
end

steps_skip = 30;
trail_x = []; trail_y = []; trail_z = [];
threat_shown = false;
replan_shown = false;
frame_count = 0;

for k = 1:steps_skip:n_steps
    t_curr = log.t(k);
    pos = log.x(:,k);
    q_k = log.q(:,k);
    R = quat_to_R(q_k);

    % Show threat when time comes
    if t_curr >= popup_time && ~threat_shown
        set(h_threat, 'Visible', 'on', 'FaceAlpha', 0.6);
        threat_shown = true;
    end

    % Show new plan when replanning happens
    if t_curr >= replan_time && replan_time > 0 && ~replan_shown
        set(h_new_plan, 'XData', new_waypoints(1,:), ...
            'YData', new_waypoints(2,:), 'ZData', new_waypoints(3,:));
        set(h_old_plan, 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
        replan_shown = true;
    end

    % Update trail
    trail_x = [trail_x, pos(1)];
    trail_y = [trail_y, pos(2)];
    trail_z = [trail_z, pos(3)];
    set(h_trail, 'XData', trail_x, 'YData', trail_y, 'ZData', trail_z);

    tip_x = pos + R(:,1) * axis_len;
    tip_y = pos + R(:,2) * axis_len;
    tip_z = pos + R(:,3) * axis_len;

    set(h_body_x, 'XData', [pos(1), tip_x(1)], 'YData', [pos(2), tip_x(2)], 'ZData', [pos(3), tip_x(3)]);
    set(h_body_y, 'XData', [pos(1), tip_y(1)], 'YData', [pos(2), tip_y(2)], 'ZData', [pos(3), tip_y(3)]);
    set(h_body_z, 'XData', [pos(1), tip_z(1)], 'YData', [pos(2), tip_z(2)], 'ZData', [pos(3), tip_z(3)]);

    % Update title
    dist_to_goal = norm(pos - goal);
    if t_curr < popup_time
        status = 'No Threats';
    elseif t_curr < replan_time || replan_time < 0
        status = 'THREAT DETECTED - Replanning...';
    else
        status = 'Following New Path';
    end
    set(h_title, 'String', sprintf('Time: %.2f s | Goal: %.2f m | %s', ...
        t_curr, dist_to_goal, status));

    drawnow;

    if save_video
        frame = getframe(fig);
        writeVideo(v, frame);
        frame_count = frame_count + 1;
    else
        pause(0.02);
    end
end

% Final frame
set(h_title, 'String', sprintf('Mission Complete! Final distance to goal: %.2f m', ...
    norm(log.x(:,end) - goal)));

if save_video
    % Capture final frame
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
    close(v);
    fprintf('Video saved: %s (%d frames)\n', video_filename, frame_count + 1);
end
end
