%% DEMO_MISSION_2 - Static Obstacle Avoidance
%
% Task 2: Navigate from start to goal while avoiding static obstacles.
% This demonstrates RRT-based path planning:
%   1. Define obstacles in the environment
%   2. Use RRT planner to find collision-free path
%   3. Smooth waypoints into flyable trajectory
%   4. Execute trajectory with existing control system
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

%% Define Mission with Obstacles
% Start position
start_pos = [0; 0; -2];  % Starting at 2m altitude

% Goal position
goal_pos = [12; 6; -2];  % Target waypoint

% Define obstacles (spheres blocking direct path)
obstacles = {};

% Obstacle 1: Large sphere in the middle of direct path
obstacles{1} = struct(...
    'type', 'sphere', ...
    'center', [5; 2; -2], ...
    'radius', 1.5, ...
    'active', true);

% Obstacle 2: Another obstacle forcing path deviation
obstacles{2} = struct(...
    'type', 'sphere', ...
    'center', [8; 4; -2.5], ...
    'radius', 1.2, ...
    'active', true);

% Obstacle 3: Obstacle near goal area
obstacles{3} = struct(...
    'type', 'sphere', ...
    'center', [10; 5; -1.5], ...
    'radius', 0.8, ...
    'active', true);

% Workspace bounds for RRT
bounds = struct(...
    'min', [-2; -2; -6], ...
    'max', [15; 10; 0]);

fprintf('=== Mission 2: Static Obstacle Avoidance ===\n');
fprintf('Start:     [%.1f, %.1f, %.1f] m\n', start_pos);
fprintf('Goal:      [%.1f, %.1f, %.1f] m\n', goal_pos);
fprintf('Obstacles: %d\n', length(obstacles));
for i = 1:length(obstacles)
    obs = obstacles{i};
    fprintf('  Obstacle %d: %s at [%.1f, %.1f, %.1f], radius %.1f m\n', ...
        i, obs.type, obs.center, obs.radius);
end
fprintf('\n');

%% Setup Waypoint Manager with Obstacles
wm = waypoint_manager();

% Configure RRT* parameters (use RRT* for better path quality in static environments)
wm.rrt_params.algorithm = 'rrt_star';  % 'rrt' or 'rrt_star'
wm.rrt_params.max_iter = 8000;
wm.rrt_params.step_size = 0.4;
wm.rrt_params.goal_bias = 0.15;
wm.rrt_params.safety_margin = 0.5;    % Stay 0.5m away from obstacles
wm.rrt_params.rewire_radius = 1.5;    % RRT* rewiring neighborhood

wm.set_mission(start_pos, goal_pos, obstacles, bounds);

% Plan path using RRT/RRT*
fprintf('Planning path with %s...\n', upper(wm.rrt_params.algorithm));
tic;
[success, waypoints] = wm.plan_path();
planning_time = toc;

if ~success
    error('RRT path planning failed! Try increasing max_iter or adjusting parameters.');
end

fprintf('Path planning completed in %.2f seconds\n', planning_time);
fprintf('Waypoints: %d\n', size(waypoints, 2));

% Calculate path length
path_length = 0;
for i = 1:size(waypoints, 2)-1
    path_length = path_length + norm(waypoints(:,i+1) - waypoints(:,i));
end
direct_distance = norm(goal_pos - start_pos);
fprintf('Direct distance: %.2f m\n', direct_distance);
fprintf('Path length:     %.2f m (%.1f%% longer)\n', path_length, (path_length/direct_distance - 1)*100);

% Mission duration based on path length
mission_speed = 2.0;  % m/s average speed
Tf = max(8.0, path_length / mission_speed + 4.0);
fprintf('Duration: %.1f s\n\n', Tf);

% Generate smooth trajectory
traj = wm.generate_trajectory(Tf);

%% Initialize Simulation State
N = round(Tf/dt);
[state, motors, filters] = init_sim_state(params, start_pos);

% Initialize logs
log = init_log(N);

%% Main Simulation Loop
fprintf('Running simulation...\n');

for k = 1:N
    t = (k-1)*dt;

    % Get reference from waypoint manager
    ref = wm.get_reference(t);

    % Create state_ref struct for yaw_planner (expects .v and .a fields)
    state_ref.v = ref.vel;
    state_ref.a = ref.acc;

    % Add yaw planning (tangent mode - point in direction of travel)
    [ref.psi, ref.psi_dot, ref.psi_ddot] = yaw_planner(t, 'tangent', state_ref, params);

    % Run simulation step
    [state, motors, filters, out] = sim_step(state, motors, filters, ref, params, dt);

    % Check for collision (safety monitoring)
    [is_collision, min_dist, ~] = collision_checker(state.x, obstacles, 0.1);
    if is_collision
        fprintf('[%.2fs] WARNING: Collision detected! Min distance: %.3f m\n', t, min_dist);
    end

    % Log data
    log = log_step(log, k, t, state, ref, out);
end

%% Results
fprintf('\n=== Results ===\n');
final_pos = state.x;
final_error = norm(final_pos - goal_pos);
fprintf('Final position: [%.3f, %.3f, %.3f] m\n', final_pos);
fprintf('Final error:    %.3f m\n', final_error);

% Check minimum clearance from obstacles during flight
min_clearance = inf;
for k = 1:size(log.x, 2)
    [~, dist, ~] = collision_checker(log.x(:,k), obstacles, 0);
    min_clearance = min(min_clearance, dist);
end
fprintf('Minimum obstacle clearance: %.3f m\n', min_clearance);

if final_error < 0.5 && min_clearance > 0
    fprintf('Mission SUCCESS: Goal reached without collision!\n');
elseif min_clearance <= 0
    fprintf('Mission FAILED: Collision occurred!\n');
else
    fprintf('Mission INCOMPLETE: Did not reach goal within tolerance.\n');
end

%% Plots
plot_mission_results(log, params, goal_pos, obstacles, waypoints, 'Mission 2: Static Obstacle Avoidance');

%% Animation
animate_mission(log, goal_pos, obstacles, waypoints);

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

function plot_mission_results(log, params, goal, obstacles, waypoints, title_str)
    t = log.t;

    % 3D Trajectory with obstacles
    figure('Name', title_str, 'Position', [100, 100, 900, 700]);

    % Draw trajectory
    plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Actual Path');
    hold on;
    plot3(log.x_ref(1,:), log.x_ref(2,:), log.x_ref(3,:), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Reference');

    % Draw waypoints
    plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ko-', 'MarkerSize', 8, ...
          'MarkerFaceColor', 'y', 'LineWidth', 1, 'DisplayName', 'RRT Waypoints');

    % Start and goal
    plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'gs', 'MarkerSize', 15, ...
          'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot3(goal(1), goal(2), goal(3), 'r^', 'MarkerSize', 15, ...
          'MarkerFaceColor', 'r', 'DisplayName', 'Goal');

    % Draw obstacles
    for i = 1:length(obstacles)
        obs = obstacles{i};
        if strcmp(obs.type, 'sphere')
            [X, Y, Z] = sphere(30);
            X = obs.radius * X + obs.center(1);
            Y = obs.radius * Y + obs.center(2);
            Z = obs.radius * Z + obs.center(3);
            surf(X, Y, Z, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.6, ...
                 'EdgeColor', 'none');
        end
    end

    xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');
    title(title_str);
    legend('Location', 'best');
    set(gca, 'ZDir', 'reverse');
    grid on; axis equal;
    view(30, 25);

    % Top-down view (N-E plane)
    figure('Name', [title_str ' - Top View']);
    plot(log.x(1,:), log.x(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Actual');
    hold on;
    plot(waypoints(1,:), waypoints(2,:), 'ko-', 'MarkerSize', 8, ...
         'MarkerFaceColor', 'y', 'DisplayName', 'Waypoints');
    plot(log.x(1,1), log.x(2,1), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(goal(1), goal(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

    % Draw obstacle circles (top view)
    theta_circle = linspace(0, 2*pi, 50);
    for i = 1:length(obstacles)
        obs = obstacles{i};
        if strcmp(obs.type, 'sphere')
            x_circle = obs.center(1) + obs.radius * cos(theta_circle);
            y_circle = obs.center(2) + obs.radius * sin(theta_circle);
            fill(x_circle, y_circle, [1 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);
        end
    end

    xlabel('North [m]'); ylabel('East [m]');
    title([title_str ' - Top View (N-E Plane)']);
    legend('Location', 'best');
    grid on; axis equal;

    % Position error
    position_error = vecnorm(log.x - log.x_ref, 2, 1);
    figure('Name', [title_str ' - Error']);
    plot(t, position_error, 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('Position Error [m]');
    title('Position Tracking Error');
    grid on;

    [peak_err, peak_idx] = max(position_error);
    fprintf('Peak position error: %.3f m at t=%.3f s\n', peak_err, t(peak_idx));
end

function animate_mission(log, goal, obstacles, waypoints, save_video, video_filename)
    % ANIMATE_MISSION Animate drone flight with optional MP4 saving
    %
    % Inputs:
    %   save_video    - (optional) true to save as MP4 (default: true)
    %   video_filename - (optional) Output filename (default: 'mission_2_animation.mp4')

    if nargin < 5
        save_video = true;
    end
    if nargin < 6
        video_filename = 'mission_2_animation.mp4';
    end

    n_steps = numel(log.t);

    fig = figure('Name', 'Mission 2 Animation', 'Color', 'w', 'Position', [100, 100, 900, 700]);
    axis_len = 0.6;
    lw = 2;

    % Draw full trajectory path
    plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'b:', 'LineWidth', 1); hold on;

    % Draw waypoints
    plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ko-', 'MarkerSize', 6, ...
          'MarkerFaceColor', 'y', 'LineWidth', 1);

    % Draw obstacles
    for i = 1:length(obstacles)
        obs = obstacles{i};
        if strcmp(obs.type, 'sphere')
            [X, Y, Z] = sphere(25);
            X = obs.radius * X + obs.center(1);
            Y = obs.radius * Y + obs.center(2);
            Z = obs.radius * Z + obs.center(3);
            surf(X, Y, Z, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
    end

    % Goal
    plot3(goal(1), goal(2), goal(3), 'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

    grid on; axis equal;

    h_body_x = plot3([0,0], [0,0], [0,0], 'r-', 'LineWidth', lw);
    h_body_y = plot3([0,0], [0,0], [0,0], 'g-', 'LineWidth', lw);
    h_body_z = plot3([0,0], [0,0], [0,0], 'b-', 'LineWidth', lw);
    h_trail = plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'c-', 'LineWidth', 1.5);

    set(gca, 'ZDir', 'reverse');
    xlabel('North'); ylabel('East'); zlabel('Down');
    h_title = title('Time: 0.00 s');
    view(45, 25);

    % Setup video writer if saving
    if save_video
        v = VideoWriter(video_filename, 'MPEG-4');
        v.FrameRate = 30;
        v.Quality = 95;
        open(v);
        fprintf('Recording animation to %s...\n', video_filename);
    end

    steps_skip = 40;
    trail_x = []; trail_y = []; trail_z = [];
    frame_count = 0;

    for k = 1:steps_skip:n_steps
        pos = log.x(:,k);
        q_k = log.q(:,k);
        R = quat_to_R(q_k);

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

        dist_to_goal = norm(pos - goal);
        [~, min_obs_dist, ~] = collision_checker(pos, obstacles, 0);
        set(h_title, 'String', sprintf('Time: %.2f s | Goal: %.2f m | Obstacle: %.2f m', ...
            log.t(k), dist_to_goal, min_obs_dist));

        drawnow;

        if save_video
            frame = getframe(fig);
            writeVideo(v, frame);
            frame_count = frame_count + 1;
        else
            pause(0.03);
        end
    end

    if save_video
        close(v);
        fprintf('Video saved: %s (%d frames)\n', video_filename, frame_count);
    end
end
