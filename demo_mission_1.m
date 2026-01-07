%% DEMO_MISSION_1 - Point-to-Point Navigation (No Obstacles)
%
% Task 1: Navigate from start position to a target point with no obstacles.
% This demonstrates the basic motion planning pipeline:
%   1. Define start and goal positions
%   2. Generate smooth trajectory using waypoint_manager
%   3. Execute trajectory with existing control system
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

%% Define Mission
% Start position (where drone initializes)
start_pos = [0; 0; 0];  % [N; E; D] in meters

% Goal position (target waypoint)
goal_pos = [10; 5; -3];  % 10m North, 5m East, 3m Up (D=-3 in NED)

% Mission duration (will be auto-calculated based on distance)
mission_speed = 2.0;  % m/s average speed
distance = norm(goal_pos - start_pos);
Tf = max(5.0, distance / mission_speed + 4.0);  % Add settling time

fprintf('=== Mission 1: Point-to-Point Navigation ===\n');
fprintf('Start:    [%.1f, %.1f, %.1f] m\n', start_pos);
fprintf('Goal:     [%.1f, %.1f, %.1f] m\n', goal_pos);
fprintf('Distance: %.2f m\n', distance);
fprintf('Duration: %.1f s\n\n', Tf);

%% Setup Waypoint Manager (No Obstacles)
wm = waypoint_manager();
wm.set_mission(start_pos, goal_pos, {}, []);  % Empty obstacles

% Plan path (direct line for no obstacles)
[success, waypoints] = wm.plan_path();
if ~success
    error('Path planning failed!');
end

% Generate smooth trajectory
traj = wm.generate_trajectory(Tf);

fprintf('Planned waypoints: %d\n', size(waypoints, 2));
fprintf('Trajectory samples: %d\n\n', length(traj.t));

%% Initialize Simulation State
N = round(Tf/dt);
[state, motors, filters] = init_sim_state(params);

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

    % Add yaw planning (constant yaw for this mission)
    [ref.psi, ref.psi_dot, ref.psi_ddot] = yaw_planner(t, 'constant', state_ref, params);

    % Run simulation step (same control as original demo)
    [state, motors, filters, out] = sim_step(state, motors, filters, ref, params, dt);

    % Log data
    log = log_step(log, k, t, state, ref, out);
end

%% Results
fprintf('\n=== Results ===\n');
final_pos = state.x;
final_error = norm(final_pos - goal_pos);
fprintf('Final position: [%.3f, %.3f, %.3f] m\n', final_pos);
fprintf('Final error:    %.3f m\n', final_error);

if final_error < 0.5
    fprintf('Mission SUCCESS: Goal reached within tolerance!\n');
else
    fprintf('Mission INCOMPLETE: Did not reach goal within tolerance.\n');
end

%% Plots
plot_mission_results(log, params, goal_pos, {}, 'Mission 1: Point-to-Point');

%% Animation
animate_mission(log, goal_pos, {});

%% ==================== LOCAL FUNCTIONS ====================

function [state, motors, filters] = init_sim_state(params)
    state.x = [0;0;0];
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

function plot_mission_results(log, params, goal, obstacles, title_str)
    t = log.t;

    % 3D Trajectory
    figure('Name', title_str);
    plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Actual');
    hold on;
    plot3(log.x_ref(1,:), log.x_ref(2,:), log.x_ref(3,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reference');

    % Start and goal markers
    plot3(log.x(1,1), log.x(2,1), log.x(3,1), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot3(goal(1), goal(2), goal(3), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Goal');

    % Draw obstacles
    for i = 1:length(obstacles)
        obs = obstacles{i};
        if strcmp(obs.type, 'sphere')
            [X, Y, Z] = sphere(20);
            X = obs.radius * X + obs.center(1);
            Y = obs.radius * Y + obs.center(2);
            Z = obs.radius * Z + obs.center(3);
            surf(X, Y, Z, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.5, ...
                 'EdgeColor', 'none', 'DisplayName', 'Obstacle');
        end
    end

    xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');
    title(title_str);
    legend('Location', 'best');
    set(gca, 'ZDir', 'reverse');
    grid on; axis equal;
    view(30, 20);

    % Position tracking
    figure('Name', [title_str ' - Tracking']);
    subplot(3,1,1);
    plot(t, log.x(1,:), 'b-', t, log.x_ref(1,:), 'r--');
    ylabel('N [m]'); legend('Actual', 'Ref'); title('North');

    subplot(3,1,2);
    plot(t, log.x(2,:), 'b-', t, log.x_ref(2,:), 'r--');
    ylabel('E [m]'); legend('Actual', 'Ref'); title('East');

    subplot(3,1,3);
    plot(t, log.x(3,:), 'b-', t, log.x_ref(3,:), 'r--');
    ylabel('D [m]'); legend('Actual', 'Ref'); title('Down');
    xlabel('Time [s]');

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

function animate_mission(log, goal, obstacles, save_video, video_filename)
    % ANIMATE_MISSION Animate drone flight with optional MP4 saving
    %
    % Inputs:
    %   log           - Simulation log struct
    %   goal          - [3x1] Goal position
    %   obstacles     - Cell array of obstacles
    %   save_video    - (optional) true to save as MP4 (default: true)
    %   video_filename - (optional) Output filename (default: 'mission_1_animation.mp4')

    if nargin < 4
        save_video = true;
    end
    if nargin < 5
        video_filename = 'mission_1_animation.mp4';
    end

    n_steps = numel(log.t);

    fig = figure('Name', 'Mission Animation', 'Color', 'w', 'Position', [100, 100, 800, 600]);
    axis_len = 0.8;
    lw = 2;

    % Draw trajectory path
    plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'k:', 'LineWidth', 1); hold on;

    % Draw obstacles
    for i = 1:length(obstacles)
        obs = obstacles{i};
        if strcmp(obs.type, 'sphere')
            [X, Y, Z] = sphere(20);
            X = obs.radius * X + obs.center(1);
            Y = obs.radius * Y + obs.center(2);
            Z = obs.radius * Z + obs.center(3);
            surf(X, Y, Z, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
    end

    % Goal marker
    plot3(goal(1), goal(2), goal(3), 'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

    grid on; axis equal;

    h_body_x = plot3([0,0], [0,0], [0,0], 'r-', 'LineWidth', lw);
    h_body_y = plot3([0,0], [0,0], [0,0], 'g-', 'LineWidth', lw);
    h_body_z = plot3([0,0], [0,0], [0,0], 'b-', 'LineWidth', lw);

    set(gca, 'ZDir', 'reverse');
    xlabel('North'); ylabel('East'); zlabel('Down');
    h_title = title('Time: 0.00 s');

    % Setup video writer if saving
    if save_video
        v = VideoWriter(video_filename, 'MPEG-4');
        v.FrameRate = 30;
        v.Quality = 95;
        open(v);
        fprintf('Recording animation to %s...\n', video_filename);
    end

    steps_skip = 50;
    frame_count = 0;

    for k = 1:steps_skip:n_steps
        pos = log.x(:,k);
        q_k = log.q(:,k);
        R = quat_to_R(q_k);

        tip_x = pos + R(:,1) * axis_len;
        tip_y = pos + R(:,2) * axis_len;
        tip_z = pos + R(:,3) * axis_len;

        set(h_body_x, 'XData', [pos(1), tip_x(1)], 'YData', [pos(2), tip_x(2)], 'ZData', [pos(3), tip_x(3)]);
        set(h_body_y, 'XData', [pos(1), tip_y(1)], 'YData', [pos(2), tip_y(2)], 'ZData', [pos(3), tip_y(3)]);
        set(h_body_z, 'XData', [pos(1), tip_z(1)], 'YData', [pos(2), tip_z(2)], 'ZData', [pos(3), tip_z(3)]);

        dist_to_goal = norm(pos - goal);
        set(h_title, 'String', sprintf('Time: %.2f s | Distance to goal: %.2f m', log.t(k), dist_to_goal));

        drawnow;

        if save_video
            frame = getframe(fig);
            writeVideo(v, frame);
            frame_count = frame_count + 1;
        else
            pause(0.05);
        end
    end

    if save_video
        close(v);
        fprintf('Video saved: %s (%d frames)\n', video_filename, frame_count);
    end
end
