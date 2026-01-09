%function demo_indi_6dof_ned()
%DEMO_INDI_6DOF_NED
%   Full 6-DOF quadrotor with:
%     - NED frame (x=[N;E;D])
%     - INDI linear acceleration control (Eqs. 16–21)
%     - Simple quaternion attitude PD
%     - Ideal motors (no dynamics, only G1 mapping)
clear; clc; close all;

params = quad_params_indi();
g   = params.g;
m   = params.m;
J   = params.J;
e3  = [0;0;1];
tau_f = params.tau_f;

dt = 0.002;       % [s]
Tf = 10.0;        % Default simulation time [s] (will be overridden for vloop)
N  = round(Tf/dt);

% State
[state, motors, filters] = init_sim_state(params);

%% Min Snap Poly Traj
% figure;
% empty = zeros(2,3);
% plot(empty(1,:), empty(2,:), 'ro', 'MarkerSize', 10, 'DisplayName', 'Waypoints');
% hold on; grid on; xlabel('North [m]'); ylabel('Down [m]');
% axis([-10, 10, -10, 10])
% title('Draw a polyline. Double-click to finish.');
% 
% h = drawpolyline('Color','b');        % user draws and double-clicks to finish
% pos = h.Position;                     % M-by-2 array of [x y] points
% 
% % waypoints = [pos(:,1)'; pos(:,2)'; zeros(1,size(pos,1))];
% waypoints = [pos(:,1)';zeros(1,size(pos,1)) ; pos(:,2)'];
% timePoints = linspace(0, Tf, size(pos,1));
% numSamples = N; % Number of samples to generate

% Trial for barrel roll waypoints
% waypoints = [0, 0, 0;
%             0, 2, 0; 
%             1, 4, 1; 
%             0, 5, 2;
%             -1, 6, 1;
%             0, 7, 0;
%             0, 9, 0]'; % Define waypoints for the trajectory
% timePoints = [0,1,2,3,4,5,6]; % Define time points corresponding to the waypoints
% numSamples = N; % Number of samples to generate
% [xref,vref,aref,jref,sref] = minsnappolytraj(waypoints, timePoints, numSamples);

% % Geometry
% R = 5;   % Loop radius [m]
% 
% % Waypoints: n-by-p (3-by-5)
% % Order: 0°, 90°, 180°, 270°, 360°
% waypoints = [ ...
%      0,   R,   0,  -R,   0 ;    % x
%      0,   0,   0,   0,   0 ;    % y (fixed)
%      0,   -R, -2*R,   -R,   0 ];   % z
% 
% % Time allocation (non-uniform)
% timePoints = [0, 1.0, 2.4, 3.8, 5.0];   % [s]
% 
% % Tangent velocity definition
% theta = deg2rad([0, 90, 180, 270, 360]);
% v     = [6, 5, 3, 5, 6];   % Speed profile [m/s]
% 
% vx = v .* cos(theta);
% vy = zeros(size(v));      % y velocity = 0
% vz = - v .* sin(theta);
% 
% % Velocity BC: n-by-p (3-by-5)
% VelocityBoundaryCondition = [ ...
%     vx ;
%     vy ;
%     vz ];
% 
% % Acceleration Boundary Condition
% AccelerationBoundaryCondition = nan(3,5);
% 
% % Enforce zero acceleration at key points
% AccelerationBoundaryCondition(:,1) = [0; 0; 0];   % start
% AccelerationBoundaryCondition(:,3) = [0; 0; 0];   % top
% AccelerationBoundaryCondition(:,5) = [0; 0; 0];   % end
% 
% % Trajectory generation
% numSamples = N;
% 
% [xref,vref,aref,jref,sref] = minsnappolytraj( ...
%     waypoints, timePoints, numSamples, ...
%     VelocityBoundaryCondition = VelocityBoundaryCondition, ...
%     AccelerationBoundaryCondition = AccelerationBoundaryCondition );




ref_params = params;
ref_params.shape = "barrel_roll";
ref_params.yaw = "constant";

% Initialize logs (must be done after N is finalized)
log = init_log(N);

for k = 1:N
    t = (k-1)*dt;
    ref.x = xref(:,k);
    ref.v = vref(:,k);
    ref.a = aref(:,k);
    ref.j = jref(:,k);
    ref.s = sref(:,k);
    ref.psi = 0; % Assuming constant yaw for simplicity
    ref.psi_dot = 0; % Assuming constant yaw rate for simplicity
    ref.psi_ddot = 0; % Assuming constant yaw acceleration for simplicity
    %ref = reference_flat_outputs(t, ref_params);
    [state, motors, filters, out] = sim_step(state, motors, filters, ref, params, dt);
    log = log_step(log, k, t, state, ref, out);
end

% Keep final state variables available in the workspace.
x = state.x;
v = state.v;
q = state.q;
omega = state.omega;
omega_mot = motors.omega;
a_f = filters.a_f;
tau_bz_f = filters.tau_bz_f;
%% Plots
plot_results(log, params);

%% DRONE ANIMATION
animate_flight(log);

%% Export Data to CSV for Unity
% Format: Time, N, E, D, Qw, Qx, Qy, Qz
% Not: Unity sol el kuralini kullanir, donusum C# tarafinda yapacagiz.
export_filename = 'C:/Users/tayfu/OneDrive/Belgeler/Makale ve Kitaplar/Uçak Mühendisliği Kaynakları/TAI/Tal&Karaman_2021/DifferentialFlatness/Assets/Resources/drone_sim_data.csv';
export_log_to_csv(log, export_filename);
%%

%end
function [state, motors, filters] = init_sim_state(params)
%INIT_SIM_STATE  Initialize state, motors, and filters.
state.x = [0;0;0];
state.v = [0;0;0];
state.q = [1;0;0;0];
state.omega = [0;0;0];

motors.omega = params.omega_hover * ones(4,1);

filters.a_f = zeros(3,1);
filters.tau_bz_f = zeros(3,1);
end

function log = init_log(N)
%INIT_LOG  Preallocate log arrays.
log.t               = zeros(1,N);
log.x               = zeros(3,N);
log.x_ref           = zeros(3,N);

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

function ref = reference_flat_outputs(t, params)
%REFERENCE_FLAT_OUTPUTS  Wrap flat_outputs_demo output in a struct.
[ref.x, ref.v, ref.a, ref.j, ref.s, ref.psi, ref.psi_dot, ref.psi_ddot] = flat_outputs_demo(t, params);
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


% INDI Linear Acceleration and Yaw Control
a_b = R' * (a_true - g*e3);
a_meas = R*a_b + g*e3 + 0.1*rand(3,1);
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

err_eq19 = a_true - (tau_bz - filters.tau_bz_f + filters.a_f);

out.a_true        = a_true;
out.a_c           = a_c;
out.a_f           = filters.a_f;
out.a_ref         = ref.a;
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
log.x_ref(:,k)    = ref.x;
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

function plot_results(log, params)
t = log.t;
n_steps = numel(t);

figure;
subplot(3,1,1);
plot(t, log.x(1,:), t, log.x_ref(1,:),'--');
ylabel('N [m]'); legend('N','N_{ref}');
title('North tracking');

subplot(3,1,2);
plot(t, log.x(2,:), t, log.x_ref(2,:),'--');
ylabel('E [m]'); legend('E','E_{ref}');
title('East tracking');

subplot(3,1,3);
plot(t, log.x(3,:), t, log.x_ref(3,:),'--');
ylabel('D [m]'); legend('D','D_{ref}');
xlabel('t [s]');
title('Down tracking (NED, down positive)');

figure;
plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'b-', 'DisplayName', 'Actual Trajectory');
hold on;
plot3(log.x_ref(1,:), log.x_ref(2,:), log.x_ref(3,:), 'r--', 'DisplayName', 'Reference Trajectory');
xlabel('N [m]');
ylabel('E [m]');
zlabel('D [m]');
title('3D Trajectory of the Quadrotor');
legend show;
set(gca, 'ZDir', 'reverse');
xlabel('North'); ylabel('East'); zlabel('Down');
grid on; axis equal;

axis equal;
xlim([-10, 10]);
ylim([-10, 10]);
zlim([-10, 10]);

figure(3);
hold on;
draw_step = 300;

for k = 1:draw_step:n_steps
    pos = log.x(:,k);
    q_k = log.q(:,k);
    R_k = quat_to_R(q_k);

    axis_len = 0.5;
    bx = R_k(:,1) * axis_len;
    by = R_k(:,2) * axis_len;
    bz = R_k(:,3) * axis_len;

    quiver3(pos(1), pos(2), pos(3), bx(1), bx(2), bx(3), 'r-', 'LineWidth', 2);
    quiver3(pos(1), pos(2), pos(3), by(1), by(2), by(3), 'g-', 'LineWidth', 2);
    quiver3(pos(1), pos(2), pos(3), bz(1), bz(2), bz(3), 'b-', 'LineWidth', 2);
end
title('3D Trajectory with Attitude (Red=Front)');

position_error = vecnorm(log.x - log.x_ref, 2, 1);
figure;
plot(t, position_error, 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Position Error Norm [m]');
title('Position Norm Error Over Time');
grid on;
[peak_err, peak_idx] = max(position_error);
fprintf('Peak position error: %.3f m at t=%.3f s\n', peak_err, t(peak_idx));

figure;
subplot(3,1,1);
plot(t, log.omega(1,:), t, log.omega_ref(1,:), t, log.omega_dot_ref(1,:),'--');
ylabel('p rotations');
xlabel('t [s]');
legend('p', 'p_{ref}', '\dot{p}');

subplot(3,1,2);
plot(t, log.omega(2,:), t, log.omega_ref(2,:), t, log.omega_dot_ref(2,:),'--');
ylabel('q rotations');
xlabel('t [s]');
legend('q', 'q_{ref}', '\dot{q}');

subplot(3,1,3);
plot(t, log.omega(3,:), t, log.omega_ref(3,:), t, log.omega_dot_ref(3,:),'--');
ylabel('r rotations');
xlabel('t [s]');
legend('r', 'r_{ref}', '\dot{r}');

figure;
subplot(3,1,1);
plot(t, log.alpha_cmd(1,:));
ylabel('x [rad/s^2]');
title('Alpha Command - X Element');

subplot(3,1,2);
plot(t, log.alpha_cmd(2,:));
ylabel('y [rad/s^2]');
title('Alpha Command - Y Element');

subplot(3,1,3);
plot(t, log.alpha_cmd(3,:));
ylabel('z [rad/s^2]');
title('Alpha Command - Z Element');
xlabel('t [s]');

figure;
plot(t, log.xi_e(1,:), 'DisplayName', 'xi_e_x');
hold on;
plot(t, log.xi_e(2,:), 'DisplayName', 'xi_e_y');
plot(t, log.xi_e(3,:), 'DisplayName', 'xi_e_z');
hold off;
xlabel('Time [s]');
ylabel('xi_e values');
title('xi_e Components Over Time');
legend show;

if isfield(log, 'tau_norm')
    figure;
    subplot(2,1,1);
    plot(t, log.tau_norm, 'LineWidth', 1.5);
    ylabel('||\tau_{bz,c}|| [m/s^2]');
    title('Specific Thrust Command Norm');
    grid on;
    subplot(2,1,2);
    plot(t, log.T_cmd, 'LineWidth', 1.5);
    ylabel('T_{cmd} [N]');
    xlabel('t [s]');
    title('Thrust Command');
    grid on;
end

if isfield(log, 'omega_mot_cmd') || isfield(log, 'omega_mot')
    figure;

    % Create a 2x1 subplot layout
    subplot(2, 1, 1);
    if isfield(log, 'omega_mot_cmd')
        plot(t, log.omega_mot_cmd);
        xlabel('t [s]');
        ylabel('\omega_{mot,cmd} [rad/s]');
        title('Motor Speed Commands');
        legend('m1','m2','m3','m4');
        grid on;

        sat_mask = log.omega_mot_cmd >= params.omega_max * 0.999 | ...
            log.omega_mot_cmd <= params.omega_min + 1e-6;
        if any(sat_mask, 'all')
            first_idx = find(any(sat_mask, 1), 1, 'first');
            fprintf('Motor saturation begins at t=%.3f s\n', t(first_idx));
        else
            fprintf('No motor saturation detected.\n');
        end

        max_cmd = max(log.omega_mot_cmd, [], 'all');
        min_cmd = min(log.omega_mot_cmd, [], 'all');
        fprintf('Motor cmd range: [%.1f, %.1f] rad/s (limits [%.1f, %.1f])\n', ...
            min_cmd, max_cmd, params.omega_min, params.omega_max);
    end

    % Plot motor speeds if available
    if isfield(log, 'omega_mot')
        subplot(2, 1, 2);
        plot(t, log.omega_mot);
        xlabel('t [s]');
        ylabel('\omega_{mot} [rad/s]');
        title('Motor Speeds');
        legend('m1','m2','m3','m4');
        grid on;
    end
end

if isfield(log, 'T_cmd')
    T_max = 4 * params.kT * params.omega_max^2;
    fprintf('Thrust cmd max: %.1f N (T_max=%.1f N)\n', max(log.T_cmd), T_max);
end
end

function animate_flight(log)
n_steps = numel(log.t);

figure('Name', 'Drone Flight Animation', 'Color', 'w');
axis_len = 0.8;
lw = 2;

plot3(log.x(1,:), log.x(2,:), log.x(3,:), 'k:', 'LineWidth', 1); hold on;
grid on; axis equal;
xlabel('North'); ylabel('East'); zlabel('Down');
view(3);

h_body_x = plot3([0,0], [0,0], [0,0], 'r-', 'LineWidth', lw);
h_body_y = plot3([0,0], [0,0], [0,0], 'g-', 'LineWidth', lw);
h_body_z = plot3([0,0], [0,0], [0,0], 'b-', 'LineWidth', lw);

set(gca, 'ZDir', 'reverse');
xlabel('North'); ylabel('East'); zlabel('Down');
grid on; axis equal;

h_title = title('Time: 0.00 s');

steps_skip = 50;

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

    set(h_title, 'String', sprintf('Time: %.2f s | Yaw: %.1f deg', log.t(k), rad2deg(atan2(R(2,1), R(1,1)))));

    drawnow limitrate;
    pause(0.1);
end
end

function export_log_to_csv(log, filename)
data_to_export = [log.t', log.x', log.q'];
writematrix(data_to_export, filename);
fprintf('Data saved to "%s" file for Unity.\n', filename);
end
function [phi, theta, psi] = rotm_to_euler_zyx(R)
%ROTM_TO_EULER_ZYX  ZYX (yaw-pitch-roll) Euler angles from rotation matrix.
%   R: body->inertial rotation (as in Tal & Karaman)
%   Returns phi (roll), theta (pitch), psi (yaw).

% Protect against numerical noise
R = R ./ max(det(R)^(1/3), 1e-9);

theta = asin(-R(3,1));
phi   = atan2(R(3,2), R(3,3));
psi   = atan2(R(2,1), R(1,1));
end

% function S = euler_rate_matrix(phi, theta)
% %EULER_RATE_MATRIX  Matrix S such that xi_dot = S * Omega.
% %   xi_dot = [phi_dot; theta_dot; psi_dot]
% %   Omega  = [p; q; r]
% 
% S = [ 1,  sin(phi)*tan(theta),  cos(phi)*tan(theta);
%       0,  cos(phi),           -sin(phi);
%       0,  sin(phi)/cos(theta), 0];
% end


function q_inc = compute_incremental_attitude_cmd(R_curr, tau_bz_c, psi_ref)
%COMPUTE_INCREMENTAL_ATTITUDE_CMD
%   Tal & Karaman (2021) - Section III-B, Eq. 22-26
%   Computes the "incremental" quaternion that will provide the desired thrust vector (tau_bz_c)
%   based on the current orientation (R_curr).
%
%   Output: q_inc (Rotation command to be applied in the Body Frame)

    % 1. Desired Thrust Direction (Inertial Frame)
    % The thrust vector (tau_bz_c) is typically in the [N, E, D] axes.
    % The drone's thrust is in the opposite direction of the Z axis (-bz).
    t_norm = norm(tau_bz_c);
    if t_norm < 1e-6
        t_des_in = [0;0;-1]; % Hover (Upward thrust)
    else
        t_des_in = tau_bz_c / t_norm; % Acceleration direction
    end
    
    % The thrust direction must be the OPPOSITE of the Body Z axis (-b_z).
    % That is, the target Body Z axis (in the Inertial frame):
    bz_des_in = -t_des_in; 

    % 2. Transform the target Body Z to the CURRENT Body Frame (according to Eq. 22)
    % bz_des_body = R_curr' * bz_des_in
    bz_des_body = R_curr' * bz_des_in;

    % 3. Tilt Rotation (Eq. 23 - Vector Alignment)
    % The difference between the current Body Z ([0;0;1]) and the target Body Z (bz_des_body).
    % This operation provides the "Tilt" correction.
    
    current_z = [0;0;1];
    
    % Quaternion for rotation between two vectors (Shortest Arc)
    % q = [dot + sqrt(len1*len2), cross]
    cross_prod = cross(current_z, bz_des_body);
    dot_prod   = dot(current_z, bz_des_body);
    
    % Singularity protection (if looking in the exact opposite direction)
    if dot_prod < -0.9999
        % 180 degree rotation (around the X axis)
        q_tilt = [0; 1; 0; 0]; 
    else
        s = sqrt(2 * (1 + dot_prod));
        q_tilt = [0.5 * s; cross_prod / s];
    end
    q_tilt = quat_normalize(q_tilt);

    % 4. Yaw (Drift) Correction (Eq. 24-25)
    % The intermediate orientation after applying the tilt correction
    R_tilt = quat_to_R(q_tilt);
    
    % Where is the X axis looking in this intermediate orientation?
    bx_inter = R_tilt(:,1); % New X in the Body frame
    
    % How much more do we need to turn in the "Body Frame" according to the desired Yaw angle (psi_ref)?
    
    % Current Heading (Inertial)
    % The operation R_curr * R_tilt gives the new stance in the Inertial frame.
    R_new_in = R_curr * R_tilt;
    
    % Current Yaw angle (from R_new_in)
    psi_curr = atan2(R_new_in(2,1), R_new_in(1,1));
    
    % Yaw Error (via the shortest path)
    psi_err = psi_ref - psi_curr;
    while psi_err > pi,  psi_err = psi_err - 2*pi; end
    while psi_err < -pi, psi_err = psi_err + 2*pi; end
    
    % Yaw rotation (around the Z axis)
    q_yaw = [cos(psi_err/2); 0; 0; sin(psi_err/2)];
    
    % 5. Total Incremental Command (Eq. 26)
    % First Tilt, then Yaw (or the order may change according to frame definition,
    % here we are doing body-intrinsic rotation: q_inc = q_tilt * q_yaw)
    
    q_inc = quat_mul(q_tilt, q_yaw);
    
    % Make the scalar part of q_inc positive (to prevent singularity)
    if q_inc(1) < 0
        q_inc = -q_inc;
    end
end

function R = quat_to_R_custom(q)
%QUAT_TO_R_CUSTOM  Quaternion to Rotation Matrix
%   q = [qw; qx; qy; qz] (Scalar first)
%   R = Body to Inertial rotation matrix

    qw = q(1); qx = q(2); qy = q(3); qz = q(4);

    R = [1 - 2*qy^2 - 2*qz^2,    2*qx*qy - 2*qz*qw,    2*qx*qz + 2*qy*qw;
         2*qx*qy + 2*qz*qw,    1 - 2*qx^2 - 2*qz^2,    2*qy*qz - 2*qx*qw;
         2*qx*qz - 2*qy*qw,    2*qy*qz + 2*qx*qw,    1 - 2*qx^2 - 2*qy^2];
end
