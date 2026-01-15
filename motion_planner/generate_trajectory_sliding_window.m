function [full_traj, debug_info] = generate_trajectory_sliding_window(waypoints, V_avg, window_size)
%GENERATE_TRAJECTORY_SLIDING_WINDOW Real-time trajectory generation
%
% Generates a minimum-snap trajectory using a sliding window approach.
% This avoids solving a large QP for the entire path, enabling real-time
% replanning capabilities.
%
% Inputs:
%   waypoints   - [3 x M] matrix of waypoints
%   V_avg       - Average speed [m/s] for time allocation
%   window_size - Number of waypoints per window (default: 3)
%
% Outputs:
%   full_traj   - Struct containing concatenated trajectory data
%                 .pp  - Piecewise polynomial struct (spline)
%                 .t   - Time vector
%                 .x   - Position [3 x N]
%                 .v   - Velocity
%                 .a   - Acceleration
%
% Author: Quadrotor Terrain Following Project

if nargin < 3
    window_size = 3;
end

% Ensure window size is at least 3 (start, mid, end) to allow optimization
window_size = max(3, window_size);

n_wp = size(waypoints, 2);

% Initialize state at start
current_vel = [0; 0; 0];
current_acc = [0; 0; 0];
current_jerk = [0; 0; 0];

% Storage for polynomial pieces
all_coefs = [];
all_breaks = [0];
total_duration = 0;

fprintf('Sliding Window Trajectory: %d waypoints, window size %d\n', n_wp, window_size);

% Process in windows
start_idx = 1;

while start_idx < n_wp
    % Define window indices
    end_idx = min(start_idx + window_size - 1, n_wp);

    % If fewer than 2 points remaining (shouldn't happen if logic is correct)
    if end_idx <= start_idx
        break;
    end

    % Extract local waypoints
    local_wps = waypoints(:, start_idx:end_idx);
    n_local = size(local_wps, 2);

    % Time allocation
    % Time allocation with Adaptive Turn Slow-Down
    dists = sqrt(sum(diff(local_wps, 1, 2).^2, 1));

    % Calculate turning angles
    turn_factors = ones(1, n_local - 1);
    if n_local > 2
        for i = 2:n_local-1
            v1 = local_wps(:, i) - local_wps(:, i-1);
            v2 = local_wps(:, i+1) - local_wps(:, i);

            nv1 = norm(v1);
            nv2 = norm(v2);

            if nv1 > 1e-3 && nv2 > 1e-3
                cos_theta = dot(v1, v2) / (nv1 * nv2);
                % Clamp to [-1, 1]
                cos_theta = max(min(cos_theta, 1), -1);
                angle = acos(cos_theta);

                % Slow down factor: (1 + K * angle)
                % e.g., 90 deg turn (pi/2) -> factor ~ 2.5
                turn_factors(i-1) = 1 + 1.0 * angle;
                turn_factors(i) = max(turn_factors(i), 1 + 1.0 * angle);
            end
        end
    end

    dt_segments = (dists ./ V_avg) .* turn_factors;
    times = [0, cumsum(dt_segments)];

    % Boundary Conditions
    velBC = nan(3, n_local);
    accBC = nan(3, n_local);
    jerkBC = nan(3, n_local);

    % Start BC (Fixed)
    velBC(:, 1) = current_vel;
    accBC(:, 1) = current_acc;
    jerkBC(:, 1) = current_jerk;

    % End BC
    if end_idx == n_wp
        % Final goal: Stop
        velBC(:, end) = [0; 0; 0];
        accBC(:, end) = [0; 0; 0];
        jerkBC(:, end) = [0; 0; 0];
    else
        % Intermediate point: Free (NaN)
        % The optimizer will choose best derivatives for smoothness within the window
    end

    % Solve Minimum Snap for Window
    try
        numSamples = n_local * 10;
        [~,~,~,~,pp_local] = minsnappolytraj(local_wps, times, numSamples, ...
            'VelocityBoundaryCondition', velBC, ...
            'AccelerationBoundaryCondition', accBC, ...
            'JerkBoundaryCondition', jerkBC);

        % Extract trajectory to get end state
        t_end = times(end);

        % Evaluate local trajectory at the end
        [pos_end, vel_end, acc_end, jerk_end] = ppval_derivative(pp_local, t_end);

        % Update state for next window
        current_vel = vel_end;
        current_acc = acc_end;
        current_jerk = jerk_end;

        % Shift time allocation for global storage
        local_breaks = pp_local.breaks + total_duration;

        % Append coefficients (handling dimensionality)
        % pp.coefs is [dim * pieces, order]
        % We need to be careful with concatenation

        if isempty(all_coefs)
            all_coefs = pp_local.coefs;
            all_breaks = local_breaks;
        else
            all_coefs = [all_coefs; pp_local.coefs];
            all_breaks = [all_breaks, local_breaks(2:end)];
        end

        total_duration = total_duration + t_end;

    catch ME
        warning('Trajectory generation failed for window %d-%d: %s', start_idx, end_idx, ME.message);
        % Fallback: Linear interpolation?
    end

    % Move window
    % If we processed P1..P3, next window starts at P3
    start_idx = end_idx;
end

% Check if we generated any data
if isempty(all_coefs)
    error('Trajectory generation failed: No collision-free segments generated.');
end

% Construct full piecewise polynomial
full_pp = mkpp(all_breaks, all_coefs, 3); % 3D trajectory by default

% Sample the result
dt = 0.01;
t_eval = 0:dt:total_duration;

[pos, vel, acc, jerk, snap] = ppval_derivative(full_pp, t_eval);

full_traj.pp = full_pp;
full_traj.t = t_eval;
full_traj.x = pos;
full_traj.v = vel;
full_traj.a = acc;
full_traj.j = jerk;
full_traj.s = snap;

debug_info.duration = total_duration;
debug_info.n_poly = size(all_coefs, 1)/3;

end

function [pos, vel, acc, jerk, snap] = ppval_derivative(pp, t)
% Evaluate PP and its derivatives

% Position (0th derivative)
pos = ppval(pp, t);

if nargout > 1
    % Velocity (1st derivative)
    pp_dot = fnder(pp, 1);
    vel = ppval(pp_dot, t);
end

if nargout > 2
    % Acceleration (2nd derivative)
    pp_ddot = fnder(pp, 2);
    acc = ppval(pp_ddot, t);
end

if nargout > 3
    % Jerk
    pp_dddot = fnder(pp, 3);
    jerk = ppval(pp_dddot, t);
end

if nargout > 4
    % Snap
    pp_ddddot = fnder(pp, 4);
    snap = ppval(pp_ddddot, t);
end
end
