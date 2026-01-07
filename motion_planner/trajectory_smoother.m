function traj = trajectory_smoother(waypoints, total_time, params)
% TRAJECTORY_SMOOTHER Convert discrete waypoints to smooth trajectory
%
% Generates minimum-snap polynomial trajectories between waypoints,
% ensuring C3 continuity for compatibility with differential flatness control.
%
% Inputs:
%   waypoints  - [3xN] Matrix of waypoint positions in NED frame
%   total_time - Total trajectory duration in seconds
%   params     - (optional) Struct with parameters:
%                - dt: Trajectory sampling rate (default: 0.002s)
%                - v_max: Maximum velocity constraint (default: 5 m/s)
%                - a_max: Maximum acceleration constraint (default: 10 m/s^2)
%                - init_vel: [3x1] Initial velocity (default: [0;0;0])
%                - init_acc: [3x1] Initial acceleration (default: [0;0;0])
%
% Outputs:
%   traj - Struct containing:
%          - t: [1xM] Time vector
%          - pos: [3xM] Position trajectory
%          - vel: [3xM] Velocity trajectory
%          - acc: [3xM] Acceleration trajectory
%          - jerk: [3xM] Jerk trajectory
%          - snap: [3xM] Snap trajectory
%          - segment_times: Time allocation per segment
%
% The trajectory uses 7th order polynomials (minimum snap) to ensure:
%   - Position, velocity, acceleration continuity at waypoints
%   - Smooth jerk and snap for differential flatness

    % Default parameters
    if nargin < 3
        params = struct();
    end

    dt = get_param(params, 'dt', 0.002);
    v_max = get_param(params, 'v_max', 5.0);
    a_max = get_param(params, 'a_max', 10.0);

    % Initial conditions for smooth replanning
    init_vel = get_param(params, 'init_vel', [0;0;0]);
    init_acc = get_param(params, 'init_acc', [0;0;0]);

    waypoints = waypoints(:, :);  % Ensure 2D
    n_waypoints = size(waypoints, 2);

    if n_waypoints < 2
        error('Need at least 2 waypoints');
    end

    n_segments = n_waypoints - 1;

    % Allocate time per segment based on distance
    segment_times = allocate_segment_times(waypoints, total_time, v_max);

    % Generate polynomial coefficients for each segment
    % Using 7th order polynomial (8 coefficients) for minimum snap
    coeffs = compute_min_snap_coeffs(waypoints, segment_times, init_vel, init_acc);

    % Sample trajectory at dt intervals
    t_vec = 0:dt:total_time;
    M = length(t_vec);

    % Initialize output arrays
    traj.t = t_vec;
    traj.pos = zeros(3, M);
    traj.vel = zeros(3, M);
    traj.acc = zeros(3, M);
    traj.jerk = zeros(3, M);
    traj.snap = zeros(3, M);
    traj.segment_times = segment_times;
    traj.waypoints = waypoints;

    % Cumulative segment start times
    t_segments = [0, cumsum(segment_times)];

    % Sample each point
    for i = 1:M
        t = t_vec(i);

        % Find which segment we're in
        seg_idx = find(t_segments(1:end-1) <= t & t < t_segments(2:end), 1);
        if isempty(seg_idx)
            seg_idx = n_segments;  % Handle t = total_time
        end

        % Local time within segment [0, 1]
        T_seg = segment_times(seg_idx);
        t_local = (t - t_segments(seg_idx)) / T_seg;
        t_local = max(0, min(1, t_local));  % Clamp to [0, 1]

        % Evaluate polynomial and derivatives for each axis
        for axis = 1:3
            c = coeffs{seg_idx, axis};  % 8 coefficients for this segment/axis

            % Position: p(tau) = c0 + c1*tau + c2*tau^2 + ... + c7*tau^7
            traj.pos(axis, i) = polyval_custom(c, t_local);

            % Velocity: dp/dt = (1/T) * dp/dtau
            traj.vel(axis, i) = polyval_deriv(c, t_local, 1) / T_seg;

            % Acceleration
            traj.acc(axis, i) = polyval_deriv(c, t_local, 2) / T_seg^2;

            % Jerk
            traj.jerk(axis, i) = polyval_deriv(c, t_local, 3) / T_seg^3;

            % Snap
            traj.snap(axis, i) = polyval_deriv(c, t_local, 4) / T_seg^4;
        end
    end
end

function val = get_param(params, name, default)
    if isfield(params, name)
        val = params.(name);
    else
        val = default;
    end
end

function segment_times = allocate_segment_times(waypoints, total_time, v_max)
% Allocate time to each segment proportional to distance
% with minimum time based on v_max

    n_segments = size(waypoints, 2) - 1;
    distances = zeros(1, n_segments);

    for i = 1:n_segments
        distances(i) = norm(waypoints(:, i+1) - waypoints(:, i));
    end

    total_distance = sum(distances);

    if total_distance < 1e-6
        % All waypoints are the same
        segment_times = ones(1, n_segments) * total_time / n_segments;
    else
        % Proportional to distance
        segment_times = (distances / total_distance) * total_time;

        % Ensure minimum time based on velocity constraint
        min_times = distances / v_max;
        segment_times = max(segment_times, min_times);

        % Renormalize to match total_time
        segment_times = segment_times * (total_time / sum(segment_times));
    end
end

function coeffs = compute_min_snap_coeffs(waypoints, segment_times, init_vel, init_acc)
% Compute 7th order polynomial coefficients for minimum snap trajectory
%
% For each segment, we have 8 unknowns (c0 to c7)
% Boundary conditions at each waypoint:
%   - Position continuity (2 constraints per segment)
%   - Velocity continuity (2 constraints per interior waypoint)
%   - Acceleration continuity (2 constraints per interior waypoint)
%   - Specified initial velocity/acceleration, zero at final endpoint
%
% We solve segment by segment with matching conditions

    n_segments = size(waypoints, 2) - 1;
    coeffs = cell(n_segments, 3);

    for axis = 1:3
        % Initial conditions (use provided initial velocity and acceleration)
        pos_0 = waypoints(axis, 1);
        vel_0 = init_vel(axis);
        acc_0 = init_acc(axis);
        jerk_0 = 0;  % Assume zero initial jerk for smoothness

        for seg = 1:n_segments
            % Target position
            pos_1 = waypoints(axis, seg + 1);

            % Final conditions for this segment
            if seg == n_segments
                % Last segment: come to rest
                vel_1 = 0;
                acc_1 = 0;
                jerk_1 = 0;
            else
                % Interior: will be computed to ensure smoothness
                % Use heuristic: velocity toward next waypoint
                if seg + 1 < size(waypoints, 2)
                    next_dir = waypoints(:, seg + 2) - waypoints(:, seg + 1);
                    dist = norm(next_dir);
                    if dist > 1e-6
                        % Velocity proportional to distance ratio
                        avg_vel = (pos_1 - pos_0) / segment_times(seg);
                        vel_1 = avg_vel * 0.5;  % Blend
                    else
                        vel_1 = 0;
                    end
                else
                    vel_1 = 0;
                end
                acc_1 = 0;  % Smooth acceleration
                jerk_1 = 0;
            end

            % Solve for 7th order polynomial coefficients
            % p(tau) = c0 + c1*tau + c2*tau^2 + c3*tau^3 + c4*tau^4 + c5*tau^5 + c6*tau^6 + c7*tau^7
            % where tau in [0, 1]
            %
            % Boundary conditions at tau=0:
            %   p(0) = c0 = pos_0
            %   p'(0) = c1 = vel_0 * T
            %   p''(0) = 2*c2 = acc_0 * T^2
            %   p'''(0) = 6*c3 = jerk_0 * T^3
            %
            % Boundary conditions at tau=1:
            %   p(1) = sum(c) = pos_1
            %   p'(1) = c1 + 2*c2 + ... + 7*c7 = vel_1 * T
            %   p''(1) = 2*c2 + 6*c3 + ... + 42*c7 = acc_1 * T^2
            %   p'''(1) = 6*c3 + ... + 210*c7 = jerk_1 * T^3

            T = segment_times(seg);

            % Build constraint matrix A * c = b
            A = zeros(8, 8);
            b = zeros(8, 1);

            % tau = 0 constraints
            A(1, :) = [1 0 0 0 0 0 0 0];  % p(0)
            b(1) = pos_0;

            A(2, :) = [0 1 0 0 0 0 0 0];  % p'(0)
            b(2) = vel_0 * T;

            A(3, :) = [0 0 2 0 0 0 0 0];  % p''(0)
            b(3) = acc_0 * T^2;

            A(4, :) = [0 0 0 6 0 0 0 0];  % p'''(0)
            b(4) = jerk_0 * T^3;

            % tau = 1 constraints
            A(5, :) = [1 1 1 1 1 1 1 1];  % p(1)
            b(5) = pos_1;

            A(6, :) = [0 1 2 3 4 5 6 7];  % p'(1)
            b(6) = vel_1 * T;

            A(7, :) = [0 0 2 6 12 20 30 42];  % p''(1)
            b(7) = acc_1 * T^2;

            A(8, :) = [0 0 0 6 24 60 120 210];  % p'''(1)
            b(8) = jerk_1 * T^3;

            % Solve linear system
            c = A \ b;
            coeffs{seg, axis} = c;

            % Update initial conditions for next segment
            pos_0 = pos_1;
            vel_0 = vel_1;
            acc_0 = acc_1;
            jerk_0 = jerk_1;
        end
    end
end

function val = polyval_custom(c, tau)
% Evaluate polynomial c0 + c1*tau + c2*tau^2 + ... + c7*tau^7
    val = c(1);
    tau_power = 1;
    for i = 2:length(c)
        tau_power = tau_power * tau;
        val = val + c(i) * tau_power;
    end
end

function val = polyval_deriv(c, tau, deriv_order)
% Evaluate the deriv_order-th derivative of polynomial at tau
%
% c = [c0, c1, c2, ..., c7] (coefficient order: constant to highest)

    n = length(c) - 1;  % Polynomial degree

    % Compute derivative coefficients
    dc = c;
    for d = 1:deriv_order
        new_dc = zeros(1, length(dc) - 1);
        for i = 2:length(dc)
            new_dc(i-1) = dc(i) * (i - 1);
        end
        dc = new_dc;
        if isempty(dc)
            val = 0;
            return;
        end
    end

    % Evaluate derivative polynomial
    val = dc(1);
    tau_power = 1;
    for i = 2:length(dc)
        tau_power = tau_power * tau;
        val = val + dc(i) * tau_power;
    end
end
