function traj = trajectory_smoother(waypoints, total_time, params)
% TRAJECTORY_SMOOTHER Convert discrete waypoints to smooth trajectory
%
% Uses MATLAB's minsnappolytraj function to generate minimum-snap polynomial
% trajectories between waypoints with free interior derivatives.
%
% Inputs:
%   waypoints  - [3xN] Matrix of waypoint positions in NED frame
%   total_time - Total trajectory duration in seconds
%   params     - (optional) Struct with parameters:
%                - dt: Trajectory sampling rate (default: 0.002s)
%                - v_max: Maximum velocity constraint (default: 5 m/s)
%                - init_vel: [3x1] Initial velocity (default: [0;0;0])
%                - init_acc: [3x1] Initial acceleration (default: [0;0;0])
%                - final_vel: [3x1] Final velocity (default: [0;0;0])
%                - final_acc: [3x1] Final acceleration (default: [0;0;0])
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
%          - waypoints: Input waypoints

    % Default parameters
    if nargin < 3
        params = struct();
    end

    dt = get_param(params, 'dt', 0.002);
    v_max = get_param(params, 'v_max', 5.0);

    % Boundary conditions for smooth replanning
    init_vel = get_param(params, 'init_vel', [0;0;0]);
    init_acc = get_param(params, 'init_acc', [0;0;0]);
    final_vel = get_param(params, 'final_vel', [0;0;0]);
    final_acc = get_param(params, 'final_acc', [0;0;0]);

    waypoints = waypoints(:, :);  % Ensure 2D
    n_waypoints = size(waypoints, 2);

    if n_waypoints < 2
        error('Need at least 2 waypoints');
    end

    % Allocate time per segment based on distance
    segment_times = allocate_segment_times(waypoints, total_time, v_max);

    % Create time points for each waypoint
    timePoints = [0, cumsum(segment_times)];

    % Number of samples
    numSamples = round(total_time / dt) + 1;

    % Setup velocity boundary conditions
    % NaN means "free" - let optimizer choose optimal value
    velBC = nan(3, n_waypoints);
    velBC(:, 1) = init_vel;          % Initial velocity: specified
    velBC(:, end) = final_vel;       % Final velocity: specified
    % Interior waypoints: NaN (free)

    % Setup acceleration boundary conditions
    accBC = nan(3, n_waypoints);
    accBC(:, 1) = init_acc;          % Initial acceleration: specified
    accBC(:, end) = final_acc;       % Final acceleration: specified
    % Interior waypoints: NaN (free)

    % Call MATLAB's minsnappolytraj
    [pos, vel, acc, jerk, snap, ~, ~, tSamples] = minsnappolytraj(...
        waypoints, ...
        timePoints, ...
        numSamples, ...
        'VelocityBoundaryCondition', velBC, ...
        'AccelerationBoundaryCondition', accBC);

    % Build output struct
    traj.t = tSamples;
    traj.pos = pos;
    traj.vel = vel;
    traj.acc = acc;
    traj.jerk = jerk;
    traj.snap = snap;
    traj.segment_times = segment_times;
    traj.waypoints = waypoints;
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
