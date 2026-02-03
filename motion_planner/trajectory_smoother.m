function traj = trajectory_smoother(waypoints, total_time, params)
% TRAJECTORY_SMOOTHER Convert discrete waypoints to smooth trajectory
%
% Uses MATLAB's minsnappolytraj function to generate minimum-snap polynomial
% trajectories between waypoints.
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
%                - vel_bc_mode: 'free' or 'bisector' (default: 'free')
%                   'free'    - Interior velocities are NaN (optimizer chooses)
%                   'bisector'- Interior velocities set to cruise_speed * bisector direction
%                - cruise_speed: Speed magnitude for bisector mode (default: v_max * 0.7)
%                - max_waypoints: Maximum waypoints before additional simplification (default: 25)
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
%          - velBC: [3xN] Velocity boundary conditions used (for visualization)

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

    % Velocity BC mode: 'free' (NaN) or 'bisector' (direction-based)
    vel_bc_mode = get_param(params, 'vel_bc_mode', 'free');
    cruise_speed = get_param(params, 'cruise_speed', v_max * 0.7);

    % Maximum waypoints before minsnappolytraj becomes poorly conditioned
    max_waypoints = get_param(params, 'max_waypoints', 25);

    waypoints = waypoints(:, :);  % Ensure 2D
    original_count = size(waypoints, 2);

    % --- STEP 1: FILTER DUPLICATE WAYPOINTS ---
    % Remove consecutive points that are too close (causes non-increasing time)
    min_dist = 10;  % Minimum distance between waypoints (meters)
    keep_mask = true(1, size(waypoints, 2));
    for i = 2:size(waypoints, 2)
        if norm(waypoints(:, i) - waypoints(:, i-1)) < min_dist
            keep_mask(i) = false;
        end
    end
    waypoints = waypoints(:, keep_mask);

    % --- STEP 2: ADDITIONAL SIMPLIFICATION IF TOO MANY WAYPOINTS ---
    % minsnappolytraj becomes poorly conditioned with >25-30 waypoints
    if size(waypoints, 2) > max_waypoints
        % Compute path length for adaptive tolerance
        path_len = 0;
        for i = 2:size(waypoints, 2)
            path_len = path_len + norm(waypoints(:, i) - waypoints(:, i-1));
        end

        % Iteratively increase tolerance until waypoints <= max_waypoints
        base_tol = path_len * 0.01;  % Start with 1% of path length
        tol = base_tol;
        max_tol = path_len * 0.1;    % Max 10% of path length

        while size(waypoints, 2) > max_waypoints && tol < max_tol
            waypoints_simplified = simplify_path_internal(waypoints, tol);
            if size(waypoints_simplified, 2) <= max_waypoints
                waypoints = waypoints_simplified;
                break;
            end
            tol = tol * 1.5;  % Increase tolerance by 50%
        end

        % If still too many, force uniform sampling
        if size(waypoints, 2) > max_waypoints
            indices = round(linspace(1, size(waypoints, 2), max_waypoints));
            waypoints = waypoints(:, indices);
            warning('Trajectory smoother: Forced uniform sampling to %d waypoints', max_waypoints);
        end

        fprintf('Trajectory smoother: %d -> %d waypoints (conditioning limit: %d)\n', ...
            original_count, size(waypoints, 2), max_waypoints);
    else
        fprintf('Trajectory smoother: %d of %d waypoints kept after filtering\n', ...
            size(waypoints, 2), original_count);
    end

    n_waypoints = size(waypoints, 2);

    if n_waypoints < 2
        error('Need at least 2 waypoints');
    end

    % Warn if still potentially problematic
    if n_waypoints > 20
        warning('Trajectory smoother: %d waypoints may cause numerical issues', n_waypoints);
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

    % Compute interior waypoint velocities based on mode
    if strcmpi(vel_bc_mode, 'bisector') && n_waypoints > 2
        % Bisector mode: velocity direction is average of incoming/outgoing directions
        for i = 2:n_waypoints-1
            % Incoming direction (from previous waypoint)
            dir_in = waypoints(:, i) - waypoints(:, i-1);
            len_in = norm(dir_in);
            if len_in > 1e-6
                dir_in = dir_in / len_in;
            else
                dir_in = [0; 0; 0];
            end

            % Outgoing direction (to next waypoint)
            dir_out = waypoints(:, i+1) - waypoints(:, i);
            len_out = norm(dir_out);
            if len_out > 1e-6
                dir_out = dir_out / len_out;
            else
                dir_out = [0; 0; 0];
            end

            % Bisector direction (average of incoming and outgoing)
            dir_bisect = dir_in + dir_out;
            len_bisect = norm(dir_bisect);
            if len_bisect > 1e-6
                dir_bisect = dir_bisect / len_bisect;
            else
                % Near 180-degree turn - use incoming direction
                dir_bisect = dir_in;
            end

            % Set velocity BC at this waypoint
            velBC(:, i) = cruise_speed * dir_bisect;
        end
        fprintf('Velocity BC mode: BISECTOR (cruise speed: %.1f m/s)\n', cruise_speed);
    else
        % Free mode: interior velocities are NaN (optimizer chooses)
        fprintf('Velocity BC mode: FREE (optimizer chooses interior velocities)\n');
    end

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
    traj.velBC = velBC;  % Velocity boundary conditions (for visualization)
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

function simplified = simplify_path_internal(path, tolerance)
% SIMPLIFY_PATH_INTERNAL Douglas-Peucker path simplification (internal version)
% Reduces waypoints while preserving path shape within tolerance

    n = size(path, 2);
    if n <= 2
        simplified = path;
        return;
    end

    % Find point with maximum distance from line connecting endpoints
    p1 = path(:, 1);
    p2 = path(:, end);
    line_vec = p2 - p1;
    line_len = norm(line_vec);

    if line_len < 1e-6
        simplified = path(:, [1, end]);
        return;
    end

    line_unit = line_vec / line_len;

    max_dist = 0;
    max_idx = 1;

    for i = 2:n-1
        v = path(:, i) - p1;
        proj_len = dot(v, line_unit);
        proj_len = max(0, min(line_len, proj_len));
        proj_point = p1 + proj_len * line_unit;
        dist = norm(path(:, i) - proj_point);

        if dist > max_dist
            max_dist = dist;
            max_idx = i;
        end
    end

    if max_dist > tolerance
        left = simplify_path_internal(path(:, 1:max_idx), tolerance);
        right = simplify_path_internal(path(:, max_idx:end), tolerance);
        simplified = [left, right(:, 2:end)];
    else
        simplified = path(:, [1, end]);
    end
end
