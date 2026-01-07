classdef waypoint_manager < handle
% WAYPOINT_MANAGER Manages waypoint sequencing and mission execution
%
% This class handles:
%   - Mission definition and waypoint sequencing
%   - Real-time trajectory generation from waypoints
%   - Dynamic replanning when obstacles are detected
%   - Progress tracking toward goal
%
% Usage:
%   wm = waypoint_manager();
%   wm.set_mission(start_pos, goal_pos, obstacles, bounds);
%   wm.plan_path();
%   wm.generate_trajectory(total_time);
%
%   % During simulation:
%   ref = wm.get_reference(t);
%   if wm.needs_replan(current_pos, new_obstacles)
%       wm.replan(current_pos, goal, new_obstacles);
%   end

    properties
        % Mission parameters
        start_pos       % [3x1] Mission start position
        goal_pos        % [3x1] Mission goal position
        obstacles       % Cell array of obstacles
        bounds          % Workspace bounds struct

        % Path and trajectory
        waypoints       % [3xN] Planned waypoints
        trajectory      % Struct from trajectory_smoother

        % Planning parameters
        rrt_params      % RRT planner parameters
        traj_params     % Trajectory smoother parameters

        % State
        is_planned      % true if path has been planned
        is_trajectory_generated  % true if smooth trajectory exists
        total_time      % Total mission time
        current_time    % Current time in mission
    end

    methods
        function obj = waypoint_manager()
            % Constructor - initialize empty mission
            obj.start_pos = zeros(3, 1);
            obj.goal_pos = zeros(3, 1);
            obj.obstacles = {};
            obj.bounds = struct('min', [-50; -50; -50], 'max', [50; 50; 0]);

            obj.waypoints = [];
            obj.trajectory = struct();

            % Default RRT parameters
            obj.rrt_params = struct(...
                'max_iter', 5000, ...
                'step_size', 0.5, ...
                'goal_bias', 0.15, ...
                'goal_tol', 0.5, ...
                'safety_margin', 0.5);

            % Default trajectory parameters
            obj.traj_params = struct(...
                'dt', 0.002, ...
                'v_max', 3.0, ...
                'a_max', 8.0);

            obj.is_planned = false;
            obj.is_trajectory_generated = false;
            obj.total_time = 10.0;
            obj.current_time = 0;
        end

        function set_mission(obj, start_pos, goal_pos, obstacles, bounds)
            % SET_MISSION Define a new mission
            %
            % Inputs:
            %   start_pos - [3x1] Starting position
            %   goal_pos  - [3x1] Goal position
            %   obstacles - Cell array of obstacle structs (optional)
            %   bounds    - Workspace bounds struct (optional)

            obj.start_pos = start_pos(:);
            obj.goal_pos = goal_pos(:);

            if nargin >= 4 && ~isempty(obstacles)
                obj.obstacles = obstacles;
            else
                obj.obstacles = {};
            end

            if nargin >= 5 && ~isempty(bounds)
                obj.bounds = bounds;
            end

            obj.is_planned = false;
            obj.is_trajectory_generated = false;
            obj.current_time = 0;
        end

        function [success, path] = plan_path(obj)
            % PLAN_PATH Generate collision-free path using RRT
            %
            % Outputs:
            %   success - true if path found
            %   path    - [3xN] Waypoint path

            if isempty(obj.obstacles)
                % No obstacles: direct path with intermediate point
                obj.waypoints = [obj.start_pos, obj.goal_pos];
                success = true;
                fprintf('Waypoint Manager: Direct path (no obstacles)\n');
            else
                % Use RRT planner
                [path, success, ~] = rrt_planner(...
                    obj.start_pos, obj.goal_pos, ...
                    obj.obstacles, obj.bounds, obj.rrt_params);

                if success
                    % Simplify path by removing collinear points
                    obj.waypoints = simplify_path(path, obj.obstacles, obj.rrt_params.safety_margin);
                    fprintf('Waypoint Manager: RRT path with %d waypoints\n', size(obj.waypoints, 2));
                else
                    obj.waypoints = obj.start_pos;
                    warning('Waypoint Manager: Path planning failed!');
                end
            end

            obj.is_planned = success;
            path = obj.waypoints;
        end

        function traj = generate_trajectory(obj, total_time)
            % GENERATE_TRAJECTORY Create smooth trajectory from waypoints
            %
            % Inputs:
            %   total_time - Total mission duration in seconds
            %
            % Outputs:
            %   traj - Trajectory struct with pos, vel, acc, jerk, snap

            if ~obj.is_planned
                error('Must plan path before generating trajectory');
            end

            obj.total_time = total_time;

            % Generate smooth trajectory
            obj.trajectory = trajectory_smoother(obj.waypoints, total_time, obj.traj_params);

            obj.is_trajectory_generated = true;
            traj = obj.trajectory;
            fprintf('Waypoint Manager: Trajectory generated (%.1fs, %d samples)\n', ...
                total_time, length(obj.trajectory.t));
        end

        function ref = get_reference(obj, t)
            % GET_REFERENCE Get reference state at time t
            %
            % Inputs:
            %   t - Current time
            %
            % Outputs:
            %   ref - Struct with pos, vel, acc, jerk, snap, psi, psi_dot, psi_ddot

            if ~obj.is_trajectory_generated
                error('Trajectory not generated. Call generate_trajectory first.');
            end

            obj.current_time = t;

            % Find trajectory index for this time
            t_vec = obj.trajectory.t;
            idx = find(t_vec <= t, 1, 'last');
            if isempty(idx)
                idx = 1;
            end
            if idx > length(t_vec)
                idx = length(t_vec);
            end

            % Extract reference
            ref.pos = obj.trajectory.pos(:, idx);
            ref.vel = obj.trajectory.vel(:, idx);
            ref.acc = obj.trajectory.acc(:, idx);
            ref.jerk = obj.trajectory.jerk(:, idx);
            ref.snap = obj.trajectory.snap(:, idx);

            % Yaw: point toward velocity direction (tangent mode)
            v_horiz = ref.vel(1:2);
            if norm(v_horiz) > 0.1
                ref.psi = atan2(v_horiz(2), v_horiz(1));
            else
                ref.psi = 0;
            end
            ref.psi_dot = 0;
            ref.psi_ddot = 0;
        end

        function [needs_replan, reason] = check_replan_needed(obj, current_pos, obstacles)
            % CHECK_REPLAN_NEEDED Check if replanning is required
            %
            % Inputs:
            %   current_pos - Current drone position
            %   obstacles   - Updated obstacle list
            %
            % Outputs:
            %   needs_replan - true if replanning needed
            %   reason       - String describing why replan is needed

            needs_replan = false;
            reason = '';

            if ~obj.is_trajectory_generated
                return;
            end

            % Check if remaining path collides with obstacles
            remaining_path = obj.get_remaining_path(current_pos);

            for i = 1:size(remaining_path, 2) - 1
                p1 = remaining_path(:, i);
                p2 = remaining_path(:, i + 1);

                % Sample along segment
                dist = norm(p2 - p1);
                n_samples = max(2, ceil(dist / 0.2));

                for j = 0:n_samples
                    alpha = j / n_samples;
                    p_sample = p1 + alpha * (p2 - p1);

                    [is_collision, ~, ~] = collision_checker(p_sample, obstacles, obj.rrt_params.safety_margin);

                    if is_collision
                        needs_replan = true;
                        reason = 'Path collision with obstacle detected';
                        return;
                    end
                end
            end
        end

        function [success, new_traj] = replan(obj, current_pos, current_time, obstacles, current_vel, current_acc)
            % REPLAN Replan trajectory from current position with smooth transition
            %
            % Inputs:
            %   current_pos  - Current drone position [3x1]
            %   current_time - Current simulation time
            %   obstacles    - Updated obstacle list
            %   current_vel  - (optional) Current velocity [3x1] for smooth transition
            %   current_acc  - (optional) Current acceleration [3x1] for smooth transition
            %
            % Outputs:
            %   success  - true if replanning succeeded
            %   new_traj - New trajectory struct

            fprintf('[%.2fs] Replanning from [%.1f, %.1f, %.1f]...\n', ...
                current_time, current_pos);

            % Handle optional velocity/acceleration arguments
            if nargin < 5 || isempty(current_vel)
                current_vel = [0; 0; 0];
            end
            if nargin < 6 || isempty(current_acc)
                current_acc = [0; 0; 0];
            end

            % Update mission with current position as start
            obj.start_pos = current_pos(:);
            obj.obstacles = obstacles;

            % Replan path
            [success, ~] = obj.plan_path();

            if success
                % Calculate remaining time (add extra time for smoother replanning)
                remaining_time = max(obj.total_time - current_time, 4.0);

                % Store initial conditions in traj_params for smooth transition
                obj.traj_params.init_vel = current_vel(:);
                obj.traj_params.init_acc = current_acc(:);

                % Generate new trajectory with smooth initial conditions
                new_traj = obj.generate_trajectory(remaining_time);

                % Reset init conditions for future normal trajectories
                obj.traj_params.init_vel = [0; 0; 0];
                obj.traj_params.init_acc = [0; 0; 0];

                % Adjust time base for continuation
                obj.trajectory.t = obj.trajectory.t + current_time;
                new_traj = obj.trajectory;

                fprintf('[%.2fs] Replanning successful (smooth transition with v=[%.2f,%.2f,%.2f])\n', ...
                    current_time, current_vel);
            else
                new_traj = obj.trajectory;
                warning('Replanning failed!');
            end
        end

        function remaining = get_remaining_path(obj, current_pos)
            % GET_REMAINING_PATH Get waypoints from current position to goal
            %
            % Finds nearest waypoint and returns remaining path

            if isempty(obj.waypoints)
                remaining = current_pos;
                return;
            end

            % Find nearest waypoint
            min_dist = inf;
            nearest_idx = 1;
            for i = 1:size(obj.waypoints, 2)
                d = norm(obj.waypoints(:, i) - current_pos);
                if d < min_dist
                    min_dist = d;
                    nearest_idx = i;
                end
            end

            % Return path from current position through remaining waypoints
            remaining = [current_pos, obj.waypoints(:, nearest_idx:end)];
        end

        function progress = get_progress(obj)
            % GET_PROGRESS Return mission progress (0 to 1)

            if obj.total_time <= 0
                progress = 0;
            else
                progress = min(1, obj.current_time / obj.total_time);
            end
        end

        function dist = distance_to_goal(obj, current_pos)
            % DISTANCE_TO_GOAL Return distance from current position to goal
            dist = norm(current_pos(:) - obj.goal_pos);
        end

        function reached = goal_reached(obj, current_pos, tolerance)
            % GOAL_REACHED Check if goal has been reached
            if nargin < 3
                tolerance = 0.5;
            end
            reached = obj.distance_to_goal(current_pos) < tolerance;
        end
    end
end

function simplified = simplify_path(path, obstacles, safety_margin)
% SIMPLIFY_PATH Remove unnecessary waypoints while maintaining collision-free path

    n_points = size(path, 2);
    if n_points <= 2
        simplified = path;
        return;
    end

    % Keep first point
    simplified = path(:, 1);
    current_idx = 1;

    while current_idx < n_points
        % Try to skip to furthest visible point
        furthest_visible = current_idx + 1;

        for check_idx = n_points:-1:(current_idx + 2)
            if is_segment_clear(path(:, current_idx), path(:, check_idx), obstacles, safety_margin)
                furthest_visible = check_idx;
                break;
            end
        end

        % Add furthest visible point
        simplified = [simplified, path(:, furthest_visible)]; %#ok<AGROW>
        current_idx = furthest_visible;
    end
end

function clear = is_segment_clear(p1, p2, obstacles, safety_margin)
% Check if line segment is collision-free

    if isempty(obstacles)
        clear = true;
        return;
    end

    dist = norm(p2 - p1);
    n_samples = max(2, ceil(dist / 0.1));

    for i = 0:n_samples
        alpha = i / n_samples;
        p_sample = p1 + alpha * (p2 - p1);

        [is_collision, ~, ~] = collision_checker(p_sample, obstacles, safety_margin);

        if is_collision
            clear = false;
            return;
        end
    end

    clear = true;
end
