classdef obstacle_manager < handle
% OBSTACLE_MANAGER Manages static and dynamic obstacles for motion planning
%
% This class maintains a list of obstacles that can be:
%   - Static (defined at initialization)
%   - Dynamic "pop-up" threats (appear at specified simulation time)
%
% Usage:
%   om = obstacle_manager();
%   om.add_obstacle('sphere', [5;2;-2], struct('radius', 1.0));
%   om.add_popup_threat('sphere', [8;3;-2], struct('radius', 1.5), 5.0);
%
%   % During simulation loop:
%   om.update(current_time);
%   obstacles = om.get_active_obstacles();
%
% Properties:
%   obstacles       - Cell array of all obstacle structs
%   popup_threats   - Struct array of threats with activation times
%   active_indices  - Indices of currently active obstacles

    properties
        obstacles       % Cell array of obstacle structs
        popup_threats   % Struct array: {obs_idx, activation_time, activated}
        current_time    % Current simulation time
    end

    methods
        function obj = obstacle_manager()
            % Constructor - initialize empty obstacle list
            obj.obstacles = {};
            obj.popup_threats = struct('obs_idx', {}, 'activation_time', {}, 'activated', {});
            obj.current_time = 0;
        end

        function idx = add_obstacle(obj, type, center, params)
            % ADD_OBSTACLE Add a static obstacle to the environment
            %
            % Inputs:
            %   type   - 'sphere' or 'box'
            %   center - [3x1] Position of obstacle center in NED
            %   params - Struct with type-specific parameters:
            %            For 'sphere': params.radius
            %            For 'box': params.half_size [3x1]
            %
            % Returns:
            %   idx - Index of the added obstacle

            obs = struct();
            obs.type = type;
            obs.center = center(:);  % Ensure column vector
            obs.active = true;       % Static obstacles are always active

            switch type
                case 'sphere'
                    if ~isfield(params, 'radius')
                        error('Sphere obstacle requires params.radius');
                    end
                    obs.radius = params.radius;

                case 'box'
                    if ~isfield(params, 'half_size')
                        error('Box obstacle requires params.half_size');
                    end
                    obs.half_size = params.half_size(:);

                otherwise
                    error('Unknown obstacle type: %s', type);
            end

            obj.obstacles{end+1} = obs;
            idx = length(obj.obstacles);
        end

        function idx = add_popup_threat(obj, type, center, params, activation_time)
            % ADD_POPUP_THREAT Add a dynamic threat that appears at specified time
            %
            % Inputs:
            %   type            - 'sphere' or 'box'
            %   center          - [3x1] Position where threat will appear
            %   params          - Type-specific parameters (same as add_obstacle)
            %   activation_time - Simulation time (seconds) when threat appears
            %
            % Returns:
            %   idx - Index of the added obstacle

            % Add as inactive obstacle first
            obs = struct();
            obs.type = type;
            obs.center = center(:);
            obs.active = false;  % Initially inactive

            switch type
                case 'sphere'
                    obs.radius = params.radius;
                case 'box'
                    obs.half_size = params.half_size(:);
            end

            obj.obstacles{end+1} = obs;
            idx = length(obj.obstacles);

            % Record as popup threat
            threat = struct();
            threat.obs_idx = idx;
            threat.activation_time = activation_time;
            threat.activated = false;

            obj.popup_threats(end+1) = threat;
        end

        function [new_threats_activated, threat_indices] = update(obj, t)
            % UPDATE Update obstacle states based on simulation time
            %
            % Inputs:
            %   t - Current simulation time
            %
            % Outputs:
            %   new_threats_activated - true if any new threats just activated
            %   threat_indices        - Indices of newly activated threats

            obj.current_time = t;
            new_threats_activated = false;
            threat_indices = [];

            % Check each popup threat
            for i = 1:length(obj.popup_threats)
                if ~obj.popup_threats(i).activated && ...
                   t >= obj.popup_threats(i).activation_time
                    % Activate this threat
                    obs_idx = obj.popup_threats(i).obs_idx;
                    obj.obstacles{obs_idx}.active = true;
                    obj.popup_threats(i).activated = true;

                    new_threats_activated = true;
                    threat_indices(end+1) = obs_idx; %#ok<AGROW>

                    fprintf('[%.2fs] Pop-up threat activated at [%.1f, %.1f, %.1f]\n', ...
                        t, obj.obstacles{obs_idx}.center);
                end
            end
        end

        function active_obs = get_active_obstacles(obj)
            % GET_ACTIVE_OBSTACLES Return cell array of currently active obstacles
            %
            % Returns:
            %   active_obs - Cell array of active obstacle structs for collision checking

            active_obs = {};
            for i = 1:length(obj.obstacles)
                if obj.obstacles{i}.active
                    active_obs{end+1} = obj.obstacles{i}; %#ok<AGROW>
                end
            end
        end

        function all_obs = get_all_obstacles(obj)
            % GET_ALL_OBSTACLES Return all obstacles (for visualization)
            all_obs = obj.obstacles;
        end

        function remove_obstacle(obj, idx)
            % REMOVE_OBSTACLE Remove an obstacle by index
            if idx > 0 && idx <= length(obj.obstacles)
                obj.obstacles(idx) = [];
                % Update popup threat indices
                for i = 1:length(obj.popup_threats)
                    if obj.popup_threats(i).obs_idx > idx
                        obj.popup_threats(i).obs_idx = obj.popup_threats(i).obs_idx - 1;
                    elseif obj.popup_threats(i).obs_idx == idx
                        obj.popup_threats(i) = [];
                    end
                end
            end
        end

        function clear_all(obj)
            % CLEAR_ALL Remove all obstacles
            obj.obstacles = {};
            obj.popup_threats = struct('obs_idx', {}, 'activation_time', {}, 'activated', {});
        end

        function n = get_obstacle_count(obj)
            % GET_OBSTACLE_COUNT Return total number of obstacles
            n = length(obj.obstacles);
        end

        function n = get_active_count(obj)
            % GET_ACTIVE_COUNT Return number of currently active obstacles
            n = 0;
            for i = 1:length(obj.obstacles)
                if obj.obstacles{i}.active
                    n = n + 1;
                end
            end
        end

        function visualize(obj, ax)
            % VISUALIZE Draw all obstacles on given axes
            %
            % Inputs:
            %   ax - Axes handle (optional, uses gca if not provided)

            if nargin < 2
                ax = gca;
            end

            hold(ax, 'on');

            for i = 1:length(obj.obstacles)
                obs = obj.obstacles{i};

                if obs.active
                    color = [1, 0.3, 0.3];  % Red for active
                    alpha = 0.6;
                else
                    color = [0.5, 0.5, 0.5];  % Gray for inactive
                    alpha = 0.2;
                end

                switch obs.type
                    case 'sphere'
                        draw_sphere(ax, obs.center, obs.radius, color, alpha);
                    case 'box'
                        draw_box(ax, obs.center, obs.half_size, color, alpha);
                end
            end
        end
    end
end

function draw_sphere(ax, center, radius, color, alpha)
% Draw a sphere obstacle
    [X, Y, Z] = sphere(20);
    X = radius * X + center(1);
    Y = radius * Y + center(2);
    Z = radius * Z + center(3);
    surf(ax, X, Y, Z, 'FaceColor', color, 'FaceAlpha', alpha, ...
         'EdgeColor', 'none', 'DisplayName', 'Obstacle');
end

function draw_box(ax, center, half_size, color, alpha)
% Draw a box obstacle
    vertices = [
        center(1) - half_size(1), center(2) - half_size(2), center(3) - half_size(3);
        center(1) + half_size(1), center(2) - half_size(2), center(3) - half_size(3);
        center(1) + half_size(1), center(2) + half_size(2), center(3) - half_size(3);
        center(1) - half_size(1), center(2) + half_size(2), center(3) - half_size(3);
        center(1) - half_size(1), center(2) - half_size(2), center(3) + half_size(3);
        center(1) + half_size(1), center(2) - half_size(2), center(3) + half_size(3);
        center(1) + half_size(1), center(2) + half_size(2), center(3) + half_size(3);
        center(1) - half_size(1), center(2) + half_size(2), center(3) + half_size(3);
    ];

    faces = [
        1 2 3 4;  % bottom
        5 6 7 8;  % top
        1 2 6 5;  % front
        3 4 8 7;  % back
        1 4 8 5;  % left
        2 3 7 6;  % right
    ];

    patch(ax, 'Vertices', vertices, 'Faces', faces, ...
          'FaceColor', color, 'FaceAlpha', alpha, 'EdgeColor', 'k');
end
