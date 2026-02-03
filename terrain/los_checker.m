classdef los_checker < handle
%LOS_CHECKER Line-of-sight checker for terrain masking analysis
%
% Determines whether there is a clear line of sight between two points,
% considering terrain blockage. Essential for radar detection modeling.
%
% Methods:
%   has_los(p1, p2)           - Check if LOS exists between two points
%   get_los_fraction(p1, p2)  - Get fraction of path that has LOS
%   find_visible_cells(pos)   - Find all visible cells from a position
%   get_masking_angle(pos, dir) - Get terrain masking angle in direction
%
% Example:
%   terrain_data = terrain_generator('ridge');
%   tm = terrain_map(terrain_data);
%   los = los_checker(tm);
%
%   radar_pos = [0; 0; 200];  % Radar on hill
%   uav_pos = [500; 100; 50]; % UAV behind ridge
%   if los.has_los(radar_pos, uav_pos)
%       disp('UAV is visible to radar');
%   else
%       disp('UAV is masked by terrain');
%   end
%
% Author: Quadrotor Terrain Following Project

    properties
        terrain         % terrain_map object
        sample_spacing  % Spacing between LOS check samples [m]
        clearance       % Minimum clearance above terrain for LOS [m]
    end

    methods
        function obj = los_checker(terrain_map_obj, sample_spacing, clearance)
            %LOS_CHECKER Constructor
            %   los = los_checker(terrain_map) - Use default spacing
            %   los = los_checker(terrain_map, spacing, clearance)

            obj.terrain = terrain_map_obj;

            if nargin < 2
                obj.sample_spacing = terrain_map_obj.resolution;
            else
                obj.sample_spacing = sample_spacing;
            end

            if nargin < 3
                obj.clearance = 1.0;  % 1m default clearance
            else
                obj.clearance = clearance;
            end
        end

        function [visible, blocked_at] = has_los(obj, p1, p2)
            %HAS_LOS Check line of sight between two points
            %   [visible, blocked_at] = has_los(p1, p2)
            %
            %   p1, p2 are [3x1] position vectors [N; E; altitude]
            %   where altitude is positive upward (not NED D!)
            %
            %   visible: true if clear LOS exists
            %   blocked_at: [3x1] position where LOS is blocked (NaN if visible)

            p1 = p1(:);
            p2 = p2(:);

            % Distance between points
            dist = norm(p2 - p1);

            if dist < obj.sample_spacing
                visible = true;
                blocked_at = nan(3, 1);
                return;
            end

            % Number of samples along the line
            n_samples = ceil(dist / obj.sample_spacing) + 1;

            % Sample points along the line
            t = linspace(0, 1, n_samples);
            points = p1 + (p2 - p1) * t;  % [3 x n_samples]

            N = points(1, :);
            E = points(2, :);
            alt = points(3, :);  % Altitude (positive up)

            % Get terrain height at each sample
            terrain_h = obj.terrain.get_height(N, E);

            % Check if any point is below terrain (with clearance)
            clearance_vec = obj.clearance * ones(size(alt));
            % Don't require clearance at endpoints (they might be on terrain)
            clearance_vec(1) = 0;
            clearance_vec(end) = 0;

            blocked_mask = alt < (terrain_h + clearance_vec);

            if any(blocked_mask)
                visible = false;
                % Find first blocked point
                idx = find(blocked_mask, 1, 'first');
                blocked_at = points(:, idx);
            else
                visible = true;
                blocked_at = nan(3, 1);
            end
        end

        function fraction = get_los_fraction(obj, p1, p2)
            %GET_LOS_FRACTION Get fraction of path with clear LOS
            %   fraction = get_los_fraction(p1, p2)
            %
            %   Returns value in [0, 1] where 1 = fully visible

            p1 = p1(:);
            p2 = p2(:);

            dist = norm(p2 - p1);
            if dist < obj.sample_spacing
                fraction = 1.0;
                return;
            end

            n_samples = ceil(dist / obj.sample_spacing) + 1;
            t = linspace(0, 1, n_samples);
            points = p1 + (p2 - p1) * t;

            N = points(1, :);
            E = points(2, :);
            alt = points(3, :);

            terrain_h = obj.terrain.get_height(N, E);
            visible_mask = alt >= (terrain_h + obj.clearance);

            fraction = sum(visible_mask) / n_samples;
        end

        function [visible_grid, N_grid, E_grid] = find_visible_area(obj, observer_pos, grid_alt, bounds)
            %FIND_VISIBLE_AREA Find all areas visible from observer position
            %   [visible, N, E] = find_visible_area(pos, alt, bounds)
            %
            %   observer_pos: [3x1] position [N; E; altitude]
            %   grid_alt: altitude to check visibility at
            %   bounds: [N_min, N_max, E_min, E_max] (default: terrain bounds)
            %
            %   Returns logical grid of visible cells

            if nargin < 4
                bounds = obj.terrain.bounds;
            end

            N_vec = bounds(1):obj.terrain.resolution:bounds(2);
            E_vec = bounds(3):obj.terrain.resolution:bounds(4);
            [N_grid, E_grid] = meshgrid(N_vec, E_vec);

            visible_grid = false(size(N_grid));

            for i = 1:numel(N_grid)
                target_pos = [N_grid(i); E_grid(i); grid_alt];
                visible_grid(i) = obj.has_los(observer_pos, target_pos);
            end
        end

        function angle = get_masking_angle(obj, pos, direction, range)
            %GET_MASKING_ANGLE Get terrain masking angle in a direction
            %   angle = get_masking_angle(pos, direction, range)
            %
            %   pos: [3x1] observer position [N; E; altitude]
            %   direction: [2x1] horizontal direction [dN; dE] (will be normalized)
            %   range: maximum range to check [m]
            %
            %   Returns elevation angle below which terrain blocks view
            %   (positive angle = above horizontal, negative = below)

            if nargin < 4
                range = 1000;
            end

            pos = pos(:);
            direction = direction(:);
            direction = direction / norm(direction);

            % Sample along the direction
            n_samples = ceil(range / obj.sample_spacing);
            distances = linspace(obj.sample_spacing, range, n_samples);

            max_angle = -90;  % Start with looking straight down

            for d = distances
                target_N = pos(1) + direction(1) * d;
                target_E = pos(2) + direction(2) * d;
                terrain_h = obj.terrain.get_height(target_N, target_E);

                if ~isnan(terrain_h)
                    % Elevation angle to terrain at this distance
                    elev_diff = terrain_h - pos(3);
                    elev_angle = atan2d(elev_diff, d);
                    max_angle = max(max_angle, elev_angle);
                end
            end

            angle = max_angle;
        end

        function fig = plot_los_analysis(obj, observer_pos, grid_alt, bounds)
            %PLOT_LOS_ANALYSIS Visualize line-of-sight coverage
            %   plot_los_analysis(pos, alt) - Plot visible area from position

            if nargin < 4
                bounds = obj.terrain.bounds;
            end

            [visible_grid, N_grid, E_grid] = obj.find_visible_area(observer_pos, grid_alt, bounds);

            fig = figure('Name', 'LOS Analysis', 'Position', [100, 100, 1000, 500]);

            % 2D visibility map
            subplot(1, 2, 1);
            imagesc(obj.terrain.N_vec, obj.terrain.E_vec, double(visible_grid));
            hold on;
            plot(observer_pos(1), observer_pos(2), 'rp', 'MarkerSize', 20, 'MarkerFaceColor', 'r');
            colormap(gca, [0.3 0.3 0.3; 0.2 0.8 0.2]);  % Dark = blocked, Green = visible
            xlabel('North [m]');
            ylabel('East [m]');
            title(sprintf('Visibility at %.0fm altitude', grid_alt));
            axis equal;
            colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Blocked', 'Visible'});

            % 3D terrain with visibility overlay
            subplot(1, 2, 2);
            [N_mesh, E_mesh] = meshgrid(obj.terrain.N_vec, obj.terrain.E_vec);
            surf(N_mesh, E_mesh, obj.terrain.Z, double(visible_grid), ...
                 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            hold on;
            plot3(observer_pos(1), observer_pos(2), observer_pos(3), ...
                  'rp', 'MarkerSize', 20, 'MarkerFaceColor', 'r');
            colormap(gca, [0.5 0.3 0.2; 0.2 0.8 0.2]);
            xlabel('North [m]');
            ylabel('East [m]');
            zlabel('Height [m]');
            title('3D Visibility');
            view(30, 45);
            axis equal;
        end

        function [angles, ranges] = compute_horizon_profile(obj, pos, n_azimuths, max_range)
            %COMPUTE_HORIZON_PROFILE Compute terrain horizon in all directions
            %   [angles, ranges] = compute_horizon_profile(pos, n_azimuths, max_range)
            %
            %   Returns masking angle for each azimuth direction

            if nargin < 3
                n_azimuths = 36;  % Every 10 degrees
            end
            if nargin < 4
                max_range = 1000;
            end

            azimuths = linspace(0, 360, n_azimuths + 1);
            azimuths = azimuths(1:end-1);  % Remove duplicate 360

            angles = zeros(size(azimuths));
            ranges = max_range * ones(size(azimuths));

            for i = 1:length(azimuths)
                az_rad = deg2rad(azimuths(i));
                direction = [cos(az_rad); sin(az_rad)];
                angles(i) = obj.get_masking_angle(pos, direction, max_range);
            end
        end

        function [intersects, intersection_pt, distance] = cast_radar_ray(obj, origin, direction, max_range)
            %CAST_RADAR_RAY Cast a radar beam ray and find terrain intersection
            %   [intersects, point, dist] = cast_radar_ray(origin, direction, max_range)
            %
            %   Casts a ray from origin in the given direction and finds where
            %   it intersects the terrain. Uses stepping with binary search
            %   refinement for precise intersection detection.
            %
            %   Inputs:
            %       origin    - [3x1] ray origin [N; E; altitude] (altitude positive up)
            %       direction - [3x1] ray direction vector (will be normalized)
            %       max_range - Maximum range to check [m]
            %
            %   Outputs:
            %       intersects      - true if ray intersects terrain
            %       intersection_pt - [3x1] intersection point (or endpoint if no hit)
            %       distance        - Distance from origin to intersection
            %
            %   Example:
            %       radar_pos = [0; 0; 200];
            %       direction = [1; 0; -0.1];  % Slightly downward beam
            %       [hit, pt, d] = los.cast_radar_ray(radar_pos, direction, 5000);

            origin = origin(:);
            direction = direction(:);

            % Normalize direction
            dir_norm = norm(direction);
            if dir_norm < 1e-10
                intersects = false;
                intersection_pt = origin;
                distance = 0;
                return;
            end
            direction = direction / dir_norm;

            % Use half the terrain resolution for stepping
            step = obj.sample_spacing / 2;
            n_steps = ceil(max_range / step);

            prev_above = true;
            precision = 0.5;  % Binary search precision [m]

            for i = 1:n_steps
                t = i * step;
                point = origin + direction * t;

                terrain_h = obj.terrain.get_height(point(1), point(2));

                if isnan(terrain_h)
                    continue;
                end

                current_above = point(3) >= terrain_h;

                if ~current_above && prev_above
                    % Ray crossed terrain - refine with binary search
                    t_low = (i - 1) * step;
                    t_high = t;

                    while (t_high - t_low) > precision
                        t_mid = (t_low + t_high) / 2;
                        p_mid = origin + direction * t_mid;
                        h_mid = obj.terrain.get_height(p_mid(1), p_mid(2));

                        if isnan(h_mid) || p_mid(3) >= h_mid
                            t_low = t_mid;
                        else
                            t_high = t_mid;
                        end
                    end

                    distance = (t_low + t_high) / 2;
                    intersection_pt = origin + direction * distance;
                    intersects = true;
                    return;
                end

                prev_above = current_above;
            end

            % No intersection
            intersects = false;
            intersection_pt = origin + direction * max_range;
            distance = max_range;
        end

        function [visible_from_radar, shadow_map] = compute_radar_shadow(obj, radar_pos, grid_alt, bounds)
            %COMPUTE_RADAR_SHADOW Compute radar shadow zones at flight altitude
            %   [visible, shadow] = compute_radar_shadow(radar_pos, alt, bounds)
            %
            %   Computes which areas at the specified flight altitude are
            %   visible to the radar and which are in terrain shadow.
            %
            %   Inputs:
            %       radar_pos - [3x1] radar position [N; E; altitude]
            %       grid_alt  - Flight altitude to check visibility [m]
            %       bounds    - [N_min, N_max, E_min, E_max] (default: terrain bounds)
            %
            %   Outputs:
            %       visible_from_radar - Logical grid (true = visible)
            %       shadow_map         - Logical grid (true = in shadow)

            if nargin < 4
                bounds = obj.terrain.bounds;
            end

            radar_pos = radar_pos(:);

            N_vec = bounds(1):obj.terrain.resolution:bounds(2);
            E_vec = bounds(3):obj.terrain.resolution:bounds(4);
            [N_grid, E_grid] = meshgrid(N_vec, E_vec);

            visible_from_radar = false(size(N_grid));

            for i = 1:numel(N_grid)
                target_pos = [N_grid(i); E_grid(i); grid_alt];
                visible_from_radar(i) = obj.has_los(radar_pos, target_pos);
            end

            shadow_map = ~visible_from_radar;
        end
    end
end
