classdef radar_ray_caster < handle
%RADAR_RAY_CASTER Ray casting for radar beam-terrain intersection analysis
%
% Performs ray casting from radar positions to determine terrain
% intersections, shadow zones, and visibility coverage. Supports two modes:
%   1. Grid-based stepping with binary search refinement (fallback)
%   2. Mesh-based Möller-Trumbore ray-triangle intersection (preferred)
%
% The mesh-based method is more accurate and faster with spatial indexing.
%
% Properties:
%   terrain         - terrain_map object for height queries
%   mesh            - terrain_mesh object for triangle intersection
%   use_mesh        - Use mesh-based intersection (bool)
%   step_size       - Ray stepping distance [m] (grid mode)
%   precision       - Binary search precision for intersection [m]
%   earth_curvature - Enable Earth curvature correction
%
% Methods:
%   cast_ray(origin, direction, max_range) - Single ray intersection
%   cast_ray_to_target(origin, target)     - Ray toward specific point
%   compute_shadow_map(radar_pos, bounds)  - 2D shadow zone map
%   compute_coverage_map(radar_pos, bounds, alt) - Visibility at altitude
%   get_horizon_profile(pos, n_azimuths)   - 360-degree horizon
%   plot_shadow_map(radar_pos, bounds)     - Visualize shadow zones
%
% Example:
%   terrain_data = terrain_generator(struct('type','hills'));
%   tm = terrain_map(terrain_data);
%   rc = radar_ray_caster(tm);  % Auto-creates mesh for fast intersection
%
%   radar_pos = [0; 0; 200];  % Radar position [N; E; altitude]
%   bounds = [-1000 1000 -1000 1000];
%   [hit, point, dist] = rc.cast_ray(radar_pos, [1; 0; -0.1], 5000);
%   if hit
%       fprintf('Beam hits terrain at [%.1f, %.1f, %.1f]\n', point);
%   end
%
%   % Compute shadow map
%   shadow = rc.compute_shadow_map(radar_pos, bounds);
%   rc.plot_shadow_map(radar_pos, bounds);
%   
%   % Compute visibility slice
%   alt = 200; % [m] 
%   azimuth = 90; %[deg]
%   max_range = 1000; %[m]
%   resolution = 30; %[m]
%   coverage = rc.compute_coverage_map(radar_pos,bounds,alt);
%   rc.plot_coverage_slice(radar_pos,azimuth,max_range,resolution);
%
% Author: Quadrotor Terrain Following Project
% See also: terrain_map, terrain_mesh, moller_trumbore, los_checker

    properties
        terrain         % terrain_map object
        mesh            % terrain_mesh object for ray-triangle intersection
        use_mesh        % Use mesh-based intersection (bool)
        step_size       % Ray stepping distance [m]
        precision       % Binary search precision [m]
        earth_curvature % Enable Earth curvature correction (bool)
        earth_radius    % Earth radius for curvature correction [m]
    end

    methods
        function obj = radar_ray_caster(terrain_map_obj, params)
            %RADAR_RAY_CASTER Constructor
            %   rc = radar_ray_caster(terrain_map)
            %   rc = radar_ray_caster(terrain_map, params)
            %
            %   params struct fields:
            %     .step_size  - Ray stepping distance [m] (default: resolution/2)
            %     .precision  - Intersection precision [m] (default: 0.5)
            %     .use_mesh   - Use mesh-based intersection (default: true)
            %     .mesh       - Pre-built terrain_mesh (optional)

            obj.terrain = terrain_map_obj;

            % Handle legacy call with step_size as second argument
            if nargin >= 2 && isnumeric(params)
                step_size = params;
                params = struct();
                params.step_size = step_size;
            elseif nargin < 2
                params = struct();
            end

            % Set step size
            if isfield(params, 'step_size') && ~isempty(params.step_size)
                obj.step_size = params.step_size;
            else
                obj.step_size = terrain_map_obj.resolution / 2;
            end

            % Set precision
            if isfield(params, 'precision') && ~isempty(params.precision)
                obj.precision = params.precision;
            else
                obj.precision = 0.5;  % 0.5m default precision
            end

            obj.earth_curvature = false;
            obj.earth_radius = 6371000;  % meters

            % Initialize mesh-based intersection
            if isfield(params, 'use_mesh')
                obj.use_mesh = params.use_mesh;
            else
                obj.use_mesh = true;  % Default to mesh-based (more accurate)
            end

            if obj.use_mesh
                if isfield(params, 'mesh') && ~isempty(params.mesh)
                    obj.mesh = params.mesh;
                else
                    % Create terrain mesh from terrain_map data
                    fprintf('Building terrain mesh for ray casting...\n');
                    try
                        terrain_data = struct();
                        terrain_data.N_vec = terrain_map_obj.N_vec;
                        terrain_data.E_vec = terrain_map_obj.E_vec;
                        terrain_data.Z = terrain_map_obj.Z;
                        terrain_data.bounds = terrain_map_obj.bounds;
                        terrain_data.resolution = terrain_map_obj.resolution;

                        obj.mesh = terrain_mesh(terrain_data);
                        fprintf('Terrain mesh created: %d triangles\n', size(obj.mesh.triangles, 1));
                    catch ME
                        warning('radar_ray_caster:MeshFailed', ...
                            'Failed to create terrain mesh: %s. Using grid-based method.', ME.message);
                        obj.use_mesh = false;
                        obj.mesh = [];
                    end
                end
            else
                obj.mesh = [];
            end
        end

        function [hit, intersection_point, distance] = cast_ray(obj, origin, direction, max_range)
            %CAST_RAY Cast a single ray and find terrain intersection
            %   [hit, point, dist] = cast_ray(origin, direction, max_range)
            %
            %   Inputs:
            %       origin    - [3x1] ray origin [N; E; altitude] (altitude positive up)
            %       direction - [3x1] ray direction (will be normalized)
            %       max_range - Maximum range to check [m]
            %
            %   Outputs:
            %       hit              - true if ray intersects terrain
            %       intersection_point - [3x1] intersection point or end point
            %       distance         - Distance to intersection or max_range

            origin = origin(:);
            direction = direction(:);

            % Normalize direction
            dir_norm = norm(direction);
            if dir_norm < 1e-10
                hit = false;
                intersection_point = origin;
                distance = 0;
                return;
            end
            direction = direction / dir_norm;

            % Use mesh-based intersection if available (faster and more accurate)
            if obj.use_mesh && ~isempty(obj.mesh)
                [hit, intersection_point, distance] = obj.cast_ray_mesh(origin, direction, max_range);
                return;
            end

            % Fallback to grid-based stepping method
            [hit, intersection_point, distance] = obj.cast_ray_grid(origin, direction, max_range);
        end

        function [hit, intersection_point, distance] = cast_ray_mesh(obj, origin, direction, max_range)
            %CAST_RAY_MESH Ray casting using mesh-based Möller-Trumbore algorithm
            %   Uses terrain_mesh.ray_intersect for precise intersection

            [hit, intersection_point, distance, ~] = obj.mesh.ray_intersect(origin, direction, max_range);

            if ~hit
                intersection_point = origin + direction * max_range;
                distance = max_range;
            end
        end

        function [hit, intersection_point, distance] = cast_ray_grid(obj, origin, direction, max_range)
            %CAST_RAY_GRID Grid-based ray casting with binary search refinement
            %   Original stepping algorithm - used as fallback

            % Number of steps
            n_steps = ceil(max_range / obj.step_size);

            % Step along ray
            prev_above = true;  % Track if previous point was above terrain

            for i = 1:n_steps
                t = i * obj.step_size;
                point = origin + direction * t;

                % Get terrain height at current position
                terrain_h = obj.terrain.get_height(point(1), point(2));

                % Apply Earth curvature correction if enabled
                if obj.earth_curvature
                    horiz_dist = sqrt((point(1) - origin(1))^2 + (point(2) - origin(2))^2);
                    curvature_drop = horiz_dist^2 / (2 * obj.earth_radius);
                    terrain_h = terrain_h + curvature_drop;
                end

                % Check for out of bounds
                if isnan(terrain_h)
                    continue;
                end

                % Check if ray is below terrain
                current_above = point(3) >= terrain_h;

                if ~current_above && prev_above
                    % Ray crossed from above to below terrain - intersection!
                    % Use binary search to refine intersection point
                    t_low = (i - 1) * obj.step_size;
                    t_high = t;
                    [intersection_point, distance] = obj.binary_search_intersection(...
                        origin, direction, t_low, t_high);
                    hit = true;
                    return;
                end

                prev_above = current_above;
            end

            % No intersection found
            hit = false;
            intersection_point = origin + direction * max_range;
            distance = max_range;
        end

        function [hit, intersection_point, distance] = cast_ray_to_target(obj, origin, target)
            %CAST_RAY_TO_TARGET Cast ray from origin toward target
            %   [hit, point, dist] = cast_ray_to_target(origin, target)
            %
            %   Checks if there is terrain between origin and target.
            %   If hit is true, intersection_point is where beam hits terrain.
            %   If hit is false, the beam reaches the target unobstructed.

            origin = origin(:);
            target = target(:);

            direction = target - origin;
            max_range = norm(direction);

            if max_range < obj.precision
                hit = false;
                intersection_point = target;
                distance = max_range;
                return;
            end

            [hit, intersection_point, distance] = obj.cast_ray(origin, direction, max_range);

            % If no hit before target, ray reaches target
            if ~hit || distance >= max_range - obj.precision
                hit = false;
                intersection_point = target;
                distance = max_range;
            end
        end

        function [shadow_map, N_grid, E_grid] = compute_shadow_map(obj, radar_pos, bounds, resolution)
            %COMPUTE_SHADOW_MAP Compute 2D shadow zone map on terrain
            %   [shadow, N, E] = compute_shadow_map(radar_pos, bounds)
            %   [shadow, N, E] = compute_shadow_map(radar_pos, bounds, resolution)
            %
            %   Inputs:
            %       radar_pos  - [3x1] radar position [N; E; altitude]
            %       bounds     - [N_min, N_max, E_min, E_max]
            %       resolution - Grid resolution [m] (default: terrain resolution)
            %
            %   Outputs:
            %       shadow_map - Logical grid (true = in shadow, false = visible)
            %       N_grid, E_grid - Meshgrid coordinates

            if nargin < 4
                resolution = obj.terrain.resolution;
            end

            radar_pos = radar_pos(:);

            N_vec = bounds(1):resolution:bounds(2);
            E_vec = bounds(3):resolution:bounds(4);
            [N_grid, E_grid] = meshgrid(N_vec, E_vec);

            shadow_map = true(size(N_grid));  % Default to shadow

            fprintf('Computing shadow map (%d x %d points)...\n', length(N_vec), length(E_vec));

            for i = 1:numel(N_grid)
                % Target point on terrain surface
                N = N_grid(i);
                E = E_grid(i);
                terrain_h = obj.terrain.get_height(N, E);

                if isnan(terrain_h)
                    continue;
                end

                target = [N; E; terrain_h + 1];  % 1m above terrain

                % Cast ray from radar to target
                [hit, ~, ~] = obj.cast_ray_to_target(radar_pos, target);

                % If ray doesn't hit terrain before target, point is visible
                shadow_map(i) = hit;
            end

            fprintf('Shadow map complete. Shadow coverage: %.1f%%\n', ...
                100 * sum(shadow_map(:)) / numel(shadow_map));
        end

        function [visibility_map, N_grid, E_grid] = compute_coverage_map(obj, radar_pos, bounds, flight_alt, resolution)
            %COMPUTE_COVERAGE_MAP Compute visibility at specified altitude
            %   [vis, N, E] = compute_coverage_map(radar_pos, bounds, alt)
            %
            %   Computes radar visibility for aircraft flying at constant altitude.
            %
            %   Inputs:
            %       radar_pos  - [3x1] radar position
            %       bounds     - [N_min, N_max, E_min, E_max]
            %       flight_alt - Flight altitude [m] (positive up)
            %       resolution - Grid resolution [m]
            %
            %   Outputs:
            %       visibility_map - Logical grid (true = visible to radar)

            if nargin < 5
                resolution = obj.terrain.resolution;
            end

            radar_pos = radar_pos(:);

            N_vec = bounds(1):resolution:bounds(2);
            E_vec = bounds(3):resolution:bounds(4);
            [N_grid, E_grid] = meshgrid(N_vec, E_vec);

            visibility_map = false(size(N_grid));

            fprintf('Computing coverage map at %.0fm altitude...\n', flight_alt);

            for i = 1:numel(N_grid)
                target = [N_grid(i); E_grid(i); flight_alt];

                % Cast ray from radar to target
                [hit, ~, ~] = obj.cast_ray_to_target(radar_pos, target);

                % Visible if ray doesn't hit terrain
                visibility_map(i) = ~hit;
            end

            fprintf('Coverage map complete. Visible area: %.1f%%\n', ...
                100 * sum(visibility_map(:)) / numel(visibility_map));
        end

        function [horizon_angles, azimuths] = get_horizon_profile(obj, pos, n_azimuths, max_range)
            %GET_HORIZON_PROFILE Compute terrain horizon in all directions
            %   [angles, azimuths] = get_horizon_profile(pos, n_azimuths, max_range)
            %
            %   Returns elevation angle to terrain horizon for each azimuth.
            %   Positive angles = terrain above horizontal plane.

            if nargin < 3
                n_azimuths = 360;
            end
            if nargin < 4
                max_range = 5000;
            end

            pos = pos(:);
            azimuths = linspace(0, 360, n_azimuths + 1);
            azimuths = azimuths(1:end-1);
            horizon_angles = zeros(size(azimuths));

            for i = 1:length(azimuths)
                az_rad = deg2rad(azimuths(i));
                % Direction in horizontal plane (N-E)
                dir_horiz = [cos(az_rad); sin(az_rad)];

                % Find maximum elevation angle to terrain
                max_elev = -90;

                for r = obj.step_size:obj.step_size:max_range
                    N = pos(1) + dir_horiz(1) * r;
                    E = pos(2) + dir_horiz(2) * r;
                    h = obj.terrain.get_height(N, E);

                    if ~isnan(h)
                        elev_angle = atan2d(h - pos(3), r);
                        if elev_angle > max_elev
                            max_elev = elev_angle;
                        end
                    end
                end

                horizon_angles(i) = max_elev;
            end
        end

        function fig = plot_shadow_map(obj, radar_pos, bounds, resolution)
            %PLOT_SHADOW_MAP Visualize shadow zones on terrain
            %   fig = plot_shadow_map(radar_pos, bounds)

            if nargin < 4
                resolution = obj.terrain.resolution * 2;  % Coarser for speed
            end

            [shadow_map, N_grid, E_grid] = obj.compute_shadow_map(radar_pos, bounds, resolution);

            fig = figure('Name', 'Radar Shadow Map', 'Position', [100, 100, 1200, 500]);

            % 2D view
            subplot(1, 2, 1);
            % Get terrain heights for visualization
            Z_grid = zeros(size(N_grid));
            for i = 1:numel(N_grid)
                Z_grid(i) = obj.terrain.get_height(N_grid(i), E_grid(i));
            end

            % Plot terrain with shadow overlay
            imagesc(bounds(1:2), bounds(3:4), double(~shadow_map));
            hold on;
            contour(N_grid, E_grid, Z_grid, 15, 'k', 'LineWidth', 0.5);
            plot(radar_pos(1), radar_pos(2), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

            colormap(gca, [0.7 0.7 0.7; 0.2 0.8 0.2]);  % Gray = shadow, Green = visible
            xlabel('North [m]');
            ylabel('East [m]');
            title(sprintf('Radar Shadow Map (Radar at %.0fm alt)', radar_pos(3)));
            axis equal;
            colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Shadow', 'Visible'});

            % 3D view
            subplot(1, 2, 2);
            surf(N_grid, E_grid, Z_grid, double(~shadow_map), ...
                 'EdgeColor', 'none', 'FaceAlpha', 0.9);
            hold on;
            plot3(radar_pos(1), radar_pos(2), radar_pos(3), ...
                  'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

            colormap(gca, [0.5 0.3 0.2; 0.2 0.8 0.2]);
            xlabel('North [m]');
            ylabel('East [m]');
            zlabel('Height [m]');
            title('3D Shadow Visualization');
            view(30, 45);
            axis equal;
            light('Position', [1 1 1]);
            lighting gouraud;
        end

        function fig = plot_coverage_slice(obj, radar_pos, azimuth, max_range, resolution)
            %PLOT_COVERAGE_SLICE Plot vertical slice of radar coverage
            %   fig = plot_coverage_slice(radar_pos, azimuth, max_range, resolution)
            %
            %   Shows terrain profile and radar beam coverage in vertical plane.

            if nargin < 4
                max_range = 5000;
            end
            if nargin < 5
                resolution = obj.terrain.resolution;
            end

            radar_pos = radar_pos(:);
            az_rad = deg2rad(azimuth);
            dir_horiz = [cos(az_rad); sin(az_rad)];

            % Sample along azimuth
            ranges = 0:resolution:max_range;
            terrain_profile = zeros(size(ranges));

            for i = 1:length(ranges)
                N = radar_pos(1) + dir_horiz(1) * ranges(i);
                E = radar_pos(2) + dir_horiz(2) * ranges(i);
                h = obj.terrain.get_height(N, E);
                if isnan(h)
                    terrain_profile(i) = 0;
                else
                    terrain_profile(i) = h;
                end
            end

            fig = figure('Name', sprintf('Coverage Slice (Az=%.0f deg)', azimuth), ...
                         'Position', [100, 100, 1000, 400]);

            % Plot terrain profile
            area(ranges, terrain_profile, 'FaceColor', [0.6 0.4 0.2], 'EdgeColor', 'k');
            hold on;

            % Plot radar position
            plot(0, radar_pos(3), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

            % Plot LOS rays at various elevations
            elevations = -30:5:30;  % degrees
            colors = jet(length(elevations));

            for i = 1:length(elevations)
                el_rad = deg2rad(elevations(i));
                dir_3d = [dir_horiz(1) * cos(el_rad); ...
                          dir_horiz(2) * cos(el_rad); ...
                          sin(el_rad)];

                [hit, int_pt, dist] = obj.cast_ray(radar_pos, dir_3d, max_range);

                % Project intersection to range
                if hit
                    r_hit = sqrt((int_pt(1) - radar_pos(1))^2 + (int_pt(2) - radar_pos(2))^2);
                    plot([0, r_hit], [radar_pos(3), int_pt(3)], '-', ...
                         'Color', colors(i,:), 'LineWidth', 1);
                    plot(r_hit, int_pt(3), 'o', 'Color', colors(i,:), 'MarkerSize', 5);
                else
                    r_end = max_range * cos(el_rad);
                    h_end = radar_pos(3) + max_range * sin(el_rad);
                    plot([0, r_end], [radar_pos(3), h_end], '--', ...
                         'Color', colors(i,:), 'LineWidth', 0.5);
                end
            end

            xlabel('Range [m]');
            ylabel('Altitude [m]');
            title(sprintf('Radar Coverage Slice (Azimuth = %.0f deg)', azimuth));
            grid on;
            ylim([min(terrain_profile) - 50, max(radar_pos(3), max(terrain_profile)) + 200]);
        end
    end

    methods (Access = private)
        function [point, distance] = binary_search_intersection(obj, origin, direction, t_low, t_high)
            %BINARY_SEARCH_INTERSECTION Refine intersection point

            while (t_high - t_low) > obj.precision
                t_mid = (t_low + t_high) / 2;
                p_mid = origin + direction * t_mid;

                terrain_h = obj.terrain.get_height(p_mid(1), p_mid(2));

                if obj.earth_curvature
                    horiz_dist = sqrt((p_mid(1) - origin(1))^2 + (p_mid(2) - origin(2))^2);
                    curvature_drop = horiz_dist^2 / (2 * obj.earth_radius);
                    terrain_h = terrain_h + curvature_drop;
                end

                if isnan(terrain_h) || p_mid(3) >= terrain_h
                    t_low = t_mid;
                else
                    t_high = t_mid;
                end
            end

            distance = (t_low + t_high) / 2;
            point = origin + direction * distance;
        end
    end
end
