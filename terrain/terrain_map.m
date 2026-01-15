classdef terrain_map < handle
%TERRAIN_MAP Class for terrain height queries and visualization
%
% Provides efficient terrain height lookups using interpolation.
% Can load from terrain_generator output or external DEM data.
%
% Properties:
%   N_vec, E_vec  - Grid vectors
%   Z             - Height matrix
%   bounds        - [N_min, N_max, E_min, E_max]
%   resolution    - Grid spacing
%   interp_method - Interpolation method ('linear', 'cubic', 'spline')
%
% Methods:
%   get_height(N, E)           - Get terrain height at point(s)
%   get_height_and_normal(N,E) - Get height and surface normal
%   get_clearance(N, E, alt)   - Get altitude above terrain
%   is_underground(N, E, alt)  - Check if point is below terrain
%   plot()                     - Visualize terrain
%   plot_with_path(path)       - Plot terrain with flight path
%
% Example:
%   terrain_data = terrain_generator('hills');
%   tm = terrain_map(terrain_data);
%   h = tm.get_height(100, 50);  % Height at N=100, E=50
%   tm.plot();
%
% Author: Quadrotor Terrain Following Project

    properties
        N_vec           % North coordinate vector
        E_vec           % East coordinate vector
        Z               % Height matrix [M x N] where M = length(E_vec), N = length(N_vec)
        bounds          % [N_min, N_max, E_min, E_max]
        resolution      % Grid resolution
        interp_method   % Interpolation method
        type            % Terrain type description
    end

    properties (Access = private)
        F_interp        % Interpolant object for fast queries
    end

    methods
        function obj = terrain_map(terrain_data, interp_method)
            %TERRAIN_MAP Constructor
            %   tm = terrain_map(terrain_data) - Load from terrain_generator output
            %   tm = terrain_map(terrain_data, 'cubic') - Specify interpolation

            if nargin < 2
                interp_method = 'linear';
            end

            obj.N_vec = terrain_data.N_vec;
            obj.E_vec = terrain_data.E_vec;
            obj.Z = terrain_data.Z;
            obj.bounds = terrain_data.bounds;
            obj.resolution = terrain_data.resolution;
            obj.interp_method = interp_method;

            if isfield(terrain_data, 'type')
                obj.type = terrain_data.type;
            else
                obj.type = 'unknown';
            end

            % Create interpolant for fast queries
            % Note: griddedInterpolant expects (row, col) = (E, N) ordering
            obj.F_interp = griddedInterpolant({obj.E_vec, obj.N_vec}, obj.Z, ...
                                               interp_method, 'nearest');
        end

        function h = get_height(obj, N, E)
            %GET_HEIGHT Get terrain height at specified coordinates
            %   h = get_height(N, E) - Returns height(s) at position(s)
            %
            %   N, E can be scalars, vectors, or matrices of same size.
            %   Returns NaN for points outside terrain bounds.

            % Handle both scalar and array inputs
            N = N(:);
            E = E(:);

            % Check bounds
            in_bounds = N >= obj.bounds(1) & N <= obj.bounds(2) & ...
                        E >= obj.bounds(3) & E <= obj.bounds(4);

            h = nan(size(N));
            if any(in_bounds)
                h(in_bounds) = obj.F_interp(E(in_bounds), N(in_bounds));
            end
        end

        function [h, normal] = get_height_and_normal(obj, N, E)
            %GET_HEIGHT_AND_NORMAL Get terrain height and surface normal
            %   [h, normal] = get_height_and_normal(N, E)
            %
            %   normal is [3 x n] matrix with unit normals pointing up

            h = obj.get_height(N, E);

            % Compute normal using finite differences
            delta = obj.resolution / 2;
            h_n = obj.get_height(N + delta, E);
            h_s = obj.get_height(N - delta, E);
            h_e = obj.get_height(N, E + delta);
            h_w = obj.get_height(N, E - delta);

            % Gradient
            dh_dN = (h_n - h_s) / (2 * delta);
            dh_dE = (h_e - h_w) / (2 * delta);

            % Surface normal: n = [-dh/dN, -dh/dE, 1] normalized
            % In NED frame, positive Z is down, but terrain height is positive up
            % So normal pointing "up" in world frame is [-dh/dN, -dh/dE, 1]
            n_N = -dh_dN(:);
            n_E = -dh_dE(:);
            n_D = -ones(size(n_N));  % In NED, up is -D

            mag = sqrt(n_N.^2 + n_E.^2 + n_D.^2);
            normal = [n_N ./ mag, n_E ./ mag, n_D ./ mag]';
        end

        function clearance = get_clearance(obj, N, E, alt)
            %GET_CLEARANCE Get altitude above ground level (AGL)
            %   clearance = get_clearance(N, E, alt)
            %
            %   alt is altitude in NED frame (negative = above ground)
            %   clearance is positive if above terrain, negative if below

            h = obj.get_height(N, E);
            % In NED: altitude = -D, terrain height = h (positive up)
            % Clearance = terrain_height - (-altitude) = h + alt
            % But alt in NED is negative when above ground, so:
            % If alt = -100 (100m above origin), terrain at h = 50
            % Clearance = 100 - 50 = 50m above terrain
            clearance = -alt(:) - h(:);  % -alt converts NED D to altitude
        end

        function underground = is_underground(obj, N, E, alt)
            %IS_UNDERGROUND Check if point is below terrain surface
            %   underground = is_underground(N, E, alt)
            %
            %   Returns true if point would be inside terrain

            clearance = obj.get_clearance(N, E, alt);
            underground = clearance < 0;
        end

        function [N_path, E_path, Z_path] = get_terrain_profile(obj, N_start, E_start, N_end, E_end, n_points)
            %GET_TERRAIN_PROFILE Get terrain profile along a line
            %   [N, E, Z] = get_terrain_profile(N1, E1, N2, E2, n_points)

            if nargin < 6
                n_points = 100;
            end

            N_path = linspace(N_start, N_end, n_points);
            E_path = linspace(E_start, E_end, n_points);
            Z_path = obj.get_height(N_path, E_path);
        end

        function fig = plot(obj, fig_handle)
            %PLOT Visualize terrain as 3D surface
            %   plot() - Create new figure
            %   plot(fig) - Plot to existing figure

            if nargin < 2
                fig = figure('Name', 'Terrain Map', 'Position', [100, 100, 900, 700]);
            else
                fig = figure(fig_handle);
            end

            % Create mesh grids
            [N_grid, E_grid] = meshgrid(obj.N_vec, obj.E_vec);

            % Plot surface
            surf(N_grid, E_grid, obj.Z, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
            hold on;

            % Formatting
            colormap(terrain_colormap());
            % colorbar('Label', 'Elevation [m]');
            % % Check if the colorbar was added successfully
            % if isempty(findall(fig, 'Type', 'ColorBar'))
            %     disp('Colorbar was not added successfully.');
            % else
            %     disp('Colorbar added successfully.');
            % end
            xlabel('North [m]');
            ylabel('East [m]');
            zlabel('Height [m]');
            title(sprintf('Terrain: %s', obj.type));
            axis equal;
            grid on;
            view(30, 45);

            % Add contour lines at base
            contour3(N_grid, E_grid, obj.Z, 15, 'k', 'LineWidth', 0.5);
        end

        function fig = plot_with_path(obj, path, fig_handle)
            %PLOT_WITH_PATH Plot terrain with flight path overlay
            %   plot_with_path(path) - path is [3 x N] in NED
            %   plot_with_path(path, fig) - Use existing figure

            if nargin < 3
                fig = obj.plot();
            else
                fig = obj.plot(fig_handle);
            end

            hold on;

            % Convert path from NED to terrain coordinates
            % In NED: D is positive down, so altitude = -D
            N = path(1, :);
            E = path(2, :);
            alt = -path(3, :);  % Convert D to altitude above origin

            % Plot path
            plot3(N, E, alt, 'r-', 'LineWidth', 2);

            % Mark start and end
            plot3(N(1), E(1), alt(1), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g');
            plot3(N(end), E(end), alt(end), 'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

            % Plot shadow on terrain
            terrain_h = obj.get_height(N, E);
            plot3(N, E, terrain_h, 'k--', 'LineWidth', 1);

            legend('Terrain', 'Contours', 'Flight Path', 'Start', 'End', 'Ground Track', ...
                   'Location', 'best');
        end

        function fig = plot_2d(obj, fig_handle)
            %PLOT_2D Plot terrain as 2D contour map
            if nargin < 2
                fig = figure('Name', 'Terrain Contour Map', 'Position', [100, 100, 800, 600]);
            else
                fig = figure(fig_handle);
            end

            [N_grid, E_grid] = meshgrid(obj.N_vec, obj.E_vec);

            contourf(N_grid, E_grid, obj.Z, 20);
            colormap(terrain_colormap());
            %colorbar('Label', 'Elevation [m]');
            xlabel('North [m]');
            ylabel('East [m]');
            title(sprintf('Terrain Contour: %s', obj.type));
            axis equal;
            grid on;
        end
    end
end

function cmap = terrain_colormap()
    % Custom colormap for terrain visualization
    % Green (low) -> Brown (mid) -> Gray (high) -> White (peaks)
    colors = [
        0.2, 0.5, 0.2;   % Dark green (valley)
        0.4, 0.6, 0.3;   % Light green
        0.6, 0.5, 0.3;   % Brown
        0.5, 0.4, 0.3;   % Dark brown
        0.6, 0.6, 0.6;   % Gray (rock)
        0.8, 0.8, 0.8;   % Light gray
        1.0, 1.0, 1.0;   % White (snow)
    ];
    n_colors = 256;
    x = linspace(0, 1, size(colors, 1));
    xi = linspace(0, 1, n_colors);
    cmap = interp1(x, colors, xi);
end
