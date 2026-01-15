classdef threat_map < handle
%THREAT_MAP 3D voxel grid of cumulative radar exposure risk
%
% Precomputes and stores radar detection probability throughout the
% operating volume. Supports multiple radars and terrain masking.
%
% Properties:
%   N_vec, E_vec, alt_vec - Grid coordinate vectors
%   risk_grid             - 3D matrix of detection probability
%   radars                - Cell array of radar_site objects
%   terrain               - terrain_map object
%   los                   - los_checker object
%   rcs                   - Target radar cross section [m^2]
%
% Methods:
%   add_radar(radar)      - Add a radar to the threat environment
%   compute_map()         - Calculate risk at all grid points
%   get_risk(N, E, alt)   - Query risk at arbitrary point
%   get_risk_along_path(path) - Get risk along a trajectory
%   plot_slice(alt)       - Visualize risk at constant altitude
%   plot_3d()             - 3D visualization of threat volume
%
% Example:
%   terrain_data = terrain_generator('ridge');
%   tm = terrain_map(terrain_data);
%   los = los_checker(tm);
%
%   threat = threat_map(tm, los, 0.1);  % RCS = 0.1 m^2
%   threat.add_radar(radar_site([0;0;100], 'SAM', 'tracking'));
%   threat.compute_map();
%   threat.plot_slice(50);  % Risk at 50m altitude
%
% Author: Quadrotor Terrain Following Project

    properties
        N_vec               % North coordinate vector
        E_vec               % East coordinate vector
        alt_vec             % Altitude vector (positive up)
        risk_grid           % [M x N x P] detection probability grid
        radars              % Cell array of radar_site objects
        terrain             % terrain_map object
        los                 % los_checker object
        rcs                 % Target RCS [m^2]
        bounds              % [N_min, N_max, E_min, E_max, alt_min, alt_max]
        resolution          % [horiz_res, vert_res] in meters
        computed            % Flag indicating map is computed
    end

    properties (Access = private)
        F_interp            % Interpolant for fast queries
    end

    methods
        function obj = threat_map(terrain_map_obj, los_checker_obj, rcs, params)
            %THREAT_MAP Constructor
            %   threat = threat_map(terrain, los, rcs)
            %   threat = threat_map(terrain, los, rcs, params)
            %
            %   params can include:
            %     .alt_range  - [alt_min, alt_max] (default: [0, 500])
            %     .resolution - [horiz, vert] in meters (default: [50, 25])
            %     .bounds     - [N_min, N_max, E_min, E_max] (default: terrain bounds)

            if nargin < 4
                params = struct();
            end

            obj.terrain = terrain_map_obj;
            obj.los = los_checker_obj;
            obj.rcs = rcs;
            obj.radars = {};
            obj.computed = false;

            % Get bounds
            if isfield(params, 'bounds')
                horiz_bounds = params.bounds;
            else
                horiz_bounds = terrain_map_obj.bounds;
            end

            alt_range = get_param(params, 'alt_range', [0, 500]);
            resolution = get_param(params, 'resolution', [50, 25]);

            obj.bounds = [horiz_bounds, alt_range];
            obj.resolution = resolution;

            % Create coordinate vectors
            obj.N_vec = horiz_bounds(1):resolution(1):horiz_bounds(2);
            obj.E_vec = horiz_bounds(3):resolution(1):horiz_bounds(4);
            obj.alt_vec = alt_range(1):resolution(2):alt_range(2);

            % Initialize empty risk grid
            obj.risk_grid = zeros(length(obj.E_vec), length(obj.N_vec), length(obj.alt_vec));
        end

        function add_radar(obj, radar)
            %ADD_RADAR Add a radar site to the threat environment
            obj.radars{end+1} = radar;
            obj.computed = false;  % Invalidate computed map
        end

        function remove_radar(obj, name_or_index)
            %REMOVE_RADAR Remove a radar by name or index
            if isnumeric(name_or_index)
                obj.radars(name_or_index) = [];
            else
                for i = length(obj.radars):-1:1
                    if strcmp(obj.radars{i}.name, name_or_index)
                        obj.radars(i) = [];
                    end
                end
            end
            obj.computed = false;
        end

        function compute_map(obj, method)
            %COMPUTE_MAP Calculate risk at all grid points
            %   compute_map()        - Use 'max' combination
            %   compute_map('max')   - Max detection probability across radars
            %   compute_map('sum')   - Sum (can exceed 1, useful for cost)
            %   compute_map('prob')  - Probabilistic: 1 - prod(1-P_i)

            if nargin < 2
                method = 'max';
            end

            if isempty(obj.radars)
                warning('No radars added. Risk map will be zero everywhere.');
                obj.computed = true;
                return;
            end

            fprintf('Computing threat map (%d radars)...\n', length(obj.radars));

            n_N = length(obj.N_vec);
            n_E = length(obj.E_vec);
            n_alt = length(obj.alt_vec);
            total_cells = n_N * n_E * n_alt;

            % Initialize grid
            obj.risk_grid = zeros(n_E, n_N, n_alt);

            % Progress counter
            progress_interval = max(1, floor(total_cells / 20));
            cell_count = 0;

            for k = 1:n_alt
                alt = obj.alt_vec(k);

                for j = 1:n_N
                    N = obj.N_vec(j);

                    for i = 1:n_E
                        E = obj.E_vec(i);

                        % Check if point is below terrain
                        terrain_h = obj.terrain.get_height(N, E);
                        if alt < terrain_h
                            obj.risk_grid(i, j, k) = 1.0;  % Inside terrain = maximum risk
                            continue;
                        end

                        target_pos = [N; E; alt];

                        % Compute detection probability from each radar
                        P_det_list = zeros(1, length(obj.radars));
                        for r = 1:length(obj.radars)
                            P_det_list(r) = obj.radars{r}.get_detection_probability(...
                                target_pos, obj.rcs, obj.los);
                        end

                        % Combine probabilities
                        switch lower(method)
                            case 'max'
                                obj.risk_grid(i, j, k) = max(P_det_list);
                            case 'sum'
                                obj.risk_grid(i, j, k) = sum(P_det_list);
                            case 'prob'
                                % Probability of detection by at least one radar
                                obj.risk_grid(i, j, k) = 1 - prod(1 - P_det_list);
                        end

                        cell_count = cell_count + 1;
                        if mod(cell_count, progress_interval) == 0
                            fprintf('  Progress: %d%%\n', round(100*cell_count/total_cells));
                        end
                    end
                end
            end

            % Create interpolant for fast queries
            obj.F_interp = griddedInterpolant({obj.E_vec, obj.N_vec, obj.alt_vec}, ...
                                               obj.risk_grid, 'linear', 'nearest');

            obj.computed = true;
            fprintf('Threat map computed.\n');
        end

        function risk = get_risk(obj, N, E, alt)
            %GET_RISK Query risk at arbitrary point(s)
            %   risk = get_risk(N, E, alt)
            %
            %   Inputs can be scalars, vectors, or arrays of same size.
            %   Returns interpolated risk value.

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            N = N(:);
            E = E(:);
            alt = alt(:);

            risk = obj.F_interp(E, N, alt);
        end

        function [risks, cumulative] = get_risk_along_path(obj, path, dt)
            %GET_RISK_ALONG_PATH Calculate risk exposure along a trajectory
            %   [risks, cumulative] = get_risk_along_path(path)
            %   [risks, cumulative] = get_risk_along_path(path, dt)
            %
            %   path: [3 x T] trajectory [N; E; D] in NED frame
            %   dt: time step (default: 1)
            %
            %   Returns:
            %     risks: detection probability at each point
            %     cumulative: integrated risk exposure

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            if nargin < 3
                dt = 1;
            end

            N = path(1, :);
            E = path(2, :);
            alt = -path(3, :);  % Convert NED D to altitude

            risks = obj.get_risk(N, E, alt);
            risks = risks(:)';

            cumulative = cumsum(risks) * dt;
        end

        function [max_risk, avg_risk, integrated] = evaluate_path(obj, path, dt)
            %EVALUATE_PATH Compute risk metrics for a path
            %   [max_risk, avg_risk, integrated] = evaluate_path(path, dt)

            if nargin < 3
                dt = 1;
            end

            [risks, ~] = obj.get_risk_along_path(path, dt);

            max_risk = max(risks);
            avg_risk = mean(risks);
            integrated = sum(risks) * dt;
        end

        function fig = plot_slice(obj, alt, fig_handle)
            %PLOT_SLICE Visualize risk at constant altitude
            %   plot_slice(alt) - Create new figure
            %   plot_slice(alt, fig) - Use existing figure

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            if nargin < 3
                fig = figure('Name', sprintf('Threat Map at %.0fm', alt), ...
                             'Position', [100, 100, 900, 700]);
            else
                fig = figure(fig_handle);
            end

            % Find nearest altitude index
            [~, alt_idx] = min(abs(obj.alt_vec - alt));
            actual_alt = obj.alt_vec(alt_idx);

            % Extract slice
            risk_slice = obj.risk_grid(:, :, alt_idx);

            % Get terrain at this altitude (show as blocked)
            [N_grid, E_grid] = meshgrid(obj.N_vec, obj.E_vec);
            terrain_h = obj.terrain.get_height(N_grid(:), E_grid(:));
            terrain_h = reshape(terrain_h, size(N_grid));

            % Mask underground areas
            underground_mask = actual_alt < terrain_h;

            % Plot
            imagesc(obj.N_vec, obj.E_vec, risk_slice);
            hold on;

            % Overlay terrain blocking
            contour(N_grid, E_grid, terrain_h, [actual_alt actual_alt], ...
                    'k-', 'LineWidth', 2);

            % Show underground areas
            if any(underground_mask(:))
                [C, h] = contourf(N_grid, E_grid, double(underground_mask), [0.5 0.5]);
                set(h, 'FaceColor', [0.4 0.3 0.2], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            end

            % Plot radar positions
            for r = 1:length(obj.radars)
                radar = obj.radars{r};
                plot(radar.position(1), radar.position(2), 'r^', ...
                     'MarkerSize', 15, 'MarkerFaceColor', 'r', 'LineWidth', 2);
                text(radar.position(1), radar.position(2) + 30, radar.name, ...
                     'Color', 'r', 'FontWeight', 'bold');
            end

            colormap(flipud(hot));
            %colorbar('Label', 'Detection Probability');
            clim([0 1]);
            xlabel('North [m]');
            ylabel('East [m]');
            title(sprintf('Threat Map at %.0fm Altitude (%.0fm actual)', alt, actual_alt));
            axis equal;
            set(gca, 'YDir', 'normal');
        end

        function fig = plot_3d(obj, threshold, fig_handle)
            %PLOT_3D 3D visualization of threat volume
            %   plot_3d() - Show isosurface at P_det = 0.5
            %   plot_3d(threshold) - Custom threshold
            %   plot_3d(threshold, fig) - Use existing figure

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            if nargin < 2
                threshold = 0.5;
            end

            if nargin < 3
                fig = figure('Name', 'Threat Volume', 'Position', [100, 100, 1000, 800]);
            else
                fig = figure(fig_handle);
            end

            % Create mesh grids
            [N_grid, E_grid, alt_grid] = meshgrid(obj.N_vec, obj.E_vec, obj.alt_vec);

            % Plot terrain surface
            [N_surf, E_surf] = meshgrid(obj.N_vec, obj.E_vec);
            terrain_h = obj.terrain.get_height(N_surf(:), E_surf(:));
            terrain_h = reshape(terrain_h, size(N_surf));

            surf(N_surf, E_surf, terrain_h, 'EdgeColor', 'none', ...
                 'FaceAlpha', 0.7, 'FaceColor', [0.5 0.4 0.3]);
            hold on;

            % Plot threat isosurface
            if max(obj.risk_grid(:)) >= threshold
                p = patch(isosurface(N_grid, E_grid, alt_grid, obj.risk_grid, threshold));
                set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                isonormals(N_grid, E_grid, alt_grid, obj.risk_grid, p);
            end

            % Plot radar positions
            for r = 1:length(obj.radars)
                radar = obj.radars{r};
                plot3(radar.position(1), radar.position(2), radar.position(3), ...
                      'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'LineWidth', 2);
            end

            xlabel('North [m]');
            ylabel('East [m]');
            zlabel('Altitude [m]');
            title(sprintf('Threat Volume (P_{det} > %.2f)', threshold));
            view(30, 30);
            axis equal;
            grid on;
            lighting gouraud;
            camlight;
        end

        function fig = plot_path_risk(obj, path, fig_handle)
            %PLOT_PATH_RISK Visualize risk along a flight path

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            if nargin < 3
                fig = figure('Name', 'Path Risk Analysis', 'Position', [100, 100, 1000, 600]);
            else
                fig = figure(fig_handle);
            end

            N = path(1, :);
            E = path(2, :);
            alt = -path(3, :);

            risks = obj.get_risk(N, E, alt);
            t = 1:length(risks);

            subplot(2, 1, 1);
            plot(t, risks, 'r-', 'LineWidth', 2);
            hold on;
            plot(t, 0.5 * ones(size(t)), 'k--', 'LineWidth', 1);
            xlabel('Path Index');
            ylabel('Detection Probability');
            title('Risk Along Path');
            ylim([0 1]);
            grid on;
            legend('P_{det}', 'Threshold', 'Location', 'best');

            subplot(2, 1, 2);
            cumulative = cumsum(risks);
            plot(t, cumulative, 'b-', 'LineWidth', 2);
            xlabel('Path Index');
            ylabel('Cumulative Risk');
            title('Integrated Risk Exposure');
            grid on;
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
