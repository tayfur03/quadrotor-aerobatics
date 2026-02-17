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

        function compute_map(obj, method, opts)
            %COMPUTE_MAP Calculate risk at all grid points
            %   compute_map()        - Use 'max' combination
            %   compute_map('max')   - Max detection probability across radars
            %   compute_map('sum')   - Sum (can exceed 1, useful for cost)
            %   compute_map('prob')  - Probabilistic: 1 - prod(1-P_i)
            %   compute_map('binary')- Deterministic visibility: 1 if visible
            %   compute_map(method, opts)
            %     opts.use_parallel - true/false, use parfor over altitude slices
            %     opts.show_progress - true/false progress prints (serial mode)

            if nargin < 2
                method = 'max';
            end
            if nargin < 3
                opts = struct();
            end

            use_parallel = get_param(opts, 'use_parallel', false);
            show_progress = get_param(opts, 'show_progress', true);
            if use_parallel && ~can_use_parfor()
                warning('threat_map:parforUnavailable', ...
                    'Parallel Computing Toolbox not available. Falling back to serial.');
                use_parallel = false;
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
            use_binary = strcmpi(method, 'binary');
            n_rad_total = length(obj.radars);
            n_rad = 0;

            % Initialize grid
            obj.risk_grid = zeros(n_E, n_N, n_alt);

            % Cache radar parameters for faster binary visibility checks
            if use_binary
                n_rad = n_rad_total;
                radar_enabled = false(1, n_rad);
                radar_pos = zeros(3, n_rad);
                radar_R2 = zeros(1, n_rad);
                radar_az_limits = cell(1, n_rad);
                radar_el_limits = cell(1, n_rad);
                for r = 1:n_rad
                    rd = obj.radars{r};
                    radar_enabled(r) = rd.enabled;
                    radar_pos(:, r) = rd.position(:);
                    radar_R2(r) = rd.R_max^2;
                    radar_az_limits{r} = rd.azimuth_limits;
                    radar_el_limits{r} = rd.elevation_limits;
                end
            end

            % Precompute terrain surface once (same for all altitudes)
            [N_grid_2d, E_grid_2d] = meshgrid(obj.N_vec, obj.E_vec);
            terrain_h_2d = obj.terrain.get_height(N_grid_2d(:), E_grid_2d(:));
            terrain_h_2d = reshape(terrain_h_2d, [n_E, n_N]);

            % Local broadcast variables for speed/parfor compatibility
            N_vec_local = obj.N_vec;
            E_vec_local = obj.E_vec;
            alt_vec_local = obj.alt_vec;
            radars_local = obj.radars;
            los_local = obj.los;
            rcs_local = obj.rcs;

            if use_parallel
                fprintf('Parallel threat-map computation enabled (parfor).\n');
                risk_local = zeros(n_E, n_N, n_alt);

                parfor k = 1:n_alt
                    alt = alt_vec_local(k);
                    risk_slice = zeros(n_E, n_N);

                    for j = 1:n_N
                        N = N_vec_local(j);

                        for i = 1:n_E
                            E = E_vec_local(i);
                            terrain_h = terrain_h_2d(i, j);
                            if alt < terrain_h
                                risk_slice(i, j) = 1.0;
                                continue;
                            end

                            target_pos = [N; E; alt];

                            if use_binary
                                visible_any = false;
                                for r = 1:n_rad
                                    if ~radar_enabled(r)
                                        continue;
                                    end

                                    rel = target_pos - radar_pos(:, r);
                                    range2 = rel(1)^2 + rel(2)^2 + rel(3)^2;
                                    if range2 > radar_R2(r)
                                        continue;
                                    end

                                    if ~is_in_coverage_fast(rel, radar_az_limits{r}, radar_el_limits{r})
                                        continue;
                                    end

                                    if ~isempty(los_local)
                                        visible_any = los_local.has_los(radar_pos(:, r), target_pos);
                                    else
                                        visible_any = true;
                                    end

                                    if visible_any
                                        break;
                                    end
                                end
                                risk_slice(i, j) = double(visible_any);
                            else
                                P_det_list = zeros(1, n_rad_total);
                                for r = 1:n_rad_total
                                    P_det_list(r) = radars_local{r}.get_detection_probability(...
                                        target_pos, rcs_local, los_local);
                                end

                                switch lower(method)
                                    case 'max'
                                        risk_slice(i, j) = max(P_det_list);
                                    case 'sum'
                                        risk_slice(i, j) = sum(P_det_list);
                                    case 'prob'
                                        risk_slice(i, j) = 1 - prod(1 - P_det_list);
                                    otherwise
                                        error('Unknown threat map method: %s', method);
                                end
                            end
                        end
                    end

                    risk_local(:, :, k) = risk_slice;
                end

                obj.risk_grid = risk_local;
            else
                progress_interval = max(1, floor(total_cells / 20));
                cell_count = 0;

                for k = 1:n_alt
                    alt = alt_vec_local(k);

                    for j = 1:n_N
                        N = N_vec_local(j);

                        for i = 1:n_E
                            E = E_vec_local(i);

                            terrain_h = terrain_h_2d(i, j);
                            if alt < terrain_h
                                obj.risk_grid(i, j, k) = 1.0;
                                continue;
                            end

                            target_pos = [N; E; alt];

                            if use_binary
                                visible_any = false;
                                for r = 1:n_rad
                                    if ~radar_enabled(r)
                                        continue;
                                    end

                                    rel = target_pos - radar_pos(:, r);
                                    range2 = rel(1)^2 + rel(2)^2 + rel(3)^2;
                                    if range2 > radar_R2(r)
                                        continue;
                                    end

                                    if ~is_in_coverage_fast(rel, radar_az_limits{r}, radar_el_limits{r})
                                        continue;
                                    end

                                    if ~isempty(los_local)
                                        visible_any = los_local.has_los(radar_pos(:, r), target_pos);
                                    else
                                        visible_any = true;
                                    end

                                    if visible_any
                                        break;
                                    end
                                end
                                obj.risk_grid(i, j, k) = double(visible_any);
                            else
                                P_det_list = zeros(1, n_rad_total);
                                for r = 1:n_rad_total
                                    P_det_list(r) = radars_local{r}.get_detection_probability(...
                                        target_pos, rcs_local, los_local);
                                end

                                switch lower(method)
                                    case 'max'
                                        obj.risk_grid(i, j, k) = max(P_det_list);
                                    case 'sum'
                                        obj.risk_grid(i, j, k) = sum(P_det_list);
                                    case 'prob'
                                        obj.risk_grid(i, j, k) = 1 - prod(1 - P_det_list);
                                    otherwise
                                        error('Unknown threat map method: %s', method);
                                end
                            end

                            cell_count = cell_count + 1;
                            if show_progress && mod(cell_count, progress_interval) == 0
                                fprintf('  Progress: %d%%\n', round(100*cell_count/total_cells));
                            end
                        end
                    end
                end
            end

            % Create interpolant for fast queries
            interp_method = 'linear';
            if strcmpi(method, 'binary')
                interp_method = 'nearest';
            end
            obj.F_interp = griddedInterpolant({obj.E_vec, obj.N_vec, obj.alt_vec}, ...
                                               obj.risk_grid, interp_method, 'nearest');

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

        function [class_slice, alt_used] = get_binary_horizontal_slice(obj, alt_query, visibility_threshold, terrain_obj, terrain_h_grid)
            %GET_BINARY_HORIZONTAL_SLICE Return 3-class horizontal binary slice.
            %   class values:
            %     1 -> below terrain
            %     2 -> hidden/safe
            %     3 -> visible

            if nargin < 3 || isempty(visibility_threshold)
                visibility_threshold = 0.5;
            end
            if nargin < 4 || isempty(terrain_obj)
                terrain_obj = obj.terrain;
            end
            if nargin < 5
                terrain_h_grid = [];
            end

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            [~, alt_idx] = min(abs(obj.alt_vec - alt_query));
            alt_used = obj.alt_vec(alt_idx);
            risk_slice = obj.risk_grid(:, :, alt_idx);

            class_slice = 2 * ones(size(risk_slice));  % hidden/safe
            class_slice(risk_slice >= visibility_threshold) = 3; % visible

            if isempty(terrain_h_grid)
                [N_grid, E_grid] = meshgrid(obj.N_vec, obj.E_vec);
                terrain_h_grid = terrain_obj.get_height(N_grid(:), E_grid(:));
                terrain_h_grid = reshape(terrain_h_grid, size(N_grid));
            end
            class_slice(alt_used < terrain_h_grid) = 1; % below terrain
        end

        function [class_grid, s_vec, alt_vec, terrain_line, N_line, E_line] = ...
                get_binary_vertical_slice(obj, anchor_NE, heading_rad, half_length, n_horiz, alt_vec, visibility_threshold, terrain_obj)
            %GET_BINARY_VERTICAL_SLICE Return 3-class vertical slice along heading.
            %   class values:
            %     1 -> below terrain
            %     2 -> hidden/safe
            %     3 -> visible
            %
            % Inputs:
            %   anchor_NE  : [2x1] [N;E] slice center
            %   heading_rad: heading angle in radians (atan2(dE,dN))

            if nargin < 4 || isempty(half_length)
                half_length = 2500;
            end
            if nargin < 5 || isempty(n_horiz)
                n_horiz = 181;
            end
            if nargin < 6 || isempty(alt_vec)
                alt_vec = obj.alt_vec;
            end
            if nargin < 7 || isempty(visibility_threshold)
                visibility_threshold = 0.5;
            end
            if nargin < 8 || isempty(terrain_obj)
                terrain_obj = obj.terrain;
            end

            if ~obj.computed
                error('Threat map not computed. Call compute_map() first.');
            end

            anchor_NE = anchor_NE(:);
            if numel(anchor_NE) ~= 2
                error('anchor_NE must be [2x1] as [N;E].');
            end

            s_vec = linspace(-half_length, half_length, n_horiz);
            N_line = anchor_NE(1) + s_vec * cos(heading_rad);
            E_line = anchor_NE(2) + s_vec * sin(heading_rad);
            terrain_line = terrain_obj.get_height(N_line, E_line);
            terrain_line = reshape(terrain_line, 1, []); % force row [1 x n_horiz]

            [S_grid, A_grid] = meshgrid(s_vec, alt_vec);
            N_query = anchor_NE(1) + S_grid * cos(heading_rad);
            E_query = anchor_NE(2) + S_grid * sin(heading_rad);

            risk = obj.get_risk(N_query(:), E_query(:), A_grid(:));
            risk = reshape(risk, size(S_grid));

            class_grid = 2 * ones(size(risk)); % hidden/safe
            class_grid(risk >= visibility_threshold) = 3; % visible
            % Build below-terrain mask with explicit broadcasting safety.
            below_mask = bsxfun(@lt, alt_vec(:), terrain_line); % [n_alt x n_horiz]
            class_grid(below_mask) = 1; % below terrain
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

function in_coverage = is_in_coverage_fast(rel_pos, az_limits, el_limits)
    % Fast equivalent of radar_site.is_in_coverage for binary map evaluation.
    in_coverage = true;

    range_horiz = hypot(rel_pos(1), rel_pos(2));
    azimuth = atan2d(rel_pos(2), rel_pos(1));
    elevation = atan2d(rel_pos(3), range_horiz);

    if ~isempty(az_limits)
        az_min = az_limits(1);
        az_max = az_limits(2);
        if az_min <= az_max
            in_coverage = in_coverage && (azimuth >= az_min && azimuth <= az_max);
        else
            in_coverage = in_coverage && (azimuth >= az_min || azimuth <= az_max);
        end
    end

    if ~isempty(el_limits)
        in_coverage = in_coverage && ...
            (elevation >= el_limits(1) && elevation <= el_limits(2));
    end
end

function tf = can_use_parfor()
    tf = license('test', 'Distrib_Computing_Toolbox') && ...
         (exist('parfor', 'builtin') > 0 || exist('parfor', 'file') > 0);
end
