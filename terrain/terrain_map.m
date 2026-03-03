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
        height_query_mode % 'lut' (default) or 'mesh'
    end

    properties (Access = private)
        F_interp        % Interpolant object for fast queries
        mesh_source     % Optional terrain_mesh handle for high-precision queries

        % Precomputed LUT (uniform grid, C++-friendly static 2D array layout)
        lut_Z           % [nE x nN] heightmap
        lut_N0          % LUT origin N
        lut_E0          % LUT origin E
        lut_dN          % LUT spacing in N
        lut_dE          % LUT spacing in E
        lut_inv_dN      % 1/lut_dN (precomputed for fast mapping)
        lut_inv_dE      % 1/lut_dE (precomputed for fast mapping)
        lut_nN          % Number of LUT columns (N dimension)
        lut_nE          % Number of LUT rows (E dimension)
        lut_num_rows    % Cached row count (= lut_nE) for manual linear indexing
        lut_bounds      % [N_min N_max E_min E_max]
        lut_ready       % True when LUT is initialized
        height_mode_code % 1=luts, 2=mesh (fast dispatch in get_height)
    end

    methods
        function obj = terrain_map(terrain_data, interp_method)
            %TERRAIN_MAP Constructor
            %   tm = terrain_map(terrain_data) - Load from terrain_generator output
            %   tm = terrain_map(terrain_data, 'cubic') - Specify interpolation

            if nargin < 2
                interp_method = 'linear';
            end

            % Ensure vectors are row vectors for consistency
            obj.N_vec = terrain_data.N_vec(:)';
            obj.E_vec = terrain_data.E_vec(:)';
            obj.Z = terrain_data.Z;
            obj.bounds = terrain_data.bounds;
            obj.resolution = terrain_data.resolution;
            obj.interp_method = interp_method;

            if isfield(terrain_data, 'type')
                obj.type = terrain_data.type;
            else
                obj.type = 'unknown';
            end

            % Ensure vectors are monotonically increasing
            if obj.N_vec(1) > obj.N_vec(end)
                obj.N_vec = flip(obj.N_vec);
                obj.Z = flip(obj.Z, 2);  % Flip along N dimension (columns)
            end
            if obj.E_vec(1) > obj.E_vec(end)
                obj.E_vec = flip(obj.E_vec);
                obj.Z = flip(obj.Z, 1);  % Flip along E dimension (rows)
            end

            % griddedInterpolant expects Z of size [length(E_vec), length(N_vec)]
            % where rows correspond to E and columns to N
            expected_size = [length(obj.E_vec), length(obj.N_vec)];
            actual_size = size(obj.Z);

            % Check if Z needs to be transposed
            if isequal(actual_size, flip(expected_size))
                obj.Z = obj.Z';
            elseif ~isequal(actual_size, expected_size)
                error('terrain_map:SizeMismatch', ...
                    'Z matrix size [%d, %d] does not match grid vectors [%d, %d]', ...
                    actual_size(1), actual_size(2), expected_size(1), expected_size(2));
            end

            % Create interpolant for fast queries
            % griddedInterpolant({E_vec, N_vec}, Z) expects:
            %   - E_vec and N_vec as row or column vectors (monotonically increasing)
            %   - Z as size [length(E_vec), length(N_vec)]
            %   - Query with F_interp(E, N) returns Z value
            try
                obj.F_interp = griddedInterpolant({obj.E_vec, obj.N_vec}, obj.Z, ...
                    interp_method, 'nearest');
            catch ME
                % Fallback: create a wrapper using interp2
                warning('terrain_map:InterpolantFallback', ...
                    'griddedInterpolant failed (%s), using interp2 fallback.', ME.message);
                obj.F_interp = [];  % Will use interp2 in get_height
            end

            obj.mesh_source = [];
            obj.height_query_mode = 'lut';
            obj.height_mode_code = 1;
            obj.lut_ready = false;
            obj.initialize_height_lut();
        end

        function h = get_height(obj, N, E, mode)
            %GET_HEIGHT Get terrain height at specified coordinates
            %   h = get_height(N, E) - Returns height(s) at position(s)
            %   h = get_height(N, E, mode) where mode is 'lut' or 'mesh'
            %
            %   N, E can be scalars, vectors, or matrices of same size.
            %   Returns NaN for points outside terrain bounds.

            if ~isequal(size(N), size(E))
                error('terrain_map:InputSizeMismatch', 'N and E must have the same size.');
            end

            if nargin < 4 || isempty(mode)
                mode_code = obj.height_mode_code;
            else
                if isstring(mode)
                    mode = char(mode);
                end
                if ~ischar(mode)
                    error('terrain_map:UnknownHeightMode', ...
                        'Unknown height mode type. Use ''lut'' or ''mesh''.');
                end
                mode_l = lower(mode);
                if strcmp(mode_l, 'lut') || strcmp(mode_l, 'grid') || strcmp(mode_l, 'heightmap')
                    mode_code = 1;
                elseif strcmp(mode_l, 'mesh') || strcmp(mode_l, 'ray') || strcmp(mode_l, 'raycast')
                    mode_code = 2;
                else
                    error('terrain_map:UnknownHeightMode', ...
                        'Unknown height mode "%s". Use ''lut'' or ''mesh''.', mode);
                end
            end

            if mode_code == 1
                h = obj.get_height_lut_impl(N, E);
            else
                h = obj.get_height_mesh_impl(N, E);
            end
        end

        function set_height_query_mode(obj, mode)
            %SET_HEIGHT_QUERY_MODE Set default mode used by get_height(N,E).
            if isstring(mode)
                mode = char(mode);
            end
            mode = lower(mode);
            if strcmp(mode, 'lut')
                obj.height_query_mode = mode;
                obj.height_mode_code = 1;
            elseif strcmp(mode, 'mesh')
                obj.height_query_mode = mode;
                obj.height_mode_code = 2;
            else
                error('terrain_map:UnknownHeightMode', ...
                    'Unknown mode "%s". Use ''lut'' or ''mesh''.', mode);
            end
        end

        function attach_mesh(obj, mesh_obj)
            %ATTACH_MESH Attach terrain_mesh for high-precision mode.
            obj.mesh_source = mesh_obj;
        end

        function initialize_height_lut(obj, mission_bounds, lut_resolution)
            %INITIALIZE_HEIGHT_LUT Precompute uniform LUT for O(1) bilinear queries.
            %   initialize_height_lut() uses full map bounds and native resolution.
            %   initialize_height_lut([Nmin Nmax Emin Emax], [dN dE]) uses custom area/grid.

            if nargin < 2 || isempty(mission_bounds)
                mission_bounds = obj.bounds(1:4);
            end
            if nargin < 3 || isempty(lut_resolution)
                dN = mean(diff(obj.N_vec));
                dE = mean(diff(obj.E_vec));
            else
                if isscalar(lut_resolution)
                    dN = lut_resolution;
                    dE = lut_resolution;
                else
                    dN = lut_resolution(1);
                    dE = lut_resolution(2);
                end
            end

            if dN <= 0 || dE <= 0
                error('terrain_map:InvalidLUTResolution', 'LUT resolution must be positive.');
            end

            N_min = max(obj.bounds(1), mission_bounds(1));
            N_max = min(obj.bounds(2), mission_bounds(2));
            E_min = max(obj.bounds(3), mission_bounds(3));
            E_max = min(obj.bounds(4), mission_bounds(4));
            if N_min >= N_max || E_min >= E_max
                error('terrain_map:InvalidLUTBounds', ...
                    'Mission bounds must overlap terrain bounds with non-zero area.');
            end

            N_lut = N_min:dN:N_max;
            E_lut = E_min:dE:E_max;
            if N_lut(end) < N_max
                N_lut = [N_lut, N_max]; %#ok<AGROW>
            end
            if E_lut(end) < E_max
                E_lut = [E_lut, E_max]; %#ok<AGROW>
            end

            [E_grid, N_grid] = ndgrid(E_lut, N_lut);
            if ~isempty(obj.F_interp)
                Z_lut = obj.F_interp(E_grid, N_grid);
            else
                Z_lut = interp2(obj.N_vec, obj.E_vec, obj.Z, N_grid, E_grid, obj.interp_method);
            end

            obj.lut_Z = Z_lut;
            obj.lut_N0 = N_lut(1);
            obj.lut_E0 = E_lut(1);
            obj.lut_dN = N_lut(2) - N_lut(1);
            obj.lut_dE = E_lut(2) - E_lut(1);
            obj.lut_inv_dN = 1 / obj.lut_dN;
            obj.lut_inv_dE = 1 / obj.lut_dE;
            obj.lut_nN = numel(N_lut);
            obj.lut_nE = numel(E_lut);
            obj.lut_num_rows = obj.lut_nE;
            obj.lut_bounds = [N_lut(1), N_lut(end), E_lut(1), E_lut(end)];
            obj.lut_ready = true;
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
            surf(N_grid, E_grid, obj.Z, 'EdgeColor', 'none', 'FaceAlpha', 1.0); % Opaque for proper lighting
            hold on;

            % Formatting
            colormap(terrain_colormap());

            % Add colorbar (safely)
            cb = colorbar;
            cb.Label.String = 'Elevation [m]';

            xlabel('North [m]');
            ylabel('East [m]');
            zlabel('Height [m]');
            title(sprintf('Terrain: %s', obj.type), 'FontSize', 12, 'FontWeight', 'bold');
            axis equal;
            grid on;
            view(30, 45);

            % Add contour lines at base
            contour3(N_grid, E_grid, obj.Z, 15, 'k', 'LineWidth', 0.5);

            % Add lighting for better 3D depth perception
            light('Position', [-1 -1 1], 'Style', 'infinite');
            lighting gouraud;
            material dull; % Non-shiny surface looks more like terrain
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

            % Plot path with shadow/curtain
            plot3(N, E, alt, 'r-', 'LineWidth', 2.5);

            % Mark start and end
            plot3(N(1), E(1), alt(1), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
            plot3(N(end), E(end), alt(end), 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

            % Plot shadow on terrain (projected path)
            terrain_h = obj.get_height(N, E);
            plot3(N, E, terrain_h + 1, 'k--', 'LineWidth', 1.5); % +1m to avoid Z-fighting

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

            % Fill contours
            contourf(N_grid, E_grid, obj.Z, 20, 'LineStyle', 'none');
            hold on;
            % Add black contour lines for definition
            contour(N_grid, E_grid, obj.Z, 20, 'k', 'LineWidth', 0.5);

            colormap(terrain_colormap());
            cb = colorbar;
            cb.Label.String = 'Elevation [m]';

            xlabel('North [m]');
            ylabel('East [m]');
            title(sprintf('Terrain Contour: %s', obj.type));
            axis equal;
            grid on;
        end
    end

    methods (Access = private)
        function h = get_height_lut_impl(obj, N, E)
            % Fully vectorized bilinear interpolation on precomputed LUT.
            if ~obj.lut_ready
                obj.initialize_height_lut();
            end

            sz = size(N);
            Nq = N(:);
            Eq = E(:);
            h = nan(size(Nq));

            % Single vectorized bounds check before interpolation math.
            valid = Nq >= obj.lut_bounds(1) & Nq <= obj.lut_bounds(2) & ...
                    Eq >= obj.lut_bounds(3) & Eq <= obj.lut_bounds(4);
            if ~any(valid)
                h = reshape(h, sz);
                return;
            end

            if obj.lut_nN < 2 || obj.lut_nE < 2
                % Degenerate LUT fallback.
                if ~isempty(obj.F_interp)
                    h(valid) = obj.F_interp(Eq(valid), Nq(valid));
                else
                    h(valid) = interp2(obj.N_vec, obj.E_vec, obj.Z, ...
                        Nq(valid), Eq(valid), obj.interp_method);
                end
                h = reshape(h, sz);
                return;
            end

            use_all = all(valid);
            if use_all
                Nv = Nq;
                Ev = Eq;
            else
                % Work only on valid points to minimize temporary allocations.
                Nv = Nq(valid);
                Ev = Eq(valid);
            end

            % Use reciprocal spacing to replace divisions with multiplies.
            u = (Nv - obj.lut_N0) .* obj.lut_inv_dN + 1;
            v = (Ev - obj.lut_E0) .* obj.lut_inv_dE + 1;

            j0 = floor(u); % column index (N dimension)
            i0 = floor(v); % row index (E dimension)
            tx = u - j0;
            ty = v - i0;

            % Clamp once for boundary-safe bilinear access.
            j0 = max(1, min(obj.lut_nN - 1, j0));
            i0 = max(1, min(obj.lut_nE - 1, i0));
            tx = max(0, min(1, tx));
            ty = max(0, min(1, ty));

            % Manual linear indexing (faster than sub2ind in this hot path):
            % idx = (col-1)*numRows + row
            base = (j0 - 1) * obj.lut_num_rows + i0;
            idx11 = base;
            idx12 = base + obj.lut_num_rows;
            idx21 = base + 1;
            idx22 = base + obj.lut_num_rows + 1;

            z11 = obj.lut_Z(idx11);
            z12 = obj.lut_Z(idx12);
            z21 = obj.lut_Z(idx21);
            z22 = obj.lut_Z(idx22);

            % Bilinear interpolation (algebraically reduced form).
            hv = z11 + tx .* (z12 - z11) + ty .* (z21 - z11) + ...
                tx .* ty .* (z11 - z12 - z21 + z22);

            if use_all
                h = hv;
            else
                h(valid) = hv;
            end

            h = reshape(h, sz);
        end

        function h = get_height_mesh_impl(obj, N, E)
            % Delegate to terrain_mesh precision query path.
            if isempty(obj.mesh_source)
                h = obj.get_height_lut_impl(N, E);
                return;
            end

            try
                h = obj.mesh_source.get_height(N, E, 'raycast');
            catch
                h = obj.mesh_source.get_height(N, E);
            end

            h = reshape(h, size(N));
        end
    end
end

function cmap = terrain_colormap()
% Custom colormap for terrain visualization
% High contrast: Blue(water) -> Green(low info) -> Brown(mid) -> Gray(high) -> White(peak)

% Control points (Normalized elevation 0.0 to 1.0)
% We assume 0 is lowest and 1 is highest in the current view

% Colors: [R G B]
colors = [
    0.2 0.4 0.8;   % Deep Blue (Water/Base)
    0.1 0.6 0.2;   % Dark Green (Low vegetation)
    0.4 0.8 0.4;   % Light Green (Grass)
    0.8 0.7 0.4;   % Sand/Dirt
    0.5 0.3 0.1;   % Dark Brown (Earth)
    0.5 0.5 0.5;   % Gray (Rock)
    0.9 0.9 0.9;   % White (Snow/Peak)
    ];

% Create strictly interpolated map
n_colors = 256;
x = linspace(0, 1, size(colors, 1));
xi = linspace(0, 1, n_colors);
cmap = interp1(x, colors, xi, 'pchip'); % 'pchip' for smooth but vibrant transition
end
