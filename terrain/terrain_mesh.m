classdef terrain_mesh < handle
%TERRAIN_MESH Triangle mesh representation of terrain for ray casting
%
% Converts a DEM grid into a triangle mesh and provides efficient
% ray-terrain intersection using the Möller-Trumbore algorithm with
% spatial acceleration via grid-based indexing.
%
% The mesh is created by converting each grid cell into two triangles.
% Spatial indexing groups triangles into grid cells for fast ray traversal.
%
% Properties:
%   vertices     - [N x 3] vertex positions [N, E, altitude]
%   triangles    - [M x 3] triangle vertex indices
%   tri_normals  - [M x 3] triangle face normals
%   bounds       - [N_min, N_max, E_min, E_max, alt_min, alt_max]
%   grid_size    - Size of spatial index grid [rows, cols]
%   cell_size    - Spatial index cell size [N_size, E_size]
%   grid_index   - Cell array mapping grid cells to triangle indices
%
% Methods:
%   ray_intersect(origin, direction, max_range)  - Find ray-terrain intersection
%   ray_intersect_all(origin, direction)         - Find all intersections
%   get_height(N, E)                             - Get terrain height at point
%   plot()                                       - Visualize mesh
%
% Example:
%   % Create mesh from DEM data
%   terrain_data = dem_loader('elevation.tif');
%   mesh = terrain_mesh(terrain_data);
%
%   % Ray casting
%   origin = [0; 0; 500];        % Radar position [N; E; altitude]
%   direction = [1; 0; -0.1];    % Ray direction
%   [hit, point, dist] = mesh.ray_intersect(origin, direction, 10000);
%
% See also: moller_trumbore, dem_loader, radar_ray_caster

    properties
        vertices        % [N x 3] vertex positions [N, E, altitude]
        triangles       % [M x 3] triangle vertex indices (1-based)
        tri_normals     % [M x 3] triangle face normals
        bounds          % [N_min, N_max, E_min, E_max, alt_min, alt_max]
        resolution      % Original DEM resolution
        n_rows          % Number of rows in original grid
        n_cols          % Number of columns in original grid
        grid_index      % Cell array for spatial acceleration
    end

    properties (Access = private)
        
        grid_size       % [rows, cols] of spatial index
        cell_size       % [N_size, E_size] of each index cell
        N_vec           % Original N coordinate vector
        E_vec           % Original E coordinate vector
        F_height        % Vectorized height interpolant for fast queries
    end

    methods
        function obj = terrain_mesh(terrain_data, index_cell_size)
            %TERRAIN_MESH Constructor - create mesh from DEM data
            %   mesh = terrain_mesh(terrain_data)
            %   mesh = terrain_mesh(terrain_data, index_cell_size)
            %
            %   terrain_data: struct with N_vec, E_vec, Z fields
            %   index_cell_size: (optional) size of spatial index cells [m]

            if nargin < 2
                % Default: 10x the DEM resolution for index cells
                index_cell_size = terrain_data.resolution * 10;
            end

            obj.N_vec = terrain_data.N_vec(:)';
            obj.E_vec = terrain_data.E_vec(:)';
            obj.resolution = terrain_data.resolution;
            obj.n_rows = length(obj.E_vec);
            obj.n_cols = length(obj.N_vec);

            % Get Z matrix and check orientation
            Z = terrain_data.Z;
            expected_size = [obj.n_rows, obj.n_cols];  % [length(E_vec), length(N_vec)]
            actual_size = size(Z);
            
            % Check if Z needs to be transposed (same logic as terrain_map)
            if isequal(actual_size, flip(expected_size))
                Z = Z';
                fprintf('terrain_mesh: Z matrix transposed to [%d, %d]\n', size(Z, 1), size(Z, 2));
            elseif ~isequal(actual_size, expected_size)
                warning('terrain_mesh:SizeMismatch', ...
                    'Z matrix size [%d, %d] does not match grid vectors [%d, %d]. Attempting to continue...', ...
                    actual_size(1), actual_size(2), expected_size(1), expected_size(2));
            end

            % Fast vectorized height lookup for batch segment checks.
            try
                obj.F_height = griddedInterpolant({obj.E_vec, obj.N_vec}, Z, 'linear', 'nearest');
            catch
                obj.F_height = [];
            end

            % Create triangle mesh from grid
            obj.triangulate_grid(Z);

            % Build spatial acceleration index
            obj.build_spatial_index(index_cell_size);

            fprintf('Terrain mesh created: %d vertices, %d triangles\n', ...
                size(obj.vertices, 1), size(obj.triangles, 1));
        end

        function [hit, intersection_point, distance, tri_idx, normal] = ray_intersect(obj, origin, direction, max_range)
            %RAY_INTERSECT Find closest ray-terrain intersection
            %   [hit, point, dist, tri, normal] = ray_intersect(origin, dir, max_range)
            %
            %   Inputs:
            %       origin    - [3x1] ray origin [N; E; altitude]
            %       direction - [3x1] ray direction (will be normalized)
            %       max_range - Maximum range to search [m]
            %
            %   Outputs:
            %       hit       - true if intersection found
            %       intersection_point - [3x1] intersection position
            %       distance  - Distance from origin to intersection
            %       tri_idx   - Index of intersected triangle
            %       normal    - [3x1] surface normal at intersection

            origin = origin(:);
            direction = direction(:);

            % Normalize direction
            dir_len = norm(direction);
            if dir_len < 1e-10
                hit = false;
                intersection_point = origin;
                distance = 0;
                tri_idx = 0;
                normal = [0; 0; 1];
                return;
            end
            direction = direction / dir_len;

            % Get candidate triangles using spatial index
            candidate_tris = obj.get_candidate_triangles(origin, direction, max_range);

            best_t = max_range;
            hit = false;
            tri_idx = 0;

            if ~isempty(candidate_tris)
                % Vectorized triangle batch against one ray.
                tri_rows = obj.triangles(candidate_tris, :); % [K x 3]
                v0 = obj.vertices(tri_rows(:, 1), :)';       % [3 x K]
                v1 = obj.vertices(tri_rows(:, 2), :)';
                v2 = obj.vertices(tri_rows(:, 3), :)';

                k = numel(candidate_tris);
                origin_batch = origin(:, ones(1, k));
                direction_batch = direction(:, ones(1, k));

                [tri_hit, t_batch, ~, ~, ~] = moller_trumbore( ...
                    origin_batch, direction_batch, v0, v1, v2);

                valid = tri_hit & isfinite(t_batch) & (t_batch > 0) & (t_batch < max_range);
                if any(valid)
                    t_valid = t_batch(valid);
                    cand_valid = candidate_tris(valid);
                    [best_t, best_loc] = min(t_valid);
                    hit = true;
                    tri_idx = cand_valid(best_loc);
                end
            end

            if hit
                intersection_point = origin + best_t * direction;
                distance = best_t;
                normal = obj.tri_normals(tri_idx, :)';
            else
                intersection_point = origin + max_range * direction;
                distance = max_range;
                normal = [0; 0; 1];
            end
        end

        function [hits, points, distances, tri_indices] = ray_intersect_all(obj, origin, direction, max_range)
            %RAY_INTERSECT_ALL Find all ray-terrain intersections
            %   [hits, points, dists, tris] = ray_intersect_all(origin, dir, max_range)
            %
            %   Returns arrays of all intersections along the ray (sorted by distance)

            origin = origin(:);
            direction = direction(:) / norm(direction(:));

            candidate_tris = obj.get_candidate_triangles(origin, direction, max_range);

            all_t = [];
            all_tri = [];

            if ~isempty(candidate_tris)
                % Vectorized triangle batch against one ray.
                tri_rows = obj.triangles(candidate_tris, :); % [K x 3]
                v0 = obj.vertices(tri_rows(:, 1), :)';       % [3 x K]
                v1 = obj.vertices(tri_rows(:, 2), :)';
                v2 = obj.vertices(tri_rows(:, 3), :)';

                k = numel(candidate_tris);
                origin_batch = origin(:, ones(1, k));
                direction_batch = direction(:, ones(1, k));

                [tri_hit, t_batch, ~, ~, ~] = moller_trumbore( ...
                    origin_batch, direction_batch, v0, v1, v2);

                valid = tri_hit & isfinite(t_batch) & (t_batch > 0) & (t_batch < max_range);
                all_t = t_batch(valid);
                all_tri = candidate_tris(valid);
            end

            % Sort by distance
            [distances, sort_idx] = sort(all_t);
            tri_indices = all_tri(sort_idx);

            hits = length(distances);
            if hits > 0
                points = origin + direction .* distances;
            else
                points = zeros(3, 0);
            end
        end

        function h = get_height(obj, N, E)
            %GET_HEIGHT Get terrain height at specified coordinates
            %   h = get_height(N, E)
            %
            %   Uses vectorized interpolation (fast path), with ray-cast fallback.
            %   Returns NaN for points outside terrain bounds.

            N = N(:);
            E = E(:);
            h = nan(size(N));

            in_bounds = N >= obj.bounds(1) & N <= obj.bounds(2) & ...
                        E >= obj.bounds(3) & E <= obj.bounds(4);

            if ~any(in_bounds)
                return;
            end

            if ~isempty(obj.F_height)
                h(in_bounds) = obj.F_height(E(in_bounds), N(in_bounds));
                return;
            end

            idx = find(in_bounds);
            for k = 1:numel(idx)
                i = idx(k);
                origin = [N(i); E(i); obj.bounds(6) + 1000];
                direction = [0; 0; -1];
                [hit, point, ~, ~, ~] = obj.ray_intersect(origin, direction, 2000);
                if hit
                    h(i) = point(3);
                end
            end
        end

        function fig = plot(obj, show_wireframe)
            %PLOT Visualize terrain mesh
            %   plot() - Solid surface
            %   plot(true) - Show wireframe

            if nargin < 2
                show_wireframe = false;
            end

            fig = figure('Name', 'Terrain Mesh', 'Position', [100, 100, 1000, 800]);

            if show_wireframe
                trisurf(obj.triangles, obj.vertices(:,1), obj.vertices(:,2), obj.vertices(:,3), ...
                    obj.vertices(:,3), 'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 0.8);
            else
                trisurf(obj.triangles, obj.vertices(:,1), obj.vertices(:,2), obj.vertices(:,3), ...
                    obj.vertices(:,3), 'EdgeColor', 'none');
            end

            colormap(terrain_colormap());
            colorbar;
            xlabel('North [m]');
            ylabel('East [m]');
            zlabel('Altitude [m]');
            title(sprintf('Terrain Mesh (%d triangles)', size(obj.triangles, 1)));
            axis equal;
            view(30, 45);
            light('Position', [1 1 1]);
            lighting gouraud;
        end

        function fig = plot_with_ray(obj, origin, direction, max_range)
            %PLOT_WITH_RAY Visualize mesh with ray and intersection
            fig = obj.plot(false);
            hold on;

            [hit, point, ~, ~, ~] = obj.ray_intersect(origin, direction, max_range);

            % Plot ray
            direction = direction(:) / norm(direction(:));
            if hit
                end_point = point;
            else
                end_point = origin(:) + direction * max_range;
            end

            plot3([origin(1), end_point(1)], [origin(2), end_point(2)], ...
                  [origin(3), end_point(3)], 'r-', 'LineWidth', 2);

            % Plot origin
            plot3(origin(1), origin(2), origin(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

            % Plot intersection
            if hit
                plot3(point(1), point(2), point(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            end

            legend('Terrain', 'Ray', 'Origin', 'Intersection');
        end
    end

    methods (Access = private)
        function triangulate_grid(obj, Z)
            %TRIANGULATE_GRID Convert DEM grid to triangle mesh
            %
            % Each grid cell with corners A-B-C-D becomes two triangles:
            %   A---B        Triangle 1: A-B-C
            %   |   |        Triangle 2: A-C-D
            %   D---C

            n_E = obj.n_rows;  % rows in Z matrix (E direction)
            n_N = obj.n_cols;  % cols in Z matrix (N direction)

            % Create vertices: [N, E, altitude]
            n_vertices = n_E * n_N;
            obj.vertices = zeros(n_vertices, 3);

            idx = 1;
            for i = 1:n_E
                for j = 1:n_N
                    obj.vertices(idx, :) = [obj.N_vec(j), obj.E_vec(i), Z(i, j)];
                    idx = idx + 1;
                end
            end

            % Create triangles (2 per grid cell)
            n_triangles = 2 * (n_E - 1) * (n_N - 1);
            obj.triangles = zeros(n_triangles, 3);
            obj.tri_normals = zeros(n_triangles, 3);

            tri_idx = 1;
            fprintf('Starting triangulation...\n');
            for i = 1:(n_E - 1)
                if mod(i, 50) == 0
                    fprintf('Triangulating row %d of %d (%.0f%%)\n', i, n_E-1, (i/(n_E-1))*100);
                end
                for j = 1:(n_N - 1)
                    % Vertex indices for corners A, B, C, D
                    % A = (i, j), B = (i, j+1), C = (i+1, j+1), D = (i+1, j)
                    A = (i - 1) * n_N + j;
                    B = (i - 1) * n_N + (j + 1);
                    C = i * n_N + (j + 1);
                    D = i * n_N + j;

                    % Triangle 1: A-B-C
                    obj.triangles(tri_idx, :) = [A, B, C];
                    obj.tri_normals(tri_idx, :) = obj.compute_normal(A, B, C);
                    tri_idx = tri_idx + 1;

                    % Triangle 2: A-C-D
                    obj.triangles(tri_idx, :) = [A, C, D];
                    obj.tri_normals(tri_idx, :) = obj.compute_normal(A, C, D);
                    tri_idx = tri_idx + 1;
                end
            end
            fprintf('Triangulation complete.\n');

            % Compute bounds
            obj.bounds = [min(obj.vertices(:,1)), max(obj.vertices(:,1)), ...
                          min(obj.vertices(:,2)), max(obj.vertices(:,2)), ...
                          min(obj.vertices(:,3)), max(obj.vertices(:,3))];
        end

        function normal = compute_normal(obj, i1, i2, i3)
            %COMPUTE_NORMAL Calculate face normal for triangle
            v0 = obj.vertices(i1, :);
            v1 = obj.vertices(i2, :);
            v2 = obj.vertices(i3, :);

            e1 = v1 - v0;
            e2 = v2 - v0;
            n = cross(e1, e2);
            normal = n / norm(n);
        end

        function build_spatial_index(obj, cell_size)
            %BUILD_SPATIAL_INDEX Create grid-based spatial index for acceleration

            obj.cell_size = [cell_size, cell_size];

            % Compute grid dimensions
            N_range = obj.bounds(2) - obj.bounds(1);
            E_range = obj.bounds(4) - obj.bounds(3);

            obj.grid_size = [ceil(E_range / cell_size), ceil(N_range / cell_size)];

            % Initialize empty cell arrays
            obj.grid_index = cell(obj.grid_size(1), obj.grid_size(2));

            % Assign each triangle to grid cells it overlaps
            n_tris = size(obj.triangles, 1);
            fprintf('Building spatial index for %d triangles...\n', n_tris);
            
            for tri = 1:n_tris
                if mod(tri, 10000) == 0
                    fprintf('Indexing triangle %d / %d (%.0f%%)\n', tri, n_tris, (tri/n_tris)*100);
                end
                
                % Get triangle bounding box
                v_idx = obj.triangles(tri, :);
                tri_verts = obj.vertices(v_idx, :);

                N_min = min(tri_verts(:, 1));
                N_max = max(tri_verts(:, 1));
                E_min = min(tri_verts(:, 2));
                E_max = max(tri_verts(:, 2));

                % Find grid cells that overlap
                j_min = max(1, floor((N_min - obj.bounds(1)) / cell_size) + 1);
                j_max = min(obj.grid_size(2), ceil((N_max - obj.bounds(1)) / cell_size));
                i_min = max(1, floor((E_min - obj.bounds(3)) / cell_size) + 1);
                i_max = min(obj.grid_size(1), ceil((E_max - obj.bounds(3)) / cell_size));

                % Add triangle to all overlapping cells
                for i = i_min:i_max
                    for j = j_min:j_max
                        obj.grid_index{i, j}(end+1) = tri;
                    end
                end
            end
            fprintf('Spatial indexing complete.\n');
        end

        function candidates = get_candidate_triangles(obj, origin, direction, max_range)
            %GET_CANDIDATE_TRIANGLES Get triangles that ray might intersect
            %
            % Uses 2D-DDA (Amanatides-Woo style) on the ray projection (N/E)
            % to visit only crossed grid cells. This avoids sampling loops,
            % unique/ismember scans, and repeated dynamic concatenation.

            candidates = zeros(1, 0);
            if max_range <= 0
                return;
            end

            origin = origin(:);
            direction = direction(:);
            dir_norm = norm(direction);
            if dir_norm < 1e-12
                return;
            end
            direction = direction / dir_norm;

            N0 = origin(1);
            E0 = origin(2);
            dN = direction(1);
            dE = direction(2);

            N_min = obj.bounds(1);
            N_max = obj.bounds(2);
            E_min = obj.bounds(3);
            E_max = obj.bounds(4);

            eps_dir = 1e-12;

            % Horizontal slab intersection against terrain XY bounds.
            if abs(dN) < eps_dir
                if N0 < N_min || N0 > N_max
                    return;
                end
                tN_min = -inf;
                tN_max = inf;
            else
                t1 = (N_min - N0) / dN;
                t2 = (N_max - N0) / dN;
                tN_min = min(t1, t2);
                tN_max = max(t1, t2);
            end

            if abs(dE) < eps_dir
                if E0 < E_min || E0 > E_max
                    return;
                end
                tE_min = -inf;
                tE_max = inf;
            else
                t1 = (E_min - E0) / dE;
                t2 = (E_max - E0) / dE;
                tE_min = min(t1, t2);
                tE_max = max(t1, t2);
            end

            t_enter = max([0, tN_min, tE_min]);
            t_exit = min([max_range, tN_max, tE_max]);
            if ~(t_enter <= t_exit)
                return;
            end

            % Initial DDA cell at horizontal entry point.
            N_entry = N0 + t_enter * dN;
            E_entry = E0 + t_enter * dE;
            cell_N = obj.cell_size(1);
            cell_E = obj.cell_size(2);

            j = floor((N_entry - N_min) / cell_N) + 1; % N-axis column
            i = floor((E_entry - E_min) / cell_E) + 1; % E-axis row

            j = min(max(j, 1), obj.grid_size(2));
            i = min(max(i, 1), obj.grid_size(1));

            if abs(dN) < eps_dir
                step_j = 0;
                tMaxN = inf;
                tDeltaN = inf;
            elseif dN > 0
                step_j = 1;
                nextN = N_min + j * cell_N;
                tMaxN = (nextN - N0) / dN;
                tDeltaN = cell_N / abs(dN);
            else
                step_j = -1;
                nextN = N_min + (j - 1) * cell_N;
                tMaxN = (nextN - N0) / dN;
                tDeltaN = cell_N / abs(dN);
            end

            if abs(dE) < eps_dir
                step_i = 0;
                tMaxE = inf;
                tDeltaE = inf;
            elseif dE > 0
                step_i = 1;
                nextE = E_min + i * cell_E;
                tMaxE = (nextE - E0) / dE;
                tDeltaE = cell_E / abs(dE);
            else
                step_i = -1;
                nextE = E_min + (i - 1) * cell_E;
                tMaxE = (nextE - E0) / dE;
                tDeltaE = cell_E / abs(dE);
            end

            % Unique candidate collection with O(1) duplicate checks.
            n_tris = size(obj.triangles, 1);
            tri_seen = false(1, n_tris);
            cap = 128;
            candidates = zeros(1, cap, 'uint32');
            n_candidates = 0;

            t_curr = t_enter;
            while i >= 1 && i <= obj.grid_size(1) && j >= 1 && j <= obj.grid_size(2) && t_curr <= t_exit
                tri_list = obj.grid_index{i, j};
                if ~isempty(tri_list)
                    tri_list = tri_list(:).';
                    add_mask = ~tri_seen(tri_list);
                    if any(add_mask)
                        new_tris = tri_list(add_mask);
                        tri_seen(new_tris) = true;

                        needed = n_candidates + numel(new_tris);
                        if needed > cap
                            while cap < needed
                                cap = cap * 2;
                            end
                            candidates(cap) = uint32(0); %#ok<AGROW>
                        end

                        candidates((n_candidates + 1):needed) = uint32(new_tris);
                        n_candidates = needed;
                    end
                end

                if tMaxN < tMaxE
                    t_curr = tMaxN;
                    tMaxN = tMaxN + tDeltaN;
                    j = j + step_j;
                elseif tMaxE < tMaxN
                    t_curr = tMaxE;
                    tMaxE = tMaxE + tDeltaE;
                    i = i + step_i;
                else
                    t_curr = tMaxN;
                    tMaxN = tMaxN + tDeltaN;
                    tMaxE = tMaxE + tDeltaE;
                    j = j + step_j;
                    i = i + step_i;
                end
            end

            candidates = double(candidates(1:n_candidates));
        end
    end
end

function cmap = terrain_colormap()
    colors = [
        0.2 0.4 0.8;
        0.1 0.6 0.2;
        0.4 0.8 0.4;
        0.8 0.7 0.4;
        0.5 0.3 0.1;
        0.5 0.5 0.5;
        0.9 0.9 0.9;
    ];
    n_colors = 256;
    x = linspace(0, 1, size(colors, 1));
    xi = linspace(0, 1, n_colors);
    cmap = interp1(x, colors, xi, 'pchip');
end
