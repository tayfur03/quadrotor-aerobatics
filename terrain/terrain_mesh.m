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

            % Test each candidate triangle
            best_t = max_range;
            hit = false;
            tri_idx = 0;

            for i = 1:length(candidate_tris)
                tri = candidate_tris(i);
                v0 = obj.vertices(obj.triangles(tri, 1), :)';
                v1 = obj.vertices(obj.triangles(tri, 2), :)';
                v2 = obj.vertices(obj.triangles(tri, 3), :)';

                [tri_hit, t, ~, ~, ~] = moller_trumbore(origin, direction, v0, v1, v2);

                if tri_hit && t > 0 && t < best_t
                    best_t = t;
                    hit = true;
                    tri_idx = tri;
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

            for i = 1:length(candidate_tris)
                tri = candidate_tris(i);
                v0 = obj.vertices(obj.triangles(tri, 1), :)';
                v1 = obj.vertices(obj.triangles(tri, 2), :)';
                v2 = obj.vertices(obj.triangles(tri, 3), :)';

                [tri_hit, t, ~, ~, ~] = moller_trumbore(origin, direction, v0, v1, v2);

                if tri_hit && t > 0 && t < max_range
                    all_t(end+1) = t;
                    all_tri(end+1) = tri;
                end
            end

            % Sort by distance
            [distances, sort_idx] = sort(all_t);
            tri_indices = all_tri(sort_idx);

            hits = length(distances);
            points = zeros(3, hits);
            for i = 1:hits
                points(:, i) = origin + distances(i) * direction;
            end
        end

        function h = get_height(obj, N, E)
            %GET_HEIGHT Get terrain height at specified coordinates
            %   h = get_height(N, E)
            %
            %   Uses ray casting from above to find terrain height.
            %   Returns NaN for points outside terrain bounds.

            N = N(:);
            E = E(:);
            h = nan(size(N));

            for i = 1:length(N)
                if N(i) >= obj.bounds(1) && N(i) <= obj.bounds(2) && ...
                   E(i) >= obj.bounds(3) && E(i) <= obj.bounds(4)

                    % Cast ray from high altitude downward
                    origin = [N(i); E(i); obj.bounds(6) + 1000];
                    direction = [0; 0; -1];

                    [hit, point, ~, ~, ~] = obj.ray_intersect(origin, direction, 2000);

                    if hit
                        h(i) = point(3);
                    end
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
            % Uses 3D-DDA (Digital Differential Analyzer) to traverse grid cells
            % along the ray path and collect candidate triangles.

            origin = origin(:)';
            direction = direction(:)' / norm(direction(:));

            % Collect unique triangles from cells along ray path
            candidates = [];
            visited_cells = [];

            % Simple approach: sample along ray and collect cells
            step = min(obj.cell_size) / 2;
            n_steps = ceil(max_range / step);

            for s = 0:n_steps
                t = s * step;
                point = origin + direction * t;

                % Check if point is in bounds
                if point(1) < obj.bounds(1) || point(1) > obj.bounds(2) || ...
                   point(2) < obj.bounds(3) || point(2) > obj.bounds(4)
                    continue;
                end

                % Get grid cell indices
                j = floor((point(1) - obj.bounds(1)) / obj.cell_size(1)) + 1;
                i = floor((point(2) - obj.bounds(3)) / obj.cell_size(2)) + 1;

                j = max(1, min(obj.grid_size(2), j));
                i = max(1, min(obj.grid_size(1), i));

                cell_key = i * 10000 + j;
                if ~ismember(cell_key, visited_cells)
                    visited_cells(end+1) = cell_key;
                    candidates = [candidates, obj.grid_index{i, j}];
                end
            end

            % Remove duplicates
            candidates = unique(candidates);
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
