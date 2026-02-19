%% DEMO_DEM_RADAR - GeoTIFF DEM Loading and Radar Ray Casting Demonstration
%
% This demo shows how to:
%   1. Load a GeoTIFF DEM file using dem_loader
%   2. Create a terrain_map from the DEM data
%   3. Build a triangle mesh from DEM grid points
%   4. Perform radar ray casting using Möller-Trumbore algorithm
%   5. Compute and visualize radar shadow zones
%   6. Plan a path through terrain-masked areas
%
% Ray Casting Methods:
%   - MESH-BASED (default): Uses Möller-Trumbore ray-triangle intersection
%     algorithm for precise and efficient terrain intersection detection.
%     DEM grid is converted to triangle mesh (2 triangles per grid cell).
%   - GRID-BASED (fallback): Steps along ray and queries terrain height,
%     uses binary search refinement. Less accurate, may miss intersections.
%
% Reference:
%   Möller, T., & Trumbore, B. (1997). "Fast, minimum storage ray-triangle
%   intersection." Journal of Graphics Tools, 2(1), 21-28.
%
% Usage:
%   - Set the 'dem_file' variable to your GeoTIFF path, OR
%   - Use synthetic terrain for demonstration
%
% Requirements:
%   - MATLAB Mapping Toolbox (for readgeoraster with real DEMs)
%   - terrain/, radar/, motion_planner/ folders in path
%
% Author: Quadrotor Terrain Following Project

clear; clc; close all;

%% Add paths
addpath('terrain');
addpath('radar');
addpath('motion_planner');

%% Configuration
% Set this to your GeoTIFF file path, or leave empty for synthetic terrain
dem_file = 'DEM/bosphorus.tif';  % e.g., 'C:/data/elevation.tif'

% Demo mode selection
fprintf('=== GeoTIFF DEM and Radar Ray Casting Demo ===\n\n');

if ~isempty(dem_file) && isfile(dem_file)
    use_real_dem = true;
    fprintf('Using real DEM: %s\n\n', dem_file);
else
    use_real_dem = false;
    fprintf('No DEM file specified. Using synthetic terrain.\n');
    fprintf('To use a real DEM, set dem_file variable to your .tif path.\n\n');
end

%% Load or Generate Terrain
if use_real_dem
    % Load GeoTIFF DEM
    fprintf('--- Loading GeoTIFF DEM ---\n');

    dem_params = struct();
    dem_params.target_resolution = 30;  % Resample to 30m if needed
    dem_params.fill_nodata = 'nearest';

    terrain_data = dem_loader(dem_file, dem_params);
else
    % Generate synthetic terrain for demonstration
    fprintf('--- Generating Synthetic Terrain ---\n');

    terrain_params = struct();
    terrain_params.bounds = [-2000, 2000, -2000, 2000];  % 4km x 4km area
    terrain_params.resolution = 20;  % 20m resolution
    terrain_params.base_height = 100;
    terrain_params.amplitude = 300;
    terrain_params.wavelength = 800;

    terrain_data = terrain_generator('ridge', terrain_params);
end

% Create terrain_map object
tm = terrain_map(terrain_data);
fprintf('Terrain loaded: %s\n', terrain_data.type);
fprintf('Bounds: N=[%.0f, %.0f], E=[%.0f, %.0f] m\n', ...
    terrain_data.bounds(1), terrain_data.bounds(2), ...
    terrain_data.bounds(3), terrain_data.bounds(4));
fprintf('Resolution: %.1f m\n', terrain_data.resolution);
fprintf('Height range: [%.1f, %.1f] m\n\n', min(terrain_data.Z(:)), max(terrain_data.Z(:)));

%% Create Radar Ray Caster with Mesh-Based Intersection
fprintf('--- Initializing Ray Caster (Möller-Trumbore Algorithm) ---\n');

% Configure ray caster parameters
ray_caster_params = struct();
ray_caster_params.use_mesh = true;  % Use mesh-based Möller-Trumbore (default)
ray_caster_params.precision = 0.5;  % Intersection precision [m]

ray_caster = radar_ray_caster(tm, ray_caster_params);

if ray_caster.use_mesh
    fprintf('Ray caster mode: MESH-BASED (Möller-Trumbore)\n');
    fprintf('Triangle mesh: %d triangles\n', size(ray_caster.mesh.triangles, 1));
    fprintf('Spatial index cells: %d x %d\n', ...
        size(ray_caster.mesh.grid_index, 1), size(ray_caster.mesh.grid_index, 2));
else
    fprintf('Ray caster mode: GRID-BASED (stepping with binary search)\n');
    fprintf('Step size: %.1f m\n', ray_caster.step_size);
end
fprintf('Precision: %.1f m\n\n', ray_caster.precision);

%% Define Radar Position
% Place radar on high ground
[max_h, max_idx] = max(terrain_data.Z(:));
[N_idx, E_idx] = ind2sub(size(terrain_data.Z), max_idx);
radar_N = terrain_data.N_vec(N_idx);
radar_E = terrain_data.E_vec(E_idx);
radar_alt = max_h + 50;  % 50m above highest terrain

radar_pos = [radar_N; radar_E; radar_alt];
fprintf('--- Radar Configuration ---\n');
fprintf('Position: [%.0f, %.0f, %.0f] m (N, E, Alt)\n', radar_pos);
fprintf('Altitude above terrain: 50 m\n\n');

%% Demonstrate Single Ray Casting
fprintf('--- Single Ray Casting Demo ---\n');

% Cast rays at various elevations
elevations = [-20, -10, 0, 10, 20];  % degrees
azimuths = [0, 45, 90, 135, 180];    % degrees

fprintf('Casting rays from radar position:\n');
fprintf('%-10s %-10s %-10s %-15s %-10s\n', 'Azimuth', 'Elevation', 'Hits', 'Distance', 'Int. Alt');
fprintf('%s\n', repmat('-', 1, 60));

for az = azimuths
    for el = elevations
        % Convert to direction vector
        az_rad = deg2rad(az);
        el_rad = deg2rad(el);
        direction = [cos(az_rad)*cos(el_rad); sin(az_rad)*cos(el_rad); sin(el_rad)];

        [hit, point, dist] = ray_caster.cast_ray(radar_pos, direction, 5000);

        if hit
            fprintf('%-10.0f %-10.0f %-10s %-15.1f %-10.1f\n', ...
                az, el, 'YES', dist, point(3));
        else
            fprintf('%-10.0f %-10.0f %-10s %-15s %-10s\n', ...
                az, el, 'NO', '-', '-');
        end
    end
end
fprintf('\n');

%% Compute Radar Shadow Map
fprintf('--- Computing Shadow Map ---\n');

% Define analysis region (subset for speed)
analysis_bounds = [terrain_data.bounds(1) + 500, terrain_data.bounds(2) - 500, ...
                   terrain_data.bounds(3) + 500, terrain_data.bounds(4) - 500];

% ADAPTIVE RESOLUTION:
% To prevent extremely slow computation, we limit the analysis grid size.
% We target roughly a 100x100 to 150x150 grid for visualization.
max_grid_dim = 120; % Target maximum dimension (pixels)
N_len = analysis_bounds(2) - analysis_bounds(1);
E_len = analysis_bounds(4) - analysis_bounds(3);
target_res = max(N_len, E_len) / max_grid_dim;

% Use the coarser of: original resolution or the target adaptive resolution
analysis_resolution = max(terrain_data.resolution, target_res);

fprintf('Adaptive analysis resolution: %.2f m (Grid approx %dx%d)\n', ...
    analysis_resolution, round(N_len/analysis_resolution), round(E_len/analysis_resolution));

[shadow_map, N_grid, E_grid] = ray_caster.compute_shadow_map(...
    radar_pos, analysis_bounds, analysis_resolution);

shadow_percent = 100 * sum(shadow_map(:)) / numel(shadow_map);
fprintf('Shadow coverage: %.1f%% of analyzed area\n\n', shadow_percent);

%% Compute Coverage at Flight Altitude
fprintf('--- Computing Coverage at Flight Altitude ---\n');

flight_alt = 150;  % 150m flight altitude
[coverage_map, ~, ~] = ray_caster.compute_coverage_map(...
    radar_pos, analysis_bounds, flight_alt, analysis_resolution);

coverage_percent = 100 * sum(coverage_map(:)) / numel(coverage_map);
fprintf('Coverage at %.0fm altitude: %.1f%% visible to radar\n\n', ...
    flight_alt, coverage_percent);

%% Visualization
fprintf('--- Generating Visualizations ---\n');

% Figure 1: Terrain Overview
figure('Name', 'Terrain and Radar Position', 'Position', [50, 50, 1000, 800]);

subplot(2, 2, [1, 3]);
tm.plot();
hold on;
plot3(radar_pos(1), radar_pos(2), radar_pos(3), 'rp', ...
    'MarkerSize', 20, 'MarkerFaceColor', 'r', 'DisplayName', 'Radar');
title('Terrain with Radar Position');
legend('Location', 'best');

subplot(2, 2, 2);
% Shadow map
% Note: imagesc(x,y,C) -> x=East, y=North
if isvector(E_grid) && isvector(N_grid)
    imagesc(E_grid, N_grid, double(~shadow_map));
else
    % If grids are matrices (meshgrid), taking valid vectors
    imagesc(E_grid(1,:), N_grid(:,1), double(~shadow_map)); 
end
hold on;
plot(radar_pos(2), radar_pos(1), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
set(gca, 'YDir', 'normal');
colormap(gca, [0.5 0.5 0.5; 0.2 0.8 0.2]);
xlabel('East [m]');
ylabel('North [m]');
title('Terrain Shadow Map (Green = Visible)');
axis equal tight;
colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Shadow', 'Visible'});

subplot(2, 2, 4);
% Coverage at flight altitude
if isvector(E_grid) && isvector(N_grid)
    imagesc(E_grid, N_grid, double(coverage_map));
else
    imagesc(E_grid(1,:), N_grid(:,1), double(coverage_map));
end
hold on;
plot(radar_pos(2), radar_pos(1), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
set(gca, 'YDir', 'normal');
colormap(gca, [0.8 0.2 0.2; 0.2 0.8 0.2]);
xlabel('East [m]');
ylabel('North [m]');
title(sprintf('Radar Coverage at %.0fm Alt (Green = Visible)', flight_alt));
axis equal tight;
colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Hidden', 'Visible'});

% Figure 2: Ray Casting Visualization
figure('Name', 'Radar Ray Casting', 'Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
% Plot terrain surface
% Match surf dimensions to Z (Rows=N, Cols=E)
[E_mesh, N_mesh] = meshgrid(terrain_data.E_vec, terrain_data.N_vec);
surf(E_mesh, N_mesh, terrain_data.Z, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
colormap(gca, terrain_colormap());

% Plot radar (X=E, Y=N)
plot3(radar_pos(2), radar_pos(1), radar_pos(3), 'rp', ...
    'MarkerSize', 15, 'MarkerFaceColor', 'r');

% Cast and plot rays in a fan pattern
fan_azimuths = 0:15:345;
fan_elevation = -5;  % 5 degrees below horizontal
ray_colors = jet(length(fan_azimuths));

for i = 1:length(fan_azimuths)
    az_rad = deg2rad(fan_azimuths(i));
    el_rad = deg2rad(fan_elevation);
    direction = [cos(az_rad)*cos(el_rad); sin(az_rad)*cos(el_rad); sin(el_rad)];

    [hit, point, ~] = ray_caster.cast_ray(radar_pos, direction, 3000);

    % Plot ray (X=E, Y=N)
    plot3([radar_pos(2), point(2)], [radar_pos(1), point(1)], ...
          [radar_pos(3), point(3)], '-', 'Color', ray_colors(i,:), 'LineWidth', 1);

    if hit
        plot3(point(2), point(1), point(3), 'o', ...
            'Color', ray_colors(i,:), 'MarkerSize', 5, 'MarkerFaceColor', ray_colors(i,:));
    end
end

xlabel('East [m]');
ylabel('North [m]');
zlabel('Altitude [m]');
title('Radar Ray Fan (Elevation = -5 deg)');
view(30, 45);
axis equal;
grid on;

subplot(1, 2, 2);
% Vertical slice along North axis
ray_caster.plot_coverage_slice(radar_pos, 0, 4000, 50);

% Figure 3: Horizon Profile
figure('Name', 'Radar Horizon Profile', 'Position', [150, 150, 800, 600]);

[horizon_angles, azimuths_out] = ray_caster.get_horizon_profile(radar_pos, 72, 4000);

% Polar plot of horizon
polarplot(deg2rad(azimuths_out), horizon_angles + 90, 'b-', 'LineWidth', 2);
hold on;
polarplot(deg2rad(azimuths_out), 90 * ones(size(azimuths_out)), 'k--', 'LineWidth', 1);
title('Terrain Horizon Profile (0 = horizontal)');
rlim([0 180]);
rticks([0 30 60 90 120 150 180]);
rticklabels({'-90', '-60', '-30', '0', '30', '60', '90'});

%% Path Planning Through Shadow Zones
fprintf('--- Path Planning Through Shadow Zones ---\n');

% Define start and goal in shadow zones
% Find points in shadow for start/goal
shadow_indices = find(shadow_map);
if length(shadow_indices) >= 2
    % Start from one shadow area
    start_idx = shadow_indices(1);
    [start_N_idx, start_E_idx] = ind2sub(size(shadow_map), start_idx);

    % Goal in another shadow area
    goal_idx = shadow_indices(end);
    [goal_N_idx, goal_E_idx] = ind2sub(size(shadow_map), goal_idx);

    N_steps = size(shadow_map, 1);
    E_steps = size(shadow_map, 2);
    N_vec_analysis = linspace(analysis_bounds(1), analysis_bounds(2), N_steps);
    E_vec_analysis = linspace(analysis_bounds(3), analysis_bounds(4), E_steps);

    start_N = N_vec_analysis(start_N_idx);
    start_E = E_vec_analysis(start_E_idx);
    start_h = tm.get_height(start_N, start_E);
    start_pos = [start_N; start_E; -(start_h + 50)];  % 50m AGL, NED D

    goal_N = N_vec_analysis(goal_N_idx);
    goal_E = E_vec_analysis(goal_E_idx);
    goal_h = tm.get_height(goal_N, goal_E);
    goal_pos = [goal_N; goal_E; -(goal_h + 50)];

    fprintf('Start: [%.0f, %.0f, %.0f] m (NED)\n', start_pos);
    fprintf('Goal:  [%.0f, %.0f, %.0f] m (NED)\n', goal_pos);

    % Simple straight-line path for visualization
    n_points = 50;
    path = [linspace(start_pos(1), goal_pos(1), n_points);
            linspace(start_pos(2), goal_pos(2), n_points);
            linspace(start_pos(3), goal_pos(3), n_points)];

    % Check visibility along path
    visibility_along_path = zeros(1, n_points);
    for i = 1:n_points
        path_pos_alt = [path(1,i); path(2,i); -path(3,i)];  % Convert to altitude
        [hit, ~, ~] = ray_caster.cast_ray_to_target(radar_pos, path_pos_alt);
        visibility_along_path(i) = ~hit;  % Visible if no terrain hit
    end

    visible_percent = 100 * sum(visibility_along_path) / n_points;
    fprintf('Path visibility to radar: %.1f%% of path is visible\n', visible_percent);

    % Plot path on terrain
    figure('Name', 'Path Through Shadow Zones', 'Position', [200, 200, 1000, 800]);

    tm.plot();
    hold on;

    % Plot radar
    plot3(radar_pos(1), radar_pos(2), radar_pos(3), 'rp', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'r', 'DisplayName', 'Radar');

    % Plot path colored by visibility
    % Dummy plots for Legend
    h_vis = plot3(nan, nan, nan, '-', 'Color', [1 0 0], 'LineWidth', 3, 'DisplayName', 'Visible Path');
    h_hid = plot3(nan, nan, nan, '-', 'Color', [0 0.8 0], 'LineWidth', 3, 'DisplayName', 'Hidden Path');

    for i = 1:n_points-1
        if visibility_along_path(i)
            color = [1 0 0];  % Red = visible
        else
            color = [0 0.8 0];  % Green = hidden
        end
        plot3([path(1,i), path(1,i+1)], [path(2,i), path(2,i+1)], ...
              [-path(3,i), -path(3,i+1)], '-', 'Color', color, 'LineWidth', 3, 'HandleVisibility', 'off');
    end

    % Start/goal markers
    plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', ...
        'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bs', ...
        'MarkerSize', 15, 'MarkerFaceColor', 'b', 'DisplayName', 'Goal');

    title('Flight Path (Green = Hidden, Red = Visible to Radar)');
    legend([h_vis, h_hid], 'Location', 'best');
else
    fprintf('Not enough shadow area for path planning demo.\n');
end

%% Summary
fprintf('\n=== Demo Complete ===\n');
fprintf('Created visualizations:\n');
fprintf('  - Terrain overview with radar position\n');
fprintf('  - Shadow map and coverage map\n');
fprintf('  - Ray casting fan visualization\n');
fprintf('  - Horizon profile polar plot\n');
fprintf('  - Path planning through shadow zones\n');

%% Helper function for terrain colormap
function cmap = terrain_colormap()
    colors = [
        0.2 0.4 0.8;   % Deep Blue
        0.1 0.6 0.2;   % Dark Green
        0.4 0.8 0.4;   % Light Green
        0.8 0.7 0.4;   % Sand
        0.5 0.3 0.1;   % Brown
        0.5 0.5 0.5;   % Gray
        0.9 0.9 0.9;   % White
    ];
    n_colors = 256;
    x = linspace(0, 1, size(colors, 1));
    xi = linspace(0, 1, n_colors);
    cmap = interp1(x, colors, xi, 'pchip');
end
