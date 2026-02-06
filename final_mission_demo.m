%% FINAL_MISSION_DEMO - "Final Boss" Integration Demo
%
% This script demonstrates the complete autonomous mission pipeline:
% 1. Terrain Loading (Adaptive Resolution)
% 2. Multi-Radar Threat Analysis
% 3. Geometric Path Planning (Radar-Aware RRT*)
% 4. Trajectory Smoothing (Minimum Snap Optimization)
% 5. Differential Flatness Control Verification
% 6. Full 3D Visualization & Animation
%
% Author: Quadrotor Terrain Following Project
% Date: 2026-01-31

clear; clc; close all;

%% 1. Configuration & Setup
fprintf('=== FINAL MISSION DEMO: STARTING ===\n\n');

% Add paths
addpath('terrain');
addpath('radar');
addpath('motion_planner');

% Demo Parameters
demo_params = struct();
demo_params.dem_file = 'DEM/agri.tif';  % Ignored when use_synthetic=true
demo_params.use_synthetic = false;           % Use real DEM
demo_params.adaptive_res_target = 150;       % Max grid dimension for viz
demo_params.animate = true;

% Scale-adaptive parameters for large terrain
demo_params.max_threat_cells = 25000;        % Max cells for threat map (~25k)
demo_params.max_radar_range = 10000;          % Max realistic radar range [m]
demo_params.min_step_size = 50;              % Min RRT* step size [m]
demo_params.max_step_size = 400;             % Max RRT* step size [m]
demo_params.min_clearance = 50;             % Minimum AGL clearance [m]
demo_params.max_flight_alt = 600;           % Maximum AGL altitude for planning [m]
demo_params.start_agl = 100;                 % Start altitude AGL [m]
demo_params.goal_agl = 120;                  % Goal altitude AGL [m]

%% 2. Environment Loading
fprintf('--- module: ENVIRONMENT ---\n');

if ~demo_params.use_synthetic && isfile(demo_params.dem_file)
    fprintf('Loading real DEM: %s\n', demo_params.dem_file);
    % Load normally
    dem_params.target_resolution = 30; 
    dem_params.fill_nodata = 'nearest';
    terrain_data = dem_loader(demo_params.dem_file, dem_params);
else
    fprintf('Using SYNTHETIC terrain: ridge\n');
    t_params = struct();
    t_params.bounds = [-1500, 1500, -1500, 1500];  % 3km x 3km area
    t_params.resolution = 15;
    t_params.type = 'mountain';       % North-South ridge in the middle
    t_params.amplitude = 150;      % Max height ~150m
    t_params.wavelength = 400;     % Ridge width
    t_params.base_height = 20;     % Valley floor height
    terrain_data = terrain_generator(t_params);
end

% Validate loaded terrain
if any(isnan(terrain_data.Z(:)))
    warning('DEM contains NaN values. Attempting to fill...');
    terrain_data.Z = fillmissing(terrain_data.Z, 'nearest');
end
fprintf('Terrain validated: No NaN values\n');

tm = terrain_map(terrain_data);
terrain_min = min(terrain_data.Z(:));
terrain_max = max(terrain_data.Z(:));
fprintf('Terrain stats: Bounds [%.0f, %.0f]x[%.0f, %.0f], Height: %.0f to %.0f m\n', ...
    tm.bounds(1), tm.bounds(2), tm.bounds(3), tm.bounds(4), ...
    terrain_min, terrain_max);

% Create Mesh for RRT* (Ray Casting)
fprintf('Creating terrain mesh...\n');
mesh = terrain_mesh(terrain_data);

%% 3. Threat Environment Parsing
fprintf('\n--- module: RADAR THREAT ANALYSIS ---\n');

% Initialize Threat Map with SCALE-ADAPTIVE resolution
% Instead of scaling with grid size, we limit total cells to control computation time
terrain_area = (tm.bounds(2)-tm.bounds(1)) * (tm.bounds(4)-tm.bounds(3));
target_horiz_cells = sqrt(demo_params.max_threat_cells / 10);  % ~10 altitude layers
map_res_horiz = sqrt(terrain_area) / target_horiz_cells;
map_res_horiz = max(map_res_horiz, 50);  % Minimum 50m horizontal resolution
map_res = [map_res_horiz, 50];  % [Horiz, Vert] - 50m vertical spacing

n_horiz = ceil((tm.bounds(2)-tm.bounds(1)) / map_res(1)) * ceil((tm.bounds(4)-tm.bounds(3)) / map_res(1));
n_alt = ceil(500 / map_res(2));
fprintf('Threat Map Resolution: %.0f m (Grid: ~%dx%d horiz, %d alt = %.0fk cells)\n', ...
    map_res(1), ceil(sqrt(n_horiz)), ceil(sqrt(n_horiz)), n_alt, n_horiz*n_alt/1000);

alt_min = terrain_min + demo_params.min_clearance;
alt_max = terrain_max + demo_params.max_flight_alt;
fprintf('Threat Map Altitude Range: [%.0f, %.0f] m (terrain %.0f..%.0f)\n', ...
    alt_min, alt_max, terrain_min, terrain_max);

los = los_checker(tm);
threat = threat_map(tm, los, 0.1, struct('resolution', map_res, 'alt_range', [alt_min, alt_max])); 

% --- Auto-calculate positions based on terrain bounds ---
terrain_width = tm.bounds(2) - tm.bounds(1);

% Radar 1: High point near center (DISABLED for now - uncomment to add)
radar1_N = mean(tm.bounds(1:2));
radar1_E = mean(tm.bounds(3:4));
radar1_h = tm.get_height(radar1_N, radar1_E);
radar1_pos = [radar1_N; radar1_E; radar1_h + 20];
radar1 = radar_site(radar1_pos, 'Central_Radar', 'tracking');
%radar1.R_max = min(demo_params.max_radar_range, terrain_width * 0.15);  % Capped range
radar1.R_max = 10000;
%threat.add_radar(radar1);  % DISABLED

% Radar 2: Quarter way into map
radar2_N = tm.bounds(1) + terrain_width * 0.25;
radar2_E = tm.bounds(3) + (tm.bounds(4) - tm.bounds(3)) * 0.75;
radar2_h = tm.get_height(radar2_N, radar2_E);
radar2_pos = [radar2_N; radar2_E; radar2_h + 15];
radar2 = radar_site(radar2_pos, 'Flank_Radar', 'search');
radar2.R_max = min(demo_params.max_radar_range * 0.7, terrain_width * 0.1);  % Smaller, capped
threat.add_radar(radar2);

fprintf('Radar 1 at [%.0f, %.0f, %.0f], Range: %.0fm (DISABLED)\n', radar1_pos, radar1.R_max);
fprintf('Radar 2 at [%.0f, %.0f, %.0f], Range: %.0fm\n', radar2_pos, radar2.R_max);

% Compute Threat Map
fprintf('Computing cumulative threat map (this may take a moment)...\n');
threat.compute_map();

%% 4. Geometric Path Planning (RRT*)
fprintf('\n--- module: PATH PLANNING (RRT*) ---\n');

% Define Start and Goal relative to terrain bounds
margin = 0.3;  % 10% margin from edges
% Switch the start/goal coordinates (swap North/East)
% Original positions (corner margins)
sN = tm.bounds(1) + (tm.bounds(2) - tm.bounds(1)) * margin;
sE = tm.bounds(3) + (tm.bounds(4) - tm.bounds(3)) * margin;

gN = tm.bounds(2) - (tm.bounds(2) - tm.bounds(1)) * margin;
gE = tm.bounds(4) - (tm.bounds(4) - tm.bounds(3)) * margin;

% Swap start and goal as requested
start_N = gN;
start_E = gE;
start_h = tm.get_height(start_N, start_E);

goal_N = sN;
goal_E = sE;
goal_h = tm.get_height(goal_N, goal_E);

% Validate terrain heights
if isnan(start_h)
    error('Start position [%.0f, %.0f] is outside terrain bounds!', start_N, start_E);
end
if isnan(goal_h)
    error('Goal position [%.0f, %.0f] is outside terrain bounds!', goal_N, goal_E);
end

% Use consistent AGL values from config
start_pos = [start_N; start_E; -(start_h + demo_params.start_agl)];
goal_pos = [goal_N; goal_E; -(goal_h + demo_params.goal_agl)];

fprintf('Start: [%.0f, %.0f] at %.0fm AGL (Terrain: %.0fm MSL)\n', ...
    start_N, start_E, demo_params.start_agl, start_h);
fprintf('Goal:  [%.0f, %.0f] at %.0fm AGL (Terrain: %.0fm MSL)\n', ...
    goal_N, goal_E, demo_params.goal_agl, goal_h);

% RRT* Parameters - SCALE-ADAPTIVE for large terrain
terrain_diagonal = sqrt((tm.bounds(2)-tm.bounds(1))^2 + (tm.bounds(4)-tm.bounds(3))^2);

rrt_params = struct();
rrt_params.max_iter = 6000;        % More iterations for larger terrain
rrt_params.step_size = min(demo_params.max_step_size, ...
                           max(demo_params.min_step_size, terrain_diagonal/15));
rrt_params.min_clearance = demo_params.min_clearance;  % AGL minimum (increased for mountain terrain)
rrt_params.goal_bias = 0.2;       % Slightly higher goal bias for faster convergence
rrt_params.beta = 30;              % Use fixed radar cost weight (from RRT* fixes)
rrt_params.shadow_bias = 0.4;      % Moderate shadow sampling bias
rrt_params.max_flight_alt = demo_params.max_flight_alt;   % Higher for mountain terrain
rrt_params.plot_interval = 500;    % Less frequent updates
rrt_params.animate = true;        % Disable animation for speed

% CRITICAL: Reduce altitude cost to prevent terrain-following oscillations
% High gamma causes RRT* to create many waypoints following terrain contours
rrt_params.gamma = 0.02;           % LOW: Let trajectory smoother handle altitude
rrt_params.preferred_alt = 150;    % HIGH: Prefer flying higher, fewer oscillations
rrt_params.max_climb = 60;         % Allow steeper climbs (reduces switchbacks)
rrt_params.max_descent = 55;       % Allow steeper descents

fprintf('RRT* step_size: %.0fm (terrain diagonal: %.0fm)\n', rrt_params.step_size, terrain_diagonal);

% Setup Figure for RRT - simple single view
fig_rrt = figure('Name', 'RRT* Path Planning', 'Position', [100, 100, 1000, 800]);
tm.plot(); hold on;
% Plot Radars with range circles
for i=1:length(threat.radars)
    p = threat.radars{i}.position;
    r = threat.radars{i}.R_max;
    % Radar marker
    plot3(p(1), p(2), p(3), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    % Range circle
    th = linspace(0, 2*pi, 60);
    plot3(p(1) + r*cos(th), p(2) + r*sin(th), ones(size(th))*p(3), 'r--', 'LineWidth', 1.5);
end
plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', 'MarkerSize', 14, 'MarkerFaceColor', 'g');
plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bs', 'MarkerSize', 14, 'MarkerFaceColor', 'b');
xlim([tm.bounds(1), tm.bounds(2)]); ylim([tm.bounds(3), tm.bounds(4)]);
view(35, 40);
title('RRT* Path Planning');
xlabel('North [m]'); ylabel('East [m]'); zlabel('Altitude [m]');
rrt_params.fig_handle = fig_rrt;

% Run Planner
fprintf('Running RRT*...\n');
[path_rrt, info] = rrt_star_radar(start_pos, goal_pos, tm, threat, rrt_params);

if info.success
    fprintf('RRT* Success! Cost: %.2f, Length: %.1f m\n', info.path_cost, info.path_length);
else
    error('RRT* failed to find a path. Try increasing iterations or adjusting start/goal.');
end

%% 5. Trajectory Optimization (MinSnap)
fprintf('\n--- module: TRAJECTORY OPTIMIZATION (MinSnap) ---\n');

% Simplify path with SCALE-ADAPTIVE tolerance
terrain_size = max(tm.bounds(2)-tm.bounds(1), tm.bounds(4)-tm.bounds(3));
simplify_tol = max(100, terrain_size * 0.01);  % 1.5% of terrain, min 100m
path_simplified = simplify_path(path_rrt, simplify_tol);
fprintf('Path simplified: %d -> %d waypoints (tolerance: %.0fm)\n', ...
    size(path_rrt, 2), size(path_simplified, 2), simplify_tol);

% RRT Path is [3 x N] in NED
% MinSnap inputs: [3 x N] Waypoints, Total Time
% Heuristic for time: Average speed 15 m/s
path_len = 0;
for i = 2:size(path_simplified, 2)
    path_len = path_len + norm(path_simplified(:,i) - path_simplified(:,i-1));
end
avg_speed = 5;  % Back to reasonable speed
total_time = path_len / avg_speed;
fprintf('Estimated flight time: %.1f s (Avg Speed: %.1f m/s)\n', total_time, avg_speed);

smooth_params = struct();
smooth_params.v_max = 18;  % Reasonable max velocity
smooth_params.dt = 0.05; % Coarse dt for fast calc, interpolation later if needed
smooth_params.vel_bc_mode = 'bisector';  % 'free' or 'bisector' - bisector gives smoother turns
smooth_params.cruise_speed = 12;  % Cruise speed for bisector mode [m/s]
smooth_params.max_waypoints = 20;  % Limit to prevent minsnappolytraj conditioning issues

traj_poly = trajectory_smoother(path_simplified, total_time, smooth_params);

fprintf('Trajectory generated: %d points\n', length(traj_poly.t));
fprintf('Max Vel: %.2f m/s, Max Acc: %.2f m/s^2\n', ...
    max(vecnorm(traj_poly.vel, 2, 1)), max(vecnorm(traj_poly.acc, 2, 1)));

%% 6. FINAL STATIC VISUALIZATION (No Animation)
fprintf('\n--- module: FINAL VISUALIZATION ---\n');

fig_final = figure('Name', 'Mission Summary', 'Position', [150, 100, 1100, 850]);

% Main 3D view
tm.plot();
hold on;

% Terrain bounds for axis
xlim([tm.bounds(1), tm.bounds(2)]);
ylim([tm.bounds(3), tm.bounds(4)]);

% Plot Start and Goal
plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bs', 'MarkerSize', 14, 'MarkerFaceColor', 'b', 'DisplayName', 'Goal');

% Plot Radars
for i=1:length(threat.radars)
    p = threat.radars{i}.position;
    r = threat.radars{i}.R_max;
    % Radar marker
    plot3(p(1), p(2), p(3), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    % Range circle
    th = linspace(0, 2*pi, 60);
    plot3(p(1) + r*cos(th), p(2) + r*sin(th), ones(size(th))*p(3), 'r--', 'LineWidth', 1.5);
end

% Plot Simplified RRT* Path (Green dashed)
plot3(path_simplified(1,:), path_simplified(2,:), -path_simplified(3,:), 'g--', 'LineWidth', 2, 'DisplayName', 'Simplified Path');

% Plot Simplified Waypoints (Magenta circles)
plot3(path_simplified(1,:), path_simplified(2,:), -path_simplified(3,:), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Waypoints');

% Plot velocity boundary condition arrows at waypoints
if isfield(traj_poly, 'velBC')
    velBC = traj_poly.velBC;
    wps = traj_poly.waypoints;  % Use filtered waypoints from trajectory

    % Scale factor: arrow length in meters per (m/s) of velocity
    % Adjust based on terrain size for visibility
    arrow_scale = terrain_size / 800;  % Scales with terrain

    for i = 1:size(wps, 2)
        if ~any(isnan(velBC(:, i)))  % Only plot non-NaN velocities
            wp = wps(:, i);
            v = velBC(:, i);

            % Plot arrow using quiver3: (x, y, z, u, v, w)
            % Note: D is positive down in NED, so negate for plotting
            quiver3(wp(1), wp(2), -wp(3), ...
                    v(1)*arrow_scale, v(2)*arrow_scale, -v(3)*arrow_scale, ...
                    0, 'Color', [1 0.5 0], 'LineWidth', 2.5, 'MaxHeadSize', 0.8, ...
                    'HandleVisibility', 'off');
        end
    end
    % Add single entry for legend
    quiver3(nan, nan, nan, nan, nan, nan, 0, 'Color', [1 0.5 0], 'LineWidth', 2.5, ...
            'DisplayName', 'Velocity BC');
end

% Plot MinSnap Trajectory (Blue solid)
plot3(traj_poly.pos(1,:), traj_poly.pos(2,:), -traj_poly.pos(3,:), 'b-', 'LineWidth', 2.5, 'DisplayName', 'Smooth Trajectory');

% Labels and formatting
xlabel('North [m]', 'FontSize', 11);
ylabel('East [m]', 'FontSize', 11);
zlabel('Altitude [m]', 'FontSize', 11);
legend('Location', 'northeast');
grid on;
view(40, 30);

% Title with stats
title(sprintf('Mission Path | Length: %.0f m | Time: %.1f s | Max Risk: %.2f', ...
    info.path_length, total_time, max(threat.get_risk_along_path(traj_poly.pos))), ...
    'FontSize', 12);

fprintf('\n=== MISSION COMPLETE ===\n');
fprintf('Path Length: %.1f m\n', info.path_length);
fprintf('Flight Time: %.1f s\n', total_time);
fprintf('Figure saved to current directory.\n');
