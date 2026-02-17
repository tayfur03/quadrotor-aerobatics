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
demo_params.dem_file = 'DEM/artvin.tif';  % Ignored when use_synthetic=true
demo_params.use_synthetic = true;           % Use synthetic terrain
demo_params.adaptive_res_target = 150;       % Max grid dimension for viz
demo_params.animate = true;
demo_params.use_parallel_threat = true;

% Scale-adaptive parameters
demo_params.max_threat_cells = 25000;
demo_params.max_radar_range = 10000;
demo_params.min_step_size = 30;              % Smaller steps for detailed ridge terrain
demo_params.max_step_size = 250;             % Reasonable max
demo_params.min_clearance = 30;              % 30m AGL minimum — closer to terrain = dramatic
demo_params.max_flight_alt = 400;            % Room above the ridge
demo_params.start_agl = 60;                  % Low start — close to terrain
demo_params.goal_agl = 60;                   % Low goal — close to terrain

%% 2. Environment Loading
fprintf('--- module: ENVIRONMENT ---\n');

if ~demo_params.use_synthetic && isfile(demo_params.dem_file)
    fprintf('Loading real DEM: %s\n', demo_params.dem_file);
    % Load normally
    dem_params.target_resolution = 30; 
    dem_params.fill_nodata = 'nearest';
    terrain_data = dem_loader(demo_params.dem_file, dem_params);
else
    fprintf('Using SYNTHETIC terrain: mountain\n');
    t_params = struct();
    t_params.bounds = [0, 3000, -1500, 1500];       % 3km x 3km — like demo_terrain_masking scaled up
    t_params.resolution = 12;                        % Fine mesh for smooth rendering
    t_params.type = 'mountain';                      % Central peak — same as demo_terrain_masking
    t_params.amplitude = 300;                        % Tall prominent peak (scaled from 250m)
    t_params.wavelength = 800;                       % Broad base (scaled from 300m)
    t_params.base_height = 5;                        % Low base for contrast
    t_params.seed = 42;
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

los = los_checker(tm, struct('use_mesh', true, 'mesh', mesh, 'los_eps', 1.0));
threat = threat_map(tm, los, 0.1, struct('resolution', map_res, 'alt_range', [alt_min, alt_max])); 

% --- Auto-calculate positions based on terrain bounds ---
terrain_width = tm.bounds(2) - tm.bounds(1);

% Radar on the mountain slope — same concept as demo_terrain_masking
% Mountain peak is at center of map: (1500, 0)
% Radar on the west slope, offset south — forces drone to flank around east side
peak_N = mean(tm.bounds(1:2));   % 1500
peak_E = mean(tm.bounds(3:4));   % 0

radar1_N = peak_N;               % Same N as peak
radar1_E = peak_E - 400;         % West slope (south of center in E)
radar1_h = tm.get_height(radar1_N, radar1_E);
radar1_pos = [radar1_N; radar1_E; radar1_h + 15];
radar1 = radar_site(radar1_pos, 'Slope_Radar', 'tracking');
radar1.R_max = min(2000, terrain_width * 0.3);  % Covers west approach + direct path
radar1.P_t = 250e3;  % Match demo_terrain_masking
threat.add_radar(radar1);

fprintf('Radar at [%.0f, %.0f, %.0f] (Mountain Slope), Range: %.0fm\n', radar1_pos, radar1.R_max);

% Compute Threat Map
fprintf('Computing cumulative threat map (this may take a moment)...\n');
threat.compute_map('max', struct('use_parallel', demo_params.use_parallel_threat, 'show_progress', true));

%% 4. Geometric Path Planning (RRT*)
fprintf('\n--- module: PATH PLANNING (RRT*) ---\n');

% Define Start and Goal relative to terrain bounds
% Start and Goal on opposite sides of the mountain (like demo_terrain_masking)
% Start: West of mountain (near radar side)
% Goal:  East of mountain (must flank around the peak)
margin = 0.1;  % 10% margin from edges

start_N = tm.bounds(1) + (tm.bounds(2) - tm.bounds(1)) * margin;  % Near west edge
start_E = mean(tm.bounds(3:4));  % Center E
start_h = tm.get_height(start_N, start_E);

goal_N = tm.bounds(2) - (tm.bounds(2) - tm.bounds(1)) * margin;   % Near east edge
goal_E = mean(tm.bounds(3:4));   % Center E
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
rrt_params.max_iter = 3000;       % Enough iterations for good optimization
% Adaptive step sizing: short near ridge, long in valley
rrt_params.base_step_size = min(60, max(demo_params.min_step_size, terrain_diagonal/80));
rrt_params.max_step_size = min(demo_params.max_step_size, rrt_params.base_step_size * 4);
rrt_params.min_clearance = demo_params.min_clearance;
rrt_params.goal_bias = 0.15;
rrt_params.beta = 50;              % Strong radar avoidance — makes path go around
rrt_params.shadow_bias = 0.4;      % Exploit terrain masking
rrt_params.max_flight_alt = demo_params.max_flight_alt;
rrt_params.plot_interval = 1000;
rrt_params.animate = false;

% Preferred AGL cost for smoother altitude profile
rrt_params.gamma = 0.3;            % Moderate preferred-AGL tracking
rrt_params.preferred_agl = 70;     % Target altitude above terrain [m]
rrt_params.max_climb = 60;         % Reasonable climb angle
rrt_params.max_descent = 55;       % Reasonable descent angle

fprintf('RRT* step: %.0f–%.0fm adaptive (terrain diagonal: %.0fm)\n', ...
    rrt_params.base_step_size, rrt_params.max_step_size, terrain_diagonal);

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
% MinSnap inputs: [3 x N] Waypoints
% Time allocation handled by 3-stage hybrid optimizer (Mellinger closed-form)
% optimizer inside trajectory_smoother — no manual avg_speed needed.
path_len = 0;
for i = 2:size(path_simplified, 2)
    path_len = path_len + norm(path_simplified(:,i) - path_simplified(:,i-1));
end
fprintf('Path length: %.0f m\n', path_len);

smooth_params = struct();
smooth_params.v_max = 15;  % Maximum velocity [m/s] — faster for exciting visuals
smooth_params.a_max = 6;   % Maximum acceleration [m/s^2]
smooth_params.dt = 0.05;   % Coarse dt for fast calc
smooth_params.vel_bc_mode = 'free';  % Optimizer chooses interior velocities
smooth_params.cruise_speed = 10;  % Cruise speed [m/s]
smooth_params.max_waypoints = 60;
smooth_params.max_seg_length = 150;  % Max segment length [m]
smooth_params.aggressiveness = 3.0;  % k_T=3 — aggressive but controlled

% total_time is a hint — the optimizer will adjust it
total_time = path_len / smooth_params.cruise_speed;
traj_poly = trajectory_smoother(path_simplified, total_time, smooth_params);

% Update total_time from the optimizer's output
total_time = traj_poly.t(end);

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
