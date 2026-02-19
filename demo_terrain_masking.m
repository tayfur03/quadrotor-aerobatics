%% DEMO_TERRAIN_MASKING - Terrain Masking Visualization
%
% This demo specifically highlights the "Terrain Masking" capability:
% The UAV exploits terrain features (like valleys and ridges) to block
% the line-of-sight from a hostile radar, flying in its "shadow".
%
% Visualizations:
%   1. Real-time RRT* optimization (Tree growth)
%   2. 3D Flight Animation (Drone perspective)
%   3. Detection Probability & Terrain Clearance Analysis
%
% Author: Quadrotor Terrain Following Project

clear; clc; close all;

%% Add paths
addpath('terrain');
addpath('radar');
addpath('safety');
addpath('motion_planner');

%% ==================== SCENARIO SETUP ====================
fprintf('=== Terrain Masking Demonstration ===\n');

% 1. Generate "Single Mountain" Terrain
% A prominent central peak to demonstrate masking by a solid obstacle.
terrain_params.bounds = [0, 800, -400, 400];
terrain_params.resolution = 10;
terrain_params.amplitude = 250;     % Very tall mountain
terrain_params.wavelength = 300;    % Broad base
terrain_params.type = 'mountain';

terrain_data = terrain_generator(terrain_params);
tm = terrain_map(terrain_data);

% terrain_data = dem_loader('DEM/artvin.tif');
% tm = terrain_map(terrain_data);
% Bunun yerine radar_ray_caster gelecek. meshli falan
los = los_checker(tm);

% 2. Setup Hostile Radar
% Radar on one slope
% We place the radar on the WEST slope.
% The central peak is at (400, 0).
% We place radar at (400, -200).
radar_x = 400;
radar_y = -200;
radar_z = tm.get_height(radar_x, radar_y) + 10;

% This radar covers the West side and the direct path (center).
% The East side (E > 0) should be masked by the peak.
rcs = 0.1;
threat = threat_map(tm, los, rcs);
threat_params.alt_range = [0, 400];

radar1 = radar_site([radar_x; radar_y; radar_z], 'Slope-Radar', 'tracking');
radar1.R_max = 700;
radar1.P_t = 250e3;
threat.add_radar(radar1);

fprintf('Computing Threat Map (Calculating Shadows behind the Mountain)...\n');
threat.compute_map('max');

%% ==================== ANIMATION 1: RRT* OPTIMIZATION ====================
fprintf('\nStep 1: Planning Path with Visualization...\n');

% Start and Goal
% Start (South) -> Goal (North)
% Direct path passes (400, 0) which is exposed or blocked.
% The West side is covered by Radar.
% Ideally, the drone should flank around the EAST side (Back of the mountain relative to radar).
start_pos = [50; 0; -(tm.get_height(50, 0) + 20)];
goal_pos = [750; 0; -(tm.get_height(750, 0) + 20)];

% Setup Figure for RRT Animation
h_rrt = figure('Name', 'RRT* Optimization: Flanking Maneuver', 'Position', [50, 50, 1000, 700]);
tm.plot(h_rrt);
hold on;
% Plot Radar
plot3(radar1.position(1), radar1.position(2), radar1.position(3), ...
    'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(radar1.position(1), radar1.position(2), radar1.position(3)+50, 'RADAR', 'Color', 'r', 'FontSize', 12);

title('RRT* Task: Flank the Mountain to Stay Hidden');
view(0, 90); % Top view clarifies the "Around" path
xlabel('North'); ylabel('East');

% RRT Parameters
planner_params.max_iter = 2500; % Need enough iterations to find the long way around
planner_params.step_size = 40;
planner_params.goal_bias = 0.2;
planner_params.alpha = 2.0;
planner_params.beta = 50;
planner_params.gamma = 0.1;         % Altitude cost
planner_params.min_clearance = 15;
planner_params.preferred_alt = 25;

% Visualization params
planner_params.animate = true;
planner_params.plot_interval = 25;
planner_params.fig_handle = h_rrt;

% Run Planner
[path_rrt, plan_info] = rrt_star_radar(start_pos, goal_pos, tm, threat, planner_params);

% Finalize RRT Plot (keep it open for reference)
figure(h_rrt);
plot3(path_rrt(1,:), path_rrt(2,:), -path_rrt(3,:), 'g-', 'LineWidth', 3);
plot3(start_pos(1), start_pos(2), -start_pos(3), 'yo', 'MarkerSize', 10, 'LineWidth', 2);
plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'mx', 'MarkerSize', 10, 'LineWidth', 2);
title('RRT* Optimization Complete');
drawnow;

%% ==================== TRAJECTORY GENERATION ====================
fprintf('\nStep 2: Smoothing Trajectory...\n');
waypoints_smooth = smooth_path_geometric(path_rrt, tm, 3);

% Filter out duplicate waypoints
wps_unique = waypoints_smooth(:, 1);
for i = 2:size(waypoints_smooth, 2)
    if norm(waypoints_smooth(:, i) - wps_unique(:, end)) > 0.1
        wps_unique = [wps_unique, waypoints_smooth(:, i)];
    end
end
waypoints_smooth = wps_unique;

% Simple time allocation
dists = sqrt(sum(diff(waypoints_smooth, 1, 2).^2, 1));
speed = 12; % m/s
dt_segments = max(dists / speed, 0.1);
times = [0, cumsum(dt_segments)];
total_duration = times(end);

% Min snap
[~, ~, ~, ~, ~, pp] = minsnappolytraj(waypoints_smooth, times, size(waypoints_smooth,2)*4, ...
    'MinSegmentTime', 0.1);

% Resample for sim
dt = 0.05;
t_sim = 0:dt:total_duration;
xref = ppval(pp, t_sim);
x_path = xref(1,:);
y_path = xref(2,:);
z_path = -xref(3,:);

%% ==================== ANALYSIS ====================

% 1. Detection Probability (Risk) & Altitude - Combined Plot
figure('Name', 'Terrain Masking Analysis', 'Position', [100, 100, 1000, 600]);

subplot(2,1,1);
[risks, ~] = threat.get_risk_along_path(xref);
area(t_sim, risks, 'FaceColor', 'r', 'FaceAlpha', 0.3);
hold on;
yline(0.1, 'k--', 'Stealth Threshold');
grid on;
ylabel('Detection Prob. P_{det}');
title('Risk Profile: Red Areas indicate Detection Exposure');
ylim([0 1]);

subplot(2,1,2);
path_dist = [0, cumsum(vecnorm(diff(xref, 1, 2), 2, 1))];
terrain_h_path = arrayfun(@(i) tm.get_height(xref(1,i), xref(2,i)), 1:length(t_sim));
path_alt = -xref(3,:);

plot(t_sim, path_alt, 'b-', 'LineWidth', 2); hold on;
plot(t_sim, terrain_h_path, 'k-', 'LineWidth', 1.5);
fill([t_sim, fliplr(t_sim)], [terrain_h_path, zeros(size(terrain_h_path))], ...
    [0.6 0.5 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xlabel('Time [s]');
ylabel('Altitude [m]');
legend('UAV Altitude', 'Terrain Profile');
title('Terrain Following Performance');
grid on;

%% ==================== VIDEO SETUP ====================
% video_filename = 'terrain_masking_demo';
% v = VideoWriter(video_filename, 'MPEG-4');
% v.FrameRate = 30;
% open(v);
% fprintf('\nRecording video to %s.mp4...\n', video_filename);

%% ==================== ANIMATION 2: FLIGHT VISUALIZATION ====================
fprintf('\nStep 3: 3D Flight Animation...\n');
h_anim = figure('Name', '3D Flight Animation', 'Position', [100, 100, 1200, 800]);

% Static Setup
tm.plot(h_anim); hold on;
plot3(radar1.position(1), radar1.position(2), radar1.position(3), ...
    'r^', 'MarkerSize', 20, 'MarkerFaceColor', 'r');

% Draw Radar Range Sphere (Wireframe, transparent)
[sx, sy, sz] = sphere(30);
radar_r = radar1.R_max;
surf(sx*radar_r + radar1.position(1), ...
    sy*radar_r + radar1.position(2), ...
    sz*radar_r + radar1.position(3), ...
    'FaceColor', 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'r', 'EdgeAlpha', 0.1);

% Draw Line of Sight rays from Radar to path (occasionally) to show masking?
% Maybe too cluttered. Let's just animate the drone.

view(45, 30);
grid on;
axis equal;
title('Simulation: Terrain Masking Flight');

% Drone Marker
h_drone = plot3(x_path(1), y_path(1), z_path(1), 'b.', 'MarkerSize', 25);
h_trail = plot3(x_path(1), y_path(1), z_path(1), 'b-', 'LineWidth', 1.5);
h_shadow_line = plot3([radar1.position(1), x_path(1)], ...
    [radar1.position(2), y_path(1)], ...
    [radar1.position(3), z_path(1)], 'k--', 'LineWidth', 0.5);

% status text
h_text = text(50, 50, 500, 'Initializing...', 'FontSize', 12, 'BackgroundColor', 'w');

fprintf('Starting animation...\n');

for k = 1:length(t_sim)
    % Update Drone
    set(h_drone, 'XData', x_path(k), 'YData', y_path(k), 'ZData', z_path(k));
    set(h_trail, 'XData', x_path(1:k), 'YData', y_path(1:k), 'ZData', z_path(1:k));

    % Check LOS status
    los_status = los.has_los(radar1.position, [x_path(k); y_path(k); z_path(k)]);
    p_det = risks(k);

    if los_status
        % Visible!
        color = 'r';
        status_str = sprintf('DETECTED! P_{det}: %.2f', p_det);
        set(h_shadow_line, 'Color', 'r', 'LineStyle', '-');
        set(h_drone, 'Color', 'r');
    else
        % Masked
        color = 'g';
        status_str = 'MASKED (Stealth)';
        set(h_shadow_line, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        set(h_drone, 'Color', 'b');
    end

    % Update LOS line
    set(h_shadow_line, 'XData', [radar1.position(1), x_path(k)], ...
        'YData', [radar1.position(2), y_path(k)], ...
        'ZData', [radar1.position(3), z_path(k)]);

    set(h_text, 'String', sprintf('Time: %.1fs\nStatus: %s\nAlt: %.1fm', t_sim(k), status_str, z_path(k)), ...
        'EdgeColor', color, 'LineWidth', 2);

    % Camera follow (optional)
    % xlim([x_path(k)-200, x_path(k)+200]);
    % ylim([y_path(k)-200, y_path(k)+200]);

    drawnow limitrate;

    % Record Frame
    frame = getframe(gcf);
    %writeVideo(v, frame);
end

close(v);
fprintf('Demo Complete. Video saved.\n');

%% ==================== REPORT FIGURE ====================
% Create a static, high-quality visualization for the report
fprintf('\nGenerating static report figure...\n');
fig_report = figure('Name', 'Mission Scenario Report', 'Position', [100, 100, 1000, 600], 'Color', 'w');

% 1. Plot Terrain with TrueColor (custom RGB) to avoid colormap conflicts
[N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);

% Define terrain colors (High contrast palette)
t_colors = [
    0.2 0.4 0.8;   % Deep Blue
    0.1 0.6 0.2;   % Dark Green
    0.4 0.8 0.4;   % Light Green
    0.8 0.7 0.4;   % Sand/Dirt
    0.5 0.3 0.1;   % Dark Brown
    0.5 0.5 0.5;   % Gray
    0.9 0.9 0.9;   % White
    ];
t_vals = linspace(0, 1, size(t_colors, 1));
t_cmap = interp1(t_vals, t_colors, linspace(0,1,256), 'pchip');

% Normalize Z to map to colors
z_min = min(tm.Z(:));
z_max = max(tm.Z(:));
Z_norm = (tm.Z - z_min) / (z_max - z_min);
% Clamp just in case
Z_norm = max(0, min(1, Z_norm));
Z_ind = round(Z_norm * 255) + 1;
% Create TrueColor CData
C_true = reshape(t_cmap(Z_ind, :), [size(tm.Z,1), size(tm.Z,2), 3]);

% Plot surface with explicit CData (ignores figure colormap)
surf(N_grid, E_grid, tm.Z, C_true, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
hold on;
contour3(N_grid, E_grid, tm.Z, 15, 'k', 'LineWidth', 0.5);

% Lighting for the terrain
light('Position', [-1 -1 1], 'Style', 'infinite');
lighting gouraud;
material dull;

% 2. Plot Radar and Range
plot3(radar1.position(1), radar1.position(2), radar1.position(3), ...
    'r^', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(radar1.position(1), radar1.position(2), radar1.position(3)+60, 'RADAR', ...
    'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Transparent Sphere
[sx, sy, sz] = sphere(50);
radar_r = radar1.R_max;
surf(sx*radar_r + radar1.position(1), ...
    sy*radar_r + radar1.position(2), ...
    sz*radar_r + radar1.position(3), ...
    'FaceColor', 'r', 'FaceAlpha', 0.03, 'EdgeColor', 'none');

% 3. Plot Path (Color-coded by Risk using Figure Colormap)
% Z-values need to be offset slightly to be visible above terrain
z_offset_vis = z_path + 5;
p_line = patch([x_path nan], [y_path nan], [z_offset_vis nan], [risks nan], ...
    'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 4);

% Set Figure Colormap for the PATH (Green -> Yellow -> Red)
% Since Terrain uses TrueColor CData, this colormap only affects the path
c_risk = [0 1 0; 1 1 0; 1 0 0];
cmap_risk = interp1([0 0.5 1], c_risk, linspace(0,1,256));
colormap(fig_report, cmap_risk);
clim([0 1]);
cbar = colorbar;
cbar.Label.String = 'Detection Probability (P_{det})';

% 4. Start and Goal
plot3(start_pos(1), start_pos(2), -start_pos(3), 'go', 'MarkerSize', 10, 'LineWidth', 3);
text(start_pos(1), start_pos(2), -start_pos(3)+40, 'START', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');

plot3(goal_pos(1), goal_pos(2), -goal_pos(3), 'bo', 'MarkerSize', 10, 'LineWidth', 3);
text(goal_pos(1), goal_pos(2), -goal_pos(3)+40, 'GOAL', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');

% 5. Aesthetics
view(-45, 45); % Isometric-like view
grid on;
axis equal;
xlabel('North [m]'); ylabel('East [m]'); zlabel('Altitude [m]');
title('Terrain Masking Mission Analysis');

% Add a text box with path metrics
dim = [0.15 0.75 0.2 0.1];
metrics_str = sprintf('Path Length: %.0f m\nCheckpoints: %d\nAvg P_{det}: %.2f', ...
    path_dist(end), length(waypoints_smooth), mean(risks));
annotation('textbox', dim, 'String', metrics_str, 'FitBoxToText','on', 'BackgroundColor', 'w');

fprintf('Report figure created. You can save this using File > Save As...\n');
