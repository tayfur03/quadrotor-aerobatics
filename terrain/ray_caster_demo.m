%% Radar Ray Caster class demo
clear; clc; close all;
%terrain_data = terrain_generator(struct('type','canyon'));
% Load normally
dem_params.target_resolution = 30; 
dem_params.fill_nodata = 'nearest';
dem_params.dem_file = 'DEM/artvin.tif';
terrain_data = dem_loader(dem_params.dem_file, dem_params);
tm = terrain_map(terrain_data);

% Initialize radar ray caster with terrain map
rc = radar_ray_caster(tm); % Auto-creates mesh for fast intersection
N = 0; E = 0;
mesh = terrain_mesh(terrain_data);
agl = tm.get_height(N,E);
height = 100; % [m]
alt = agl + height;
% radar_pos = [N; E; alt];  % Radar position [N; E; altitude]
radar_pos = [0; 0; alt];  % Radar position [N; E; altitude]
bounds = 2000 * [-1 1 -1 1];

[hit, point, dist] = rc.cast_ray(radar_pos, [1; 0; -0.1], 5000);
if hit
    fprintf('Beam hits terrain at [%.1f, %.1f, %.1f]\n', point);
end
% Compute shadow map
shadow = rc.compute_shadow_map(radar_pos, bounds);
rc.plot_shadow_map(radar_pos, bounds);
% Compute visibility slice
azimuth = -90; %[deg]
max_range = 1000; %[m]
resolution = 30; %[m]
coverage = rc.compute_coverage_map(radar_pos,bounds,alt);
rc.plot_coverage_slice(radar_pos,azimuth,max_range,resolution);

%%