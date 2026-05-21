%% DEMO_DEBUG
clear; clc;
clear ait_star_mex;
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
addpath(fullfile(fileparts(this_dir), 'ait_star', 'build', 'Release'));
setup_project_paths();

cfg = struct();
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 30;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 2500;
cfg.visibility_threshold = 1.1;

tm = terrain_map(cfg);

disp('tm.bounds:');
disp(tm.bounds);

% Load threat map exactly as in demo
los = load(fullfile(fileparts(this_dir), 'DEM/agri_los.mat')).los;
threat_params = struct();
terrain_max = tm.bounds(6);
terrain_min = tm.bounds(5);
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];
threat_params.num_alt_bins = 20;

threat = threat_map(tm, los, 0.1, threat_params);
radar_obj = RadarNode(1000, -500, tm.get_height(1000, -500));
radar_obj.R_max = 2000;
radar_obj.P_t = 150e3;
threat.add_radar(radar_obj);
threat.compute_map('binary', struct('use_parallel', false, 'show_progress', false));

planner_params = struct();
planner_params.risk_limit = 1.1;
planner_params.min_clearance = 50;

ait_star_mex('init_map', double(tm.Z), double(tm.resolution), double([tm.bounds(1), tm.bounds(3)]), double(planner_params.min_clearance));
ait_star_mex('init_radar3d', single(threat.risk_grid), double(threat.N_vec), double(threat.E_vec), double(threat.alt_vec), 'asl', double(3.0));

disp('Testing C++ memory:');
[t_elev, r_risk] = debug_eval_mex(double([-2000, -2000, 3880]));
fprintf('Start: (-2000, -2000) -> Elev: %.2f m, Risk: %.3f\n', t_elev, r_risk);

[t_elev, r_risk] = debug_eval_mex(double([2000, 1500, 2500]));
fprintf('Goal:  (2000,  1500)  -> Elev: %.2f m, Risk: %.3f\n', t_elev, r_risk);

fprintf('\nRANDOM SAMPLING SIMULATION:\n');
valid_count = 0;
for i = 1:100000
    px = -3000 + 6000 * rand();
    py = -3000 + 6000 * rand();
    pz = 2000 + 2800 * rand();
    
    [t, r] = debug_eval_mex(double([px, py, pz]));
    if pz >= t + 50 && r < 1.1
        valid_count = valid_count + 1;
    end
end
fprintf('Found %d valid points out of 100000!\n', valid_count);
