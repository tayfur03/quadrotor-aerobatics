%% DEMO_AIT_STAR_MEX_TEST
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
cfg.radar_alt_offset = 20;
cfg.radar_range = 2000;
cfg.rng_seed = 7;
cfg.radar_position_mode = 'manual_fraction';
cfg.num_radars = 1;
cfg.radar_manual_ne = [];
cfg.radar_manual_fraction = [0.80; 0.20];

tm = terrain_map(cfg);
los = load('DEM/agri_los.mat').los;

threat_params = struct();
terrain_max = tm.bounds(6);
terrain_min = tm.bounds(5);
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];
threat_params.num_alt_bins = 20;
threat = threat_map(tm, los, 0.1, threat_params);

if ~isempty(cfg.radar_manual_fraction)
    [rE, rN] = ndgrid(threat.E_vec, threat.N_vec);
    rE = rE(:); rN = rN(:);
    rH = tm.get_height(rN, rE);
    valid = (rH > tm.bounds(5) + 50) & (rH < tm.bounds(6) - 50);
    rN = rN(valid); rE = rE(valid); rH = rH(valid);
    for r = 1:cfg.num_radars
        f_idx = min(length(rN), max(1, round(cfg.radar_manual_fraction(r) * length(rN))));
        r_pos = [rN(f_idx); rE(f_idx); rH(f_idx) + cfg.radar_alt_offset];
        radar_obj = RadarNode(r_pos(1), r_pos(2), r_pos(3));
        radar_obj.R_max = cfg.radar_range;
        radar_obj.P_t = 150e3;
        threat.add_radar(radar_obj);
    end
end
threat.compute_map('binary', struct('use_parallel', false, 'show_progress', false));

disp('Max risk_grid:');
disp(max(threat.risk_grid(:)));

disp('Checking get_binary_horizontal_slice with 1.1');
[class_map, alt_used] = threat.get_binary_horizontal_slice(3700, 1.1);
disp('Unique class_map values (should be 1 or 2 since max risk is 1.0 and thresh is 1.1):');
disp(unique(class_map(:)));

disp('Checking get_binary_horizontal_slice with 0.5');
[class_map, alt_used] = threat.get_binary_horizontal_slice(3700, 0.5);
disp('Unique class_map values with thresh 0.5:');
disp(unique(class_map(:)));
