%% DEMO_FMM_BIT_STAR
% FMM-guided BIT* demo on terrain-aware binary radar masking.

clear; clc; close all;
clear fmm_bit_star_planner compute_stealth_corridor compute_fmm_heuristic;
rehash;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
setup_project_paths();

fprintf('=== FMM-BIT* Hybrid Demo ===\n');
fprintf('Using fmm_bit_star_planner: %s\n', which('fmm_bit_star_planner'));
fprintf('Using compute_fmm_heuristic: %s\n', which('compute_fmm_heuristic'));

cfg = struct();
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 30;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 2500;
cfg.visibility_threshold = 0.5;
cfg.start_goal_agl = 80;
cfg.corridor_alpha = 2.3;
cfg.corridor_lambda = 1.8;
cfg.corridor_r_min = 120;
cfg.corridor_r_max = 450;
cfg.corridor_max_sweeps = 60;
cfg.threat_vert_res = 40;
cfg.use_parallel_threat = false;
cfg.radar_alt_offset = 20;
cfg.radar_range = 3000;

%% 1) Terrain
terrain_data = [];
try
    if exist(cfg.dem_file, 'file')
        dem_params = struct();
        dem_params.target_resolution = cfg.dem_target_resolution;
        dem_params.fill_nodata = cfg.dem_fill_nodata;
        if ~isempty(cfg.dem_crop_half_size) && cfg.dem_crop_half_size > 0
            h = cfg.dem_crop_half_size;
            dem_params.crop_bounds = [-h, h, -h, h];
        end
        terrain_data = dem_loader(cfg.dem_file, dem_params);
    end
catch ME
    warning('DEM load failed (%s). Using synthetic terrain fallback.', ME.message);
    terrain_data = [];
end

if isempty(terrain_data)
    [X, Y] = meshgrid(linspace(0, 5000, 100), linspace(0, 5000, 100));
    Z = 200 + 150 * sin(X / 800) .* cos(Y / 600) + 80 * randn(size(X)) * 0.05;
    terrain_data = struct();
    terrain_data.N_vec = X(1, :);
    terrain_data.E_vec = Y(:, 1)';
    terrain_data.Z = Z;
    terrain_data.bounds = [min(terrain_data.N_vec), max(terrain_data.N_vec), ...
                           min(terrain_data.E_vec), max(terrain_data.E_vec)];
    terrain_data.resolution = mean(diff(terrain_data.N_vec));
    terrain_data.type = 'synthetic_fmm_bit';
end

tm = terrain_map(terrain_data);
los = los_checker(tm);

%% 2) Binary radar LOS map
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));

threat_params = struct();
threat_params.resolution = [max(40, diag_len / 120), cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];

threat = threat_map(tm, los, 0.1, threat_params);
radar_ne = [tm.bounds(1) + 0.22 * N_span, tm.bounds(1) + 0.78 * N_span; ...
            tm.bounds(3) + 0.25 * E_span, tm.bounds(3) + 0.72 * E_span];

for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + cfg.radar_alt_offset;
    radar_obj = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar_obj.R_max = cfg.radar_range;
    radar_obj.P_t = 150e3;
    threat.add_radar(radar_obj);
end

fprintf('Computing binary threat map...\n');
t_threat = tic;
threat.compute_map('binary', struct('use_parallel', cfg.use_parallel_threat, 'show_progress', false));
fprintf('Threat map compute time: %.2f s\n', toc(t_threat));

radar_los_map = threat.risk_grid >= cfg.visibility_threshold;

%% 3) Start / goal
[start_NE, start_h, start_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(1) + 0.08 * N_span; tm.bounds(3) + 0.10 * E_span], cfg.start_goal_agl, cfg.visibility_threshold);
[goal_NE, goal_h, goal_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(2) - 0.08 * N_span; tm.bounds(4) - 0.10 * E_span], cfg.start_goal_agl, cfg.visibility_threshold);

start_world = [start_NE; start_h + cfg.start_goal_agl];
goal_world = [goal_NE; goal_h + cfg.start_goal_agl];

fprintf('Start risk: %.2f | Goal risk: %.2f\n', start_risk, goal_risk);

%% 4) Strategic corridor
corr_params = struct();
corr_params.N_vec = threat.N_vec;
corr_params.E_vec = threat.E_vec;
corr_params.alt_vec = threat.alt_vec;
corr_params.max_sweeps = cfg.corridor_max_sweeps;
corr_params.lambda = cfg.corridor_lambda;
corr_params.r_min = cfg.corridor_r_min;
corr_params.r_max = cfg.corridor_r_max;
corr_params.max_backtrack_points = 5000;

fprintf('Computing stealth corridor...\n');
t_corr = tic;
corridor = compute_stealth_corridor( ...
    threat.risk_grid, start_world, goal_world, cfg.visibility_threshold, cfg.corridor_alpha, corr_params);
fprintf('Corridor build time: %.2f s | skeleton points: %d\n', toc(t_corr), size(corridor.skeleton, 2));

goal_idx = round(corridor.grid.goal_idx);
if ~isfinite(corridor.T(goal_idx(1), goal_idx(2), goal_idx(3)))
    relaxed_threshold = min(0.7, cfg.visibility_threshold + 0.15);
    warning('Strict corridor disconnected. Retrying at threshold %.2f.', relaxed_threshold);
    corridor = compute_stealth_corridor( ...
        threat.risk_grid, start_world, goal_world, relaxed_threshold, cfg.corridor_alpha, corr_params);
end

%% 5) FMM-BIT*
planner_params = struct();
planner_params.max_batches = 120;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.08;
planner_params.stealth_weight = 3.0;
planner_params.strict_zero_risk_tol = 1e-6;
planner_params.edge_check_samples = 9;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.eta_radius = 1.8;
planner_params.max_time = 45;
planner_params.min_clearance = 10;
planner_params.w_mask = 1.0;
planner_params.w_alt = 0.08;
planner_params.h_opt = cfg.start_goal_agl;
planner_params.h_scale = 200;
planner_params.dz = mean(diff(threat.alt_vec));
planner_params.mem_limit_gb = 1.5;
planner_params.fmm_point_clearance = 10;
planner_params.fmm_bounds_pad_xy = 180;
planner_params.fmm_bounds_pad_z = 80;
[Nq_fmm, Eq_fmm] = meshgrid(threat.N_vec, threat.E_vec);
Z_fmm = tm.get_height(Nq_fmm, Eq_fmm);
planner_params.fmm_terrain_map = struct( ...
    'Z', Z_fmm, ...
    'N_vec', threat.N_vec, ...
    'E_vec', threat.E_vec, ...
    'x_vec', threat.N_vec, ...
    'y_vec', threat.E_vec, ...
    'dx', mean(diff(threat.N_vec)), ...
    'dy', mean(diff(threat.E_vec)), ...
    'alt_vec', threat.alt_vec);
planner_params.fmm_radar_los_map = radar_los_map;

[path_world, final_cost, info] = fmm_bit_star_planner(start_world, goal_world, corridor, planner_params);

%% 6) Figures
figure('Name', 'FMM-BIT* Path on Terrain', 'Position', [100, 80, 1180, 820], 'Color', 'w');
tm.plot();
hold on;
for r = 1:length(threat.radars)
    rr = threat.radars{r};
    plot3(rr.position(1), rr.position(2), rr.position(3), 'kd', ...
        'MarkerFaceColor', [1.0, 0.55, 0.1], 'MarkerSize', 9);
end
if ~isempty(path_world)
    plot3(path_world(1, :), path_world(2, :), path_world(3, :), ...
        '-', 'Color', [0.92, 0.15, 0.18], 'LineWidth', 2.6);
end
plot3(start_world(1), start_world(2), start_world(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot3(goal_world(1), goal_world(2), goal_world(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
title('Terrain Mesh With FMM-BIT* Path');
view(38, 30);
axis tight;
grid on;

figure('Name', 'FMM Heuristic Diagnostics', 'Position', [130, 100, 1240, 500], 'Color', 'w');
subplot(1, 2, 1);
if info.fmm_used && ~isempty(info.fmm_T_bwd)
    if ~isempty(path_world)
        slice_alt = median(path_world(3, :));
    else
        slice_alt = median([start_world(3), goal_world(3)]);
    end
    [~, k_slice] = min(abs(info.fmm_meta.grid_z - slice_alt));
    T_slice = squeeze(info.fmm_T_bwd(:, :, k_slice));
    imagesc(info.fmm_meta.grid_x, info.fmm_meta.grid_y, T_slice');
    set(gca, 'YDir', 'normal');
    axis equal tight;
    hold on;
    if ~isempty(path_world)
        plot(path_world(1, :), path_world(2, :), 'w-', 'LineWidth', 2.0);
    end
    plot(start_world(1), start_world(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(goal_world(1), goal_world(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    title(sprintf('T_{bwd} Slice at %.1f m', info.fmm_meta.grid_z(k_slice)));
    xlabel('North [m]');
    ylabel('East [m]');
    colormap(gca, turbo);
    colorbar;
else
    text(0.1, 0.5, 'FMM heuristic unavailable (fallback used).', 'FontSize', 12);
    axis off;
end

subplot(1, 2, 2);
warm_cost_plot = info.warmstart_cost;
if ~isfinite(warm_cost_plot)
    warm_cost_plot = NaN;
end
final_cost_plot = final_cost;
if ~isfinite(final_cost_plot)
    final_cost_plot = NaN;
end
bar([warm_cost_plot, final_cost_plot], 0.6, 'FaceColor', [0.18, 0.55, 0.82]);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Warm-start', 'Final'});
ylabel('Cost');
title('Warm-start vs Final Cost');
grid on;
yl = ylim;
text(1.05, yl(2) * 0.85, sprintf('Volume ratio \\rho = %.2f', info.volume_ratio), 'FontSize', 11);

%% 7) Console summary
fprintf('\nFMM preprocessing: %.2f s\n', info.fmm_time_s);
fprintf('BIT* iterations:   %d\n', info.iterations);
fprintf('Warm-start cost:   %.2f\n', info.warmstart_cost);
fprintf('Final cost:        %.2f\n', final_cost);
fprintf('Volume ratio rho:  %.2f\n', info.volume_ratio);
fprintf('FMM used:          %s\n', string(info.fmm_used));

%% ---------------- Local helpers ----------------
function [NE, h, risk] = choose_endpoint(threat, tm, target_NE, agl, visibility_threshold)
offset_vals = -700:140:700;
offsets = zeros(numel(offset_vals)^2, 2);
cursor = 1;
for dn = offset_vals
    for de = offset_vals
        offsets(cursor, :) = [dn, de];
        cursor = cursor + 1;
    end
end
best_cost = inf;
best_NE = target_NE(:);
best_h = tm.get_height(best_NE(1), best_NE(2));
best_risk = 1.0;

for i = 1:size(offsets, 1)
    cand_NE = target_NE(:) + offsets(i, :)';
    cand_NE(1) = min(max(cand_NE(1), tm.bounds(1)), tm.bounds(2));
    cand_NE(2) = min(max(cand_NE(2), tm.bounds(3)), tm.bounds(4));
    cand_h = tm.get_height(cand_NE(1), cand_NE(2));
    if ~isfinite(cand_h)
        continue;
    end
    cand_alt = cand_h + agl;
    cand_risk = threat.get_risk(cand_NE(1), cand_NE(2), cand_alt);
    score = 1000 * max(0, cand_risk - visibility_threshold) + norm(offsets(i, :));
    if score < best_cost
        best_cost = score;
        best_NE = cand_NE;
        best_h = cand_h;
        best_risk = cand_risk;
    end
end

NE = best_NE;
h = best_h;
risk = best_risk;
end
