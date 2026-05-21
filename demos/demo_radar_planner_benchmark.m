%% DEMO_RADAR_PLANNER_BENCHMARK
% Planner benchmark for AIT*, standard BIT*, and Informed RRT*.

clear; clc; close all;
clear fsm_bit_star_planner_direct rrt_star_radar ait_star_planner_direct;
rehash;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
setup_project_paths();

fprintf('=== Planner-Focused Radar Benchmark ===\n');

cfg = struct();
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 30;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 2500;
cfg.visibility_threshold = 0.5;
cfg.start_goal_agl = 80;
cfg.start_fraction_ne = [0.08; 0.10];
cfg.goal_fraction_ne = [0.92; 0.90];
cfg.threat_vert_res = 40;
cfg.use_parallel_threat = false;
cfg.radar_alt_offset = 20;
cfg.radar_range = 3000;
cfg.radar_power = 150e3;
cfg.radar_ne_frac = [0.22, 0.45, 0.80;
                     0.25, 0.55, 0.20];
cfg.radar_position_mode = 'manual_fraction';

%cfg.radar_ne_frac = [0.448240, 0.207284, 0.489079; ...
%                     0.625442, 0.495914, 0.293070];
cfg.num_radars = size(cfg.radar_ne_frac, 2);
cfg.radar_manual_ne = [];
cfg.radar_manual_fraction = cfg.radar_ne_frac;
cfg.radar_min_separation = 800;
cfg.radar_edge_margin_frac = 0.08;
cfg.common_stealth_alpha = 2.0;
cfg.strict_zero_risk_tol = 1e-6;
cfg.bit_strict_zero_risk = true;
cfg.bit_stealth_weight = inf;
cfg.preferred_agl = 40;
cfg.preferred_agl_cost_weight = 0.2;
cfg.enable_preferred_agl_cost = false;
cfg.common_time_budget_s = 45;
cfg.common_min_clearance = 20;
cfg.planner_tuning_policy = 'planner_best_effort';
cfg.ait_profile_runtime = false;
cfg.ait_use_parallel_reverse = true;
cfg.ait_use_parallel_edge_checks = true;
cfg.ait_parallel_threshold_neighbors = 96;
cfg.ait_parallel_threshold_edges = 48;
cfg.rrt_use_parallel_rewire = false;
cfg.rng_seed = 7;
cfg.seed_list = cfg.rng_seed;
cfg.planner_seeds = make_planner_seeds(cfg.seed_list(1));
cfg.num_trials = numel(cfg.seed_list);
cfg.scenarios = make_benchmark_scenarios();
cfg.radar_fraction_trials = [];
cfg.radar_ne_trials = [];
cfg.figure_click_use_contour = true;
cfg.figure_click_close_after = true;
cfg.edge_check_spacing = 10;
cfg.common_validation_spacing = 10;
cfg.report_export = true;
cfg.report_output_dir = fullfile(fileparts(this_dir), 'latex_report', 'Interim Report', 'v2', 'figures');
cfg.report_overlay_scenario_index = 1;
cfg.show_native_figures = true;
cfg.colors = struct( ...
    'ait_star', [0.90, 0.16, 0.54], ...
    'bit_star', [0.08, 0.08, 0.08], ...
    'rrt_star', [0.12, 0.92, 0.96], ...
    'radar', [1.00, 0.52, 0.10], ...
    'terrain', [0.67, 0.61, 0.54]);

rng(cfg.rng_seed, 'twister');

%% 1) Terrain
t_dem = tic;
[terrain_data, dem_source] = load_benchmark_terrain(cfg);
t_dem_s = toc(t_dem);

tm = terrain_map(terrain_data);
los = los_checker(tm);

if strcmpi(cfg.radar_position_mode, 'figure_click')
    cfg.scenarios = assign_scenario_radar_clicks(cfg.scenarios, tm, cfg);
elseif strcmpi(cfg.radar_position_mode, 'manual_fraction') && isfield(cfg, 'radar_ne_frac') && ~isempty(cfg.radar_ne_frac)
    cfg.scenarios = assign_scenario_radar_fractions(cfg.scenarios, cfg.radar_ne_frac);
end

terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
scenario_benchmark_cells = cell(1, numel(cfg.scenarios));
for scenario_idx = 1:numel(cfg.scenarios)
    scenario_cfg = apply_scenario_config(cfg, cfg.scenarios(scenario_idx));
    fprintf('\n=== Scenario %d/%d: %s ===\n', scenario_idx, numel(cfg.scenarios), scenario_cfg.scenario_name);

    radar_ne_trials = build_radar_ne_trials(scenario_cfg, tm);
    num_trials = size(radar_ne_trials, 3);
    trial_benchmark_cells = cell(1, num_trials);

    for trial_idx = 1:num_trials
        fprintf('\n=== Single Run %d/%d ===\n', trial_idx, num_trials);
        trial_cfg = scenario_cfg;
        trial_cfg.trial_seed = cfg.seed_list(trial_idx);
        trial_cfg.seed_index = trial_idx;
        trial_cfg.planner_seeds = make_planner_seeds(trial_cfg.trial_seed);
        trial_cfg.radar_ne = radar_ne_trials(:, :, trial_idx);
        trial_benchmark_cells{trial_idx} = run_single_trial_benchmark( ...
            trial_cfg, tm, los, dem_source, t_dem_s, terrain_min, terrain_max, diag_len);
    end

    trial_benchmarks = [trial_benchmark_cells{:}];
    scenario_benchmark = aggregate_trial_benchmarks(scenario_cfg, trial_benchmarks);
    scenario_benchmark.scenario_name = scenario_cfg.scenario_name;
    scenario_benchmark.scenario_slug = scenario_cfg.scenario_slug;
    scenario_benchmark.scenario_notes = getfield_with_default(scenario_cfg, 'scenario_notes', '');
    scenario_benchmark_cells{scenario_idx} = scenario_benchmark;
end

planner_benchmark = aggregate_scenario_benchmarks(cfg, [scenario_benchmark_cells{:}]);
planner_benchmark.output_files = export_report_assets(planner_benchmark, tm);
if getfield_with_default(cfg, 'show_native_figures', false)
    show_native_report_figures(planner_benchmark, tm);
end
print_scenario_benchmark_summary(planner_benchmark);

assignin('base', 'planner_benchmark', planner_benchmark);
assignin('base', 'planner_benchmark_scenarios', planner_benchmark.scenarios);
assignin('base', 'planner_benchmark_representative', planner_benchmark.representative_scenario.representative_trial);

%% ---------------- Local helpers ----------------
function planner_seeds = make_planner_seeds(base_seed)
planner_seeds = struct( ...
    'ait_star', base_seed + 101, ...
    'bit_star', base_seed + 202, ...
    'rrt_star', base_seed + 303);
end

function scenarios = assign_scenario_radar_clicks(scenarios, tm, cfg)
num_scenarios = numel(scenarios);
fig = figure('Name', 'Select Scenario Radar Positions', ...
    'Position', [120, 120, 1050, 800], 'Color', 'w');
imagesc(tm.N_vec, tm.E_vec, tm.Z);
set(gca, 'YDir', 'normal');
axis equal tight;
xlabel('North [m]');
ylabel('East [m]');
title(sprintf('Click %d radar points in order: Scenario 1, Scenario 2, Scenario 3', num_scenarios));
colormap(benchmark_terrain_colormap());
cb = colorbar;
cb.Label.String = 'Elevation [m]';
hold on;

if isfield(cfg, 'figure_click_use_contour') && cfg.figure_click_use_contour
    [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
    contour(N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
end

fprintf('Figure click mode: select %d radar points in order for Scenario 1, Scenario 2, Scenario 3.\n', num_scenarios);
fprintf('  Left click once per scenario, then press Enter.\n');
[x, y] = ginput(num_scenarios);
if numel(x) < num_scenarios
    error('Expected %d radar points, but only %d were selected.', num_scenarios, numel(x));
end

for s = 1:num_scenarios
    radar_ne = [x(s); y(s)];
    radar_ne(1) = min(max(radar_ne(1), tm.bounds(1)), tm.bounds(2));
    radar_ne(2) = min(max(radar_ne(2), tm.bounds(3)), tm.bounds(4));
    scenarios(s).radar_position_mode = 'manual_ne';
    scenarios(s).radar_manual_ne = radar_ne;
    scenarios(s).radar_ne_trials = reshape(radar_ne, 2, 1, 1);

    plot(radar_ne(1), radar_ne(2), 'd', 'MarkerSize', 12, ...
        'MarkerFaceColor', cfg.colors.radar, 'MarkerEdgeColor', 'k');
    text(radar_ne(1) + 40, radar_ne(2) + 40, sprintf('Scenario %d', s), ...
        'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10);
end

drawnow;
if isfield(cfg, 'figure_click_close_after') && cfg.figure_click_close_after
    close(fig);
end
end

function scenarios = assign_scenario_radar_fractions(scenarios, radar_ne_frac)
num_scenarios = numel(scenarios);
if size(radar_ne_frac, 1) ~= 2 || size(radar_ne_frac, 2) ~= num_scenarios
    error('cfg.radar_ne_frac must be [2 x %d], one [north_frac; east_frac] column per scenario.', num_scenarios);
end

for s = 1:num_scenarios
    frac = radar_ne_frac(:, s);
    frac = min(max(frac, 0), 1);
    scenarios(s).radar_position_mode = 'manual_fraction';
    scenarios(s).radar_fraction_ne = frac;
    scenarios(s).radar_manual_fraction = frac;
    scenarios(s).radar_ne_trials = [];
    scenarios(s).scenario_notes = sprintf('Radar location entered through radar\\_ne\\_frac for %s.', scenarios(s).scenario_name);
end
end

function scenarios = make_benchmark_scenarios()
scenarios = repmat(struct(), 1, 3);

scenarios(1).scenario_name = 'Scenario 1';
scenarios(1).scenario_slug = 'scenario_1';
scenarios(1).scenario_notes = 'Radar location for Scenario 1.';
scenarios(1).radar_fraction_ne = [0.22; 0.25];
scenarios(1).start_fraction_ne = [0.08; 0.10];
scenarios(1).goal_fraction_ne = [0.92; 0.90];
scenarios(1).radar_range = 3000;

scenarios(2).scenario_name = 'Scenario 2';
scenarios(2).scenario_slug = 'scenario_2';
scenarios(2).scenario_notes = 'Radar location for Scenario 2.';
scenarios(2).radar_fraction_ne = [0.46; 0.48];
scenarios(2).start_fraction_ne = [0.08; 0.12];
scenarios(2).goal_fraction_ne = [0.90; 0.88];
scenarios(2).radar_range = 2900;

scenarios(3).scenario_name = 'Scenario 3';
scenarios(3).scenario_slug = 'scenario_3';
scenarios(3).scenario_notes = 'Radar location for Scenario 3.';
scenarios(3).radar_fraction_ne = [0.68; 0.36];
scenarios(3).start_fraction_ne = [0.10; 0.16];
scenarios(3).goal_fraction_ne = [0.90; 0.82];
scenarios(3).radar_range = 3000;
end

function scenario_cfg = apply_scenario_config(cfg, scenario)
scenario_cfg = cfg;
fields = fieldnames(scenario);
for i = 1:numel(fields)
    scenario_cfg.(fields{i}) = scenario.(fields{i});
end
scenario_cfg.radar_manual_fraction = scenario_cfg.radar_fraction_ne;
scenario_cfg.num_radars = size(scenario_cfg.radar_fraction_ne, 2);
end

function [terrain_data, dem_source] = load_benchmark_terrain(cfg)
terrain_data = [];
dem_source = 'synthetic_fallback';

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
        dem_source = cfg.dem_file;
    end
catch ME
    warning('%s', sprintf('DEM load failed (%s). Using synthetic terrain fallback.', ME.message));
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
    terrain_data.type = 'synthetic_benchmark';
end
end

function planner_benchmark = run_single_trial_benchmark(cfg, tm, los, dem_source, dem_load_s, terrain_min, terrain_max, diag_len)
setup_tic = tic;
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);

threat_params = struct();
threat_params.resolution = [max(40, diag_len / 120), cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + 1, terrain_max + 350];

threat = threat_map(tm, los, 0.1, threat_params);
if isfield(cfg, 'radar_ne') && ~isempty(cfg.radar_ne)
    radar_ne = cfg.radar_ne;
else
    radar_ne = [tm.bounds(1) + cfg.radar_fraction_ne(1, :) * N_span; ...
                tm.bounds(3) + cfg.radar_fraction_ne(2, :) * E_span];
end

for r = 1:size(radar_ne, 2)
    rN = radar_ne(1, r);
    rE = radar_ne(2, r);
    rH = tm.get_height(rN, rE) + cfg.radar_alt_offset;
    radar_obj = radar_site([rN; rE; rH], sprintf('Radar-%d', r), 'tracking');
    radar_obj.R_max = cfg.radar_range;
    radar_obj.P_t = cfg.radar_power;
    threat.add_radar(radar_obj);
end

fprintf('Computing binary threat map...\n');
t_threat = tic;
threat.compute_map('binary', struct('use_parallel', cfg.use_parallel_threat, 'show_progress', false));
t_threat_s = toc(t_threat);
fprintf('Threat map compute time: %.2f s\n', t_threat_s);

radar_los_map = threat.risk_grid >= cfg.visibility_threshold;

[start_NE, start_h, start_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(1) + cfg.start_fraction_ne(1) * N_span; tm.bounds(3) + cfg.start_fraction_ne(2) * E_span], ...
    cfg.start_goal_agl, cfg.visibility_threshold);
[goal_NE, goal_h, goal_risk] = choose_endpoint(threat, tm, ...
    [tm.bounds(1) + cfg.goal_fraction_ne(1) * N_span; tm.bounds(3) + cfg.goal_fraction_ne(2) * E_span], ...
    cfg.start_goal_agl, cfg.visibility_threshold);

start_world = [start_NE; start_h + cfg.start_goal_agl];
goal_world = [goal_NE; goal_h + cfg.start_goal_agl];
start_ned = world_to_ned(start_world);
goal_ned = world_to_ned(goal_world);

fprintf('Start risk: %.2f | Goal risk: %.2f\n', start_risk, goal_risk);

shared_times = struct();
shared_times.dem_load_s = dem_load_s;
shared_times.threat_compute_s = t_threat_s;
shared_times.total_setup_s = toc(setup_tic);

planner_benchmark = struct();
planner_benchmark.cfg = cfg;
planner_benchmark.trial_seed = getfield_with_default(cfg, 'trial_seed', NaN);
planner_benchmark.seed_index = getfield_with_default(cfg, 'seed_index', NaN);
planner_benchmark.environment = build_environment_summary( ...
    cfg, tm, threat, dem_source, radar_los_map, radar_ne, ...
    start_world, goal_world, start_ned, goal_ned, start_risk, goal_risk);
planner_benchmark.conditions = build_common_benchmark_conditions(cfg, diag_len);
planner_benchmark.shared_times = shared_times;
planner_benchmark.planners = struct();
planner_benchmark.runtime = struct('threat', threat);

fprintf('\n[1/3] Running AIT*...\n');
rng(cfg.planner_seeds.ait_star, 'twister');
planner_benchmark.planners.ait_star = run_ait_star_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len);

fprintf('\n[2/3] Running BIT* (Standard)...\n');
rng(cfg.planner_seeds.bit_star, 'twister');
planner_benchmark.planners.bit_star = run_bit_star_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len);

fprintf('\n[3/3] Running Informed RRT*...\n');
rng(cfg.planner_seeds.rrt_star, 'twister');
planner_benchmark.planners.rrt_star = run_rrt_star_benchmark(cfg, tm, threat, start_ned, goal_ned, terrain_min, terrain_max, diag_len);

planner_benchmark.comparison = build_comparison_struct(planner_benchmark);
end

function radar_fraction_trials = build_radar_fraction_trials(cfg, n_radars)
if isfield(cfg, 'radar_fraction_trials') && ~isempty(cfg.radar_fraction_trials)
    radar_fraction_trials = cfg.radar_fraction_trials;
    if ndims(radar_fraction_trials) == 2
        radar_fraction_trials = reshape(radar_fraction_trials, size(radar_fraction_trials, 1), size(radar_fraction_trials, 2), 1);
    end
    return;
end

num_trials = max(1, round(getfield_with_default(cfg, 'num_trials', 1)));
margin = getfield_with_default(cfg, 'radar_trial_margin_frac', 0.12);
margin = min(max(margin, 0.0), 0.45);
radar_fraction_trials = zeros(2, n_radars, num_trials);
radar_fraction_trials(:, :, 1) = cfg.radar_fraction_ne;
if num_trials == 1
    return;
end

rng(cfg.rng_seed + 700, 'twister');
for k = 2:num_trials
    radar_fraction_trials(1, :, k) = margin + (1 - 2 * margin) * rand(1, n_radars);
    radar_fraction_trials(2, :, k) = margin + (1 - 2 * margin) * rand(1, n_radars);
end
end

function radar_ne_trials = build_radar_ne_trials(cfg, tm)
n_radars = max(1, round(getfield_with_default(cfg, 'num_radars', size(cfg.radar_fraction_ne, 2))));
num_trials = max(1, round(getfield_with_default(cfg, 'num_trials', 1)));

if isfield(cfg, 'radar_ne_trials') && ~isempty(cfg.radar_ne_trials)
    radar_ne_trials = cfg.radar_ne_trials;
    if ndims(radar_ne_trials) == 2
        radar_ne_trials = reshape(radar_ne_trials, size(radar_ne_trials, 1), size(radar_ne_trials, 2), 1);
    end
    return;
end

mode_name = lower(getfield_with_default(cfg, 'radar_position_mode', 'trial_fraction'));
switch mode_name
    case 'trial_fraction'
        radar_fraction_trials = build_radar_fraction_trials(cfg, n_radars);
        num_trials = size(radar_fraction_trials, 3);
        radar_ne_trials = zeros(2, n_radars, num_trials);
        for k = 1:num_trials
            frac = radar_fraction_trials(:, :, k);
            radar_ne_trials(1, :, k) = tm.bounds(1) + frac(1, :) * (tm.bounds(2) - tm.bounds(1));
            radar_ne_trials(2, :, k) = tm.bounds(3) + frac(2, :) * (tm.bounds(4) - tm.bounds(3));
        end
    otherwise
        radar_ne = select_radar_positions(tm, cfg, n_radars);
        radar_ne_trials = repmat(radar_ne, 1, 1, num_trials);
end
end

function radar_ne = select_radar_positions(tm, cfg, n_radars)
n_radars = max(1, round(n_radars));
mode_name = lower(getfield_with_default(cfg, 'radar_position_mode', 'trial_fraction'));

switch mode_name
    case 'manual_ne'
        radar_ne = cfg.radar_manual_ne;
        if size(radar_ne, 1) ~= 2
            error('cfg.radar_manual_ne must be [2 x M] as [N;E].');
        end

    case 'manual_fraction'
        frac = cfg.radar_manual_fraction;
        if size(frac, 1) ~= 2
            error('cfg.radar_manual_fraction must be [2 x M] in [0,1].');
        end
        frac(1, :) = min(max(frac(1, :), 0), 1);
        frac(2, :) = min(max(frac(2, :), 0), 1);
        radar_ne = [tm.bounds(1) + frac(1, :) * (tm.bounds(2) - tm.bounds(1)); ...
                    tm.bounds(3) + frac(2, :) * (tm.bounds(4) - tm.bounds(3))];

    case 'figure_click'
        fig = figure('Name', 'Select Radar Positions', 'Position', [120, 120, 1000, 760], 'Color', 'w');
        imagesc(tm.N_vec, tm.E_vec, tm.Z);
        set(gca, 'YDir', 'normal');
        axis equal tight;
        xlabel('North [m]');
        ylabel('East [m]');
        scenario_name = getfield_with_default(cfg, 'scenario_name', 'Scenario');
        title(sprintf('%s: click %d radar position(s), then press Enter', scenario_name, n_radars));
        colormap(turbo);
        colorbar;
        hold on;
        if isfield(cfg, 'figure_click_use_contour') && cfg.figure_click_use_contour
            [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
            contour(N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
        end

        fprintf('Figure click mode: %s -> select %d radar position(s) on the terrain map.\n', ...
            scenario_name, n_radars);
        fprintf('  Left click to add points, Enter to finish.\n');
        [x, y] = ginput(n_radars);
        if isempty(x)
            error('No radar point selected in figure_click mode.');
        end

        radar_ne = [x(:)'; y(:)'];
        if size(radar_ne, 2) < n_radars
            warning('Selected %d/%d points. Repeating last point for remaining radars.', ...
                size(radar_ne, 2), n_radars);
            radar_ne = [radar_ne, repmat(radar_ne(:, end), 1, n_radars - size(radar_ne, 2))];
        end

        plot(radar_ne(1, :), radar_ne(2, :), 'd', 'MarkerSize', 12, ...
            'MarkerFaceColor', cfg.colors.radar, 'MarkerEdgeColor', 'k');
        for idx = 1:size(radar_ne, 2)
            text(radar_ne(1, idx) + 40, radar_ne(2, idx) + 40, sprintf('Radar %d', idx), ...
                'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10);
        end
        drawnow;
        if isfield(cfg, 'figure_click_close_after') && cfg.figure_click_close_after
            close(fig);
        end

    case 'terrain_peaks'
        margin_n = cfg.radar_edge_margin_frac * (tm.bounds(2) - tm.bounds(1));
        margin_e = cfg.radar_edge_margin_frac * (tm.bounds(4) - tm.bounds(3));
        N_min = tm.bounds(1) + margin_n;
        N_max = tm.bounds(2) - margin_n;
        E_min = tm.bounds(3) + margin_e;
        E_max = tm.bounds(4) - margin_e;

        [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
        N_flat = N_grid(:);
        E_flat = E_grid(:);
        Z_flat = tm.Z(:);

        valid = N_flat >= N_min & N_flat <= N_max & E_flat >= E_min & E_flat <= E_max;
        N_flat = N_flat(valid);
        E_flat = E_flat(valid);
        Z_flat = Z_flat(valid);

        [~, order] = sort(Z_flat, 'descend');
        selected = zeros(2, 0);

        for idx = order'
            cand = [N_flat(idx); E_flat(idx)];
            if isempty(selected)
                selected(:, end+1) = cand; %#ok<AGROW>
            else
                d = vecnorm(selected - cand, 2, 1);
                if all(d >= cfg.radar_min_separation)
                    selected(:, end+1) = cand; %#ok<AGROW>
                end
            end
            if size(selected, 2) >= n_radars
                break;
            end
        end

        if isempty(selected)
            selected = [mean(tm.bounds(1:2)); mean(tm.bounds(3:4))];
        end
        radar_ne = selected;

    otherwise
        error('Radar position mode "%s" does not support direct manual selection here.', cfg.radar_position_mode);
end

if size(radar_ne, 2) < n_radars
    radar_ne = repmat(radar_ne(:, end), 1, n_radars);
elseif size(radar_ne, 2) > n_radars
    radar_ne = radar_ne(:, 1:n_radars);
end

radar_ne(1, :) = min(max(radar_ne(1, :), tm.bounds(1)), tm.bounds(2));
radar_ne(2, :) = min(max(radar_ne(2, :), tm.bounds(3)), tm.bounds(4));
end

function planner_benchmark = aggregate_trial_benchmarks(cfg, trial_benchmarks)
comparison = build_aggregate_comparison_struct(cfg, trial_benchmarks);
rep_idx = choose_representative_trial_index(trial_benchmarks, comparison);

planner_benchmark = struct();
planner_benchmark.cfg = cfg;
planner_benchmark.num_trials = numel(trial_benchmarks);
planner_benchmark.seed_list = arrayfun(@(t) getfield_with_default(t, 'trial_seed', NaN), trial_benchmarks);
planner_benchmark.trials = trial_benchmarks;
planner_benchmark.representative_trial_index = rep_idx;
planner_benchmark.representative_trial = trial_benchmarks(rep_idx);
planner_benchmark.conditions = trial_benchmarks(1).conditions;
planner_benchmark.shared_times = struct( ...
    'dem_load_s', trial_benchmarks(1).shared_times.dem_load_s, ...
    'threat_compute_mean_s', mean(arrayfun(@(t) t.shared_times.threat_compute_s, trial_benchmarks), 'omitnan'), ...
    'threat_compute_std_s', std(arrayfun(@(t) t.shared_times.threat_compute_s, trial_benchmarks), 0, 'omitnan'), ...
    'total_setup_mean_s', mean(arrayfun(@(t) t.shared_times.total_setup_s, trial_benchmarks), 'omitnan'), ...
    'total_setup_std_s', std(arrayfun(@(t) t.shared_times.total_setup_s, trial_benchmarks), 0, 'omitnan'));
planner_benchmark.comparison = comparison;
planner_benchmark.output_files = struct();
end

function planner_benchmark = aggregate_scenario_benchmarks(cfg, scenario_benchmarks)
planner_benchmark = struct();
planner_benchmark.cfg = cfg;
planner_benchmark.num_scenarios = numel(scenario_benchmarks);
planner_benchmark.scenarios = scenario_benchmarks;
planner_benchmark.scenario_names = arrayfun(@(s) s.scenario_name, scenario_benchmarks, 'UniformOutput', false);
planner_benchmark.scenario_slugs = arrayfun(@(s) s.scenario_slug, scenario_benchmarks, 'UniformOutput', false);
rep_idx = round(getfield_with_default(cfg, 'report_overlay_scenario_index', 1));
rep_idx = min(max(rep_idx, 1), numel(scenario_benchmarks));
planner_benchmark.representative_scenario_index = rep_idx;
planner_benchmark.representative_scenario = scenario_benchmarks(rep_idx);
planner_benchmark.representative_trial = scenario_benchmarks(rep_idx).representative_trial;
planner_benchmark.output_files = struct();
end

function comparison = build_aggregate_comparison_struct(cfg, trial_benchmarks)
base_cmp = trial_benchmarks(1).comparison;
keys = base_cmp.keys;
labels = base_cmp.labels;
colors = base_cmp.colors;
n = numel(keys);
n_trials = numel(trial_benchmarks);
trial_seed = arrayfun(@(t) getfield_with_default(t, 'trial_seed', NaN), trial_benchmarks);

success_rate = zeros(1, n);
metric_sample_count = zeros(1, n);
total_time_mean_s = nan(1, n);
total_time_std_s = nan(1, n);
search_time_mean_s = nan(1, n);
search_time_std_s = nan(1, n);
path_length_mean_m = nan(1, n);
path_length_std_m = nan(1, n);
risk_integrated_mean = nan(1, n);
risk_integrated_std = nan(1, n);
per_seed = struct();

for i = 1:n
    total_time_vals = nan(1, n_trials);
    search_time_vals = nan(1, n_trials);
    length_vals = nan(1, n_trials);
    risk_integrated_vals = nan(1, n_trials);
    max_risk_vals = nan(1, n_trials);
    success_vals = false(1, n_trials);
    for t = 1:n_trials
        item = trial_benchmarks(t).planners.(keys{i});
        total_time_vals(t) = item.total_time_s;
        search_time_vals(t) = item.search_time_s;
        length_vals(t) = item.path_length_m;
        risk_integrated_vals(t) = item.risk_integrated;
        max_risk_vals(t) = item.dense_max_risk;
        success_vals(t) = logical(item.success);
    end

    good_metric_mask = success_vals & isfinite(total_time_vals) & isfinite(search_time_vals) & ...
        isfinite(length_vals) & isfinite(risk_integrated_vals);

    success_rate(i) = mean(success_vals);
    metric_sample_count(i) = sum(good_metric_mask);
    total_time_mean_s(i) = mean(total_time_vals(good_metric_mask), 'omitnan');
    total_time_std_s(i) = std(total_time_vals(good_metric_mask), 0, 'omitnan');
    search_time_mean_s(i) = mean(search_time_vals(good_metric_mask), 'omitnan');
    search_time_std_s(i) = std(search_time_vals(good_metric_mask), 0, 'omitnan');
    path_length_mean_m(i) = mean(length_vals(good_metric_mask), 'omitnan');
    path_length_std_m(i) = std(length_vals(good_metric_mask), 0, 'omitnan');
    risk_integrated_mean(i) = mean(risk_integrated_vals(good_metric_mask), 'omitnan');
    risk_integrated_std(i) = std(risk_integrated_vals(good_metric_mask), 0, 'omitnan');

    per_seed.(keys{i}) = struct( ...
        'trial_seed', trial_seed, ...
        'success', success_vals, ...
        'total_time_s', total_time_vals, ...
        'search_time_s', search_time_vals, ...
        'path_length_m', length_vals, ...
        'risk_integrated', risk_integrated_vals, ...
        'dense_max_risk', max_risk_vals);
end

comparison = struct();
comparison.keys = keys;
comparison.labels = labels;
comparison.colors = colors;
comparison.trial_seed = trial_seed;
comparison.per_seed = per_seed;
comparison.success_rate = success_rate;
comparison.metric_sample_count = metric_sample_count;
comparison.total_time_mean_s = total_time_mean_s;
comparison.total_time_std_s = total_time_std_s;
comparison.search_time_mean_s = search_time_mean_s;
comparison.search_time_std_s = search_time_std_s;
comparison.path_length_mean_m = path_length_mean_m;
comparison.path_length_std_m = path_length_std_m;
comparison.risk_integrated_mean = risk_integrated_mean;
comparison.risk_integrated_std = risk_integrated_std;
end

function representative_idx = choose_representative_trial_index(trial_benchmarks, comparison)
num_trials = numel(trial_benchmarks);
trial_scores = inf(1, num_trials);

for t = 1:num_trials
    trial_cmp = trial_benchmarks(t).comparison;
    if ~all(trial_cmp.success)
        continue;
    end

    score = 0.0;
    for i = 1:numel(trial_cmp.keys)
        score = score + normalized_metric_delta(trial_cmp.total_time_s(i), comparison.total_time_mean_s(i));
        score = score + normalized_metric_delta(trial_cmp.path_length_m(i), comparison.path_length_mean_m(i));
        score = score + normalized_metric_delta(trial_cmp.risk_integrated(i), comparison.risk_integrated_mean(i));
    end
    trial_scores(t) = score;
end

[best_score, representative_idx] = min(trial_scores);
if isfinite(best_score)
    return;
end

success_counts = arrayfun(@(t) sum(t.comparison.success), trial_benchmarks);
[~, representative_idx] = max(success_counts);
end

function delta = normalized_metric_delta(value, mean_value)
delta = 1.0;
if isfinite(value) && isfinite(mean_value)
    delta = abs(value - mean_value) / max(abs(mean_value), 1e-9);
end
end

function print_scenario_benchmark_summary(planner_benchmark)
fprintf('\n=== Scenario Benchmark Summary ===\n');
for s = 1:planner_benchmark.num_scenarios
    scenario = planner_benchmark.scenarios(s);
    cmp = scenario.comparison;
    fprintf('\n--- %s ---\n', scenario.scenario_name);
    if ~isempty(getfield_with_default(scenario, 'scenario_notes', ''))
        fprintf('%s\n', scenario.scenario_notes);
    end
    print_shared_conditions_summary(scenario.conditions);
    fprintf(['Seeds: %d | DEM load: %.2f s | Threat compute: %.2f +/- %.2f s | ', ...
             'Representative seed: %d\n'], ...
        scenario.num_trials, ...
        scenario.shared_times.dem_load_s, ...
        scenario.shared_times.threat_compute_mean_s, ...
        scenario.shared_times.threat_compute_std_s, ...
        getfield_with_default(scenario.representative_trial, 'trial_seed', NaN));

    for i = 1:numel(cmp.keys)
        fprintf(['%-15s | success_rate=%.2f | total=%.2f +/- %.2f s | search=%.2f +/- %.2f s | ', ...
                 'length=%.1f +/- %.1f m | risk=%.4f +/- %.4f\n'], ...
            cmp.labels{i}, cmp.success_rate(i), ...
            cmp.total_time_mean_s(i), cmp.total_time_std_s(i), ...
            cmp.search_time_mean_s(i), cmp.search_time_std_s(i), ...
            cmp.path_length_mean_m(i), cmp.path_length_std_m(i), ...
            cmp.risk_integrated_mean(i), cmp.risk_integrated_std(i));
    end
end

if isfield(planner_benchmark, 'output_files') && isstruct(planner_benchmark.output_files)
    fprintf('\nScenario data:  %s\n', getfield_with_default(planner_benchmark.output_files, 'scenario_csv', ''));
    fprintf('Path figure:    %s\n', getfield_with_default(planner_benchmark.output_files, 'overlay_figure', ''));
    fprintf('Summary file:   %s\n', getfield_with_default(planner_benchmark.output_files, 'summary_text', ''));
end
fprintf('planner_benchmark and planner_benchmark_scenarios exported to workspace.\n');
end

function env_summary = build_environment_summary(cfg, tm, threat, dem_source, radar_los_map, radar_ne, start_world, goal_world, start_ned, goal_ned, start_risk, goal_risk)
radar_sites = repmat(struct('name', '', 'position', [], 'range_m', NaN), 1, numel(threat.radars));
for i = 1:numel(threat.radars)
    rr = threat.radars{i};
    radar_sites(i).name = rr.name;
    radar_sites(i).position = rr.position;
    radar_sites(i).range_m = rr.R_max;
end

env_summary = struct();
env_summary.dem_source = dem_source;
env_summary.terrain_type = getfield_with_default(tm, 'type', '');
env_summary.terrain_bounds = tm.bounds;
env_summary.terrain_resolution = getfield_with_default(tm, 'resolution', mean(diff(tm.N_vec)));
env_summary.threat_resolution = threat.resolution;
env_summary.threat_alt_range = [threat.alt_vec(1), threat.alt_vec(end)];
env_summary.threat_grid_size = size(threat.risk_grid);
env_summary.visibility_threshold = cfg.visibility_threshold;
env_summary.radar_fraction_ne = cfg.radar_fraction_ne;
env_summary.radar_positions_ne = radar_ne;
env_summary.radar_sites = radar_sites;
env_summary.radar_los_map_size = size(radar_los_map);
env_summary.start_world = start_world;
env_summary.goal_world = goal_world;
env_summary.start_ned = start_ned;
env_summary.goal_ned = goal_ned;
env_summary.start_risk = start_risk;
env_summary.goal_risk = goal_risk;
end

function conditions = build_common_benchmark_conditions(cfg, diag_len)
conditions = struct();
conditions.strict_zero_risk = logical(cfg.bit_strict_zero_risk);
conditions.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
conditions.visibility_threshold = cfg.visibility_threshold;
conditions.planner_tuning_policy = getfield_with_default(cfg, 'planner_tuning_policy', 'planner_best_effort');
conditions.preferred_agl_enabled = logical(getfield_with_default(cfg, 'enable_preferred_agl_cost', false));
conditions.preferred_agl = cfg.preferred_agl;
conditions.preferred_agl_cost_weight = 0.0;
conditions.objective = 'length_only';
if conditions.preferred_agl_enabled
    conditions.preferred_agl_cost_weight = cfg.preferred_agl_cost_weight;
    conditions.objective = 'length_plus_agl_preference';
end
conditions.time_budget_s = getfield_with_default(cfg, 'common_time_budget_s', 45);
conditions.min_clearance_m = getfield_with_default(cfg, 'common_min_clearance', 20);
conditions.edge_check_spacing_m = getfield_with_default(cfg, 'edge_check_spacing', 10);
conditions.validation_spacing_m = getfield_with_default(cfg, 'common_validation_spacing', 10);
conditions.diag_len_m = diag_len;
conditions.ait_profile_runtime = logical(getfield_with_default(cfg, 'ait_profile_runtime', false));
conditions.ait_use_parallel_reverse = logical(getfield_with_default(cfg, 'ait_use_parallel_reverse', false));
conditions.ait_use_parallel_edge_checks = logical(getfield_with_default(cfg, 'ait_use_parallel_edge_checks', false));
conditions.planner_seeds = getfield_with_default(cfg, 'planner_seeds', struct());
end

function result = run_ait_star_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len)
[Nq, Eq] = meshgrid(threat.N_vec, threat.E_vec);
Zq = tm.get_height(Nq, Eq);
common = build_common_benchmark_conditions(cfg, diag_len);

env = struct();
env.terrain = struct( ...
    'Z', Zq, ...
    'N_vec', threat.N_vec, ...
    'E_vec', threat.E_vec, ...
    'x_vec', threat.N_vec, ...
    'y_vec', threat.E_vec, ...
    'dx', mean(diff(threat.N_vec)), ...
    'dy', mean(diff(threat.E_vec)), ...
    'alt_vec', threat.alt_vec);
env.radar_los_map = radar_los_map;
env.N_vec = threat.N_vec;
env.E_vec = threat.E_vec;
env.alt_vec = threat.alt_vec;
env.visibility_threshold = cfg.visibility_threshold;
env.bounds_world = [threat.N_vec(1), threat.N_vec(end), threat.E_vec(1), threat.E_vec(end), ...
                    threat.alt_vec(1), threat.alt_vec(end)];

planner_params = struct();
planner_params.max_batches = 120;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.12;
planner_params.stealth_weight = cfg.bit_stealth_weight;
planner_params.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
planner_params.gamma = common.preferred_agl_cost_weight;
planner_params.preferred_agl = cfg.preferred_agl;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.edge_check_samples = compute_edge_check_samples(planner_params.max_connection_radius, cfg.edge_check_spacing);
planner_params.eta_radius = 1.8;
planner_params.max_time = common.time_budget_s;
planner_params.min_clearance = common.min_clearance_m;
planner_params.reverse_k_neighbors = 24;
planner_params.forward_k_neighbors = 16;
planner_params.goal_sample_fraction = 0.24;
planner_params.line_sample_fraction = 0.32;
planner_params.goal_probe_count = 14;
planner_params.profile_ait = common.ait_profile_runtime;
planner_params.neighbor_precompute_mode = 'batch';
planner_params.use_parallel_reverse = common.ait_use_parallel_reverse;
planner_params.use_parallel_edge_checks = common.ait_use_parallel_edge_checks;
planner_params.parallel_threshold_neighbors = getfield_with_default(cfg, 'ait_parallel_threshold_neighbors', 96);
planner_params.parallel_threshold_edges = getfield_with_default(cfg, 'ait_parallel_threshold_edges', 48);
planner_params.reverse_batch_repair_limit = 256;
planner_params.forward_expand_chunk = 32;
planner_params.kd_rebuild_interval = 1;
planner_params.reverse_debug = false;

result = init_planner_result('AIT*');
result.strict_zero_risk = cfg.bit_strict_zero_risk;
t_total = tic;
try
    [path_world, native_cost, info] = ait_star_planner_direct(start_world, goal_world, env, planner_params);
    result = finalize_world_path_result(result, path_world, threat, cfg.common_stealth_alpha, cfg.visibility_threshold, cfg.common_validation_spacing);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = 0.0;
    result.heuristic_time_s = read_info_field(info, 'reverse_time_s');
    result.search_time_s = read_info_field(info, 'planning_time');
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = read_info_text(info, 'stop_reason');
    result.native_path_cost = finite_or_nan(native_cost);
    result.samples_generated = read_info_field(info, 'samples_generated');
    result.samples_accepted = read_info_field(info, 'samples_accepted');
    result.connection_radius_last = read_info_field(info, 'connection_radius_last');
    result.first_solution_batch = read_info_field(info, 'first_solution_batch');
    result.native_summary = extract_native_summary(info, ...
        {'risk_limit', 'reverse_repairs', 'invalid_edge_count', 'reverse_queue_pops', ...
         'forward_edge_checks', 'adaptive_heuristic_ready_count', 'reverse_rebuild_time_s', ...
         'reverse_repair_time_s', 'forward_queue_build_time_s', 'goal_probe_time_s', ...
         'neighbor_precompute_time_s', 'edge_validation_time_s', 'parallel_used_reverse', ...
         'parallel_used_edge_checks'});
    if ~result.success
        result.error_message = derive_failure_message(info, result);
    end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end


function result = run_bit_star_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len)
[Nq_fsm, Eq_fsm] = meshgrid(threat.N_vec, threat.E_vec);
Z_fsm = tm.get_height(Nq_fsm, Eq_fsm);
common = build_common_benchmark_conditions(cfg, diag_len);

env = struct();
env.terrain = struct( 'Z', Z_fsm, 'N_vec', threat.N_vec, 'E_vec', threat.E_vec, 'x_vec', threat.N_vec, 'y_vec', threat.E_vec, 'dx', mean(diff(threat.N_vec)), 'dy', mean(diff(threat.E_vec)), 'alt_vec', threat.alt_vec);
env.radar_los_map = radar_los_map;
env.N_vec = threat.N_vec;
env.E_vec = threat.E_vec;
env.alt_vec = threat.alt_vec;
env.visibility_threshold = cfg.visibility_threshold;
env.bounds_world = [threat.N_vec(1), threat.N_vec(end), threat.E_vec(1), threat.E_vec(end), threat.alt_vec(1), threat.alt_vec(end)];

planner_params = struct();
planner_params.max_batches = 200;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.12;
planner_params.stealth_weight = cfg.bit_stealth_weight;
planner_params.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
planner_params.gamma = common.preferred_agl_cost_weight;
planner_params.preferred_agl = cfg.preferred_agl;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.edge_check_samples = compute_edge_check_samples(planner_params.max_connection_radius, cfg.edge_check_spacing);
planner_params.eta_radius = 1.8;
planner_params.max_time = common.time_budget_s;
planner_params.min_clearance = common.min_clearance_m;
planner_params.force_euclidean_heuristic = true;

result = init_planner_result('BIT*');
result.strict_zero_risk = cfg.bit_strict_zero_risk;
t_total = tic;
try
    [path_world, native_cost, info] = fsm_bit_star_planner_direct(start_world, goal_world, env, planner_params);
    result = finalize_world_path_result(result, path_world, threat, cfg.common_stealth_alpha, cfg.visibility_threshold, cfg.common_validation_spacing);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = 0;
    result.heuristic_time_s = 0;
    result.search_time_s = finite_or_nan(info.planning_time);
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = read_info_text(info, 'stop_reason');
    result.native_path_cost = finite_or_nan(native_cost);
    if ~result.success, result.error_message = derive_failure_message(info, result); end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = run_fsm_heuristic_bit_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len)
[Nq_fsm, Eq_fsm] = meshgrid(threat.N_vec, threat.E_vec);
Z_fsm = tm.get_height(Nq_fsm, Eq_fsm);
common = build_common_benchmark_conditions(cfg, diag_len);

env = struct();
env.terrain = struct( 'Z', Z_fsm, 'N_vec', threat.N_vec, 'E_vec', threat.E_vec, 'x_vec', threat.N_vec, 'y_vec', threat.E_vec, 'dx', mean(diff(threat.N_vec)), 'dy', mean(diff(threat.E_vec)), 'alt_vec', threat.alt_vec);
env.radar_los_map = radar_los_map;
env.N_vec = threat.N_vec;
env.E_vec = threat.E_vec;
env.alt_vec = threat.alt_vec;
env.visibility_threshold = cfg.visibility_threshold;
env.bounds_world = [threat.N_vec(1), threat.N_vec(end), threat.E_vec(1), threat.E_vec(end), threat.alt_vec(1), threat.alt_vec(end)];

planner_params = struct();
planner_params.max_batches = 120;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.12;
planner_params.stealth_weight = cfg.bit_stealth_weight;
planner_params.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
planner_params.gamma = common.preferred_agl_cost_weight;
planner_params.preferred_agl = cfg.preferred_agl;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.edge_check_samples = compute_edge_check_samples(planner_params.max_connection_radius, cfg.edge_check_spacing);
planner_params.eta_radius = 1.8;
planner_params.max_time = common.time_budget_s;
planner_params.min_clearance = common.min_clearance_m;
planner_params.w_mask = 1.0;
planner_params.w_alt = 0.08;
planner_params.h_opt = cfg.preferred_agl;
planner_params.max_sweeps = cfg.fsm_max_sweeps;
planner_params.tol = cfg.fsm_tol;

result = init_planner_result('FSM-Heur BIT*');
result.strict_zero_risk = cfg.bit_strict_zero_risk;
t_total = tic;
try
    [path_world, native_cost, info] = fsm_bit_star_planner_direct(start_world, goal_world, env, planner_params);
    result = finalize_world_path_result(result, path_world, threat, cfg.common_stealth_alpha, cfg.visibility_threshold, cfg.common_validation_spacing);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = finite_or_nan(read_info_field(info, 'fsm_time_s'));
    result.heuristic_time_s = finite_or_nan(read_info_field(info, 'fsm_time_s'));
    result.search_time_s = finite_or_nan(read_info_field(info, 'planning_time'));
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = read_info_text(info, 'stop_reason');
    result.native_path_cost = finite_or_nan(native_cost);
    if ~result.success, result.error_message = derive_failure_message(info, result); end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = run_fsm_corridor_bit_benchmark(cfg, tm, threat, radar_los_map, start_world, goal_world, diag_len)
result = init_planner_result('FSM-Corr BIT*');
result.strict_zero_risk = cfg.bit_strict_zero_risk;
common = build_common_benchmark_conditions(cfg, diag_len);

corr_params = struct();
corr_params.N_vec = threat.N_vec;
corr_params.E_vec = threat.E_vec;
corr_params.alt_vec = threat.alt_vec;
corr_params.max_sweeps = cfg.corridor_max_sweeps;
corr_params.lambda = cfg.corridor_lambda;
corr_params.r_min = cfg.corridor_r_min;
corr_params.r_max = cfg.corridor_r_max;

planner_params = struct();
planner_params.max_batches = 120;
planner_params.batch_size = 260;
planner_params.max_nodes = 18000;
planner_params.eps_global = 0.08;
planner_params.stealth_weight = cfg.bit_stealth_weight;
planner_params.strict_zero_risk_tol = cfg.strict_zero_risk_tol;
planner_params.visibility_query_mode = 'binary';
planner_params.gamma = common.preferred_agl_cost_weight;
planner_params.preferred_agl = cfg.preferred_agl;
planner_params.terrain_map = tm;
planner_params.post_solution_batches = 12;
planner_params.goal_radius = max(40, 0.02 * diag_len);
planner_params.min_connection_radius = max(30, 0.012 * diag_len);
planner_params.max_connection_radius = max(120, 0.08 * diag_len);
planner_params.edge_check_samples = compute_edge_check_samples(planner_params.max_connection_radius, cfg.edge_check_spacing);
planner_params.eta_radius = 1.8;
planner_params.max_time = common.time_budget_s;
planner_params.min_clearance = common.min_clearance_m;

t_total = tic;
try
    t_corr = tic;
    corridor = compute_stealth_corridor(threat.risk_grid, start_world, goal_world, cfg.visibility_threshold, cfg.corridor_alpha, corr_params);
    result.corridor_time_s = toc(t_corr);

    [path_world, info] = bit_star_planner(start_world, goal_world, corridor, planner_params);
    result = finalize_world_path_result(result, path_world, threat, cfg.common_stealth_alpha, cfg.visibility_threshold, cfg.common_validation_spacing);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = finite_or_nan(result.corridor_time_s);
    result.search_time_s = finite_or_nan(read_info_field(info, 'planning_time'));
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.native_path_cost = read_info_field(info, 'path_cost');
    if ~result.success, result.error_message = derive_failure_message(info, result); end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = run_rrt_star_benchmark(cfg, tm, threat, start_ned, goal_ned, terrain_min, terrain_max, diag_len)
common = build_common_benchmark_conditions(cfg, diag_len);
planner_params = struct();
planner_params.max_iter = cfg.dem_crop_half_size * 6;
planner_params.base_step_size = max(40, min(75, diag_len / 85));
planner_params.max_step_size = planner_params.base_step_size * 3.5;
planner_params.rewire_radius = planner_params.base_step_size * 3;
planner_params.rewire_mode = 'fixed';
planner_params.rewire_k_const = 3.2 * exp(1) * (1 + 1 / 3);
planner_params.rewire_k_min = 50;
planner_params.rewire_k_max = 350;
planner_params.goal_bias = 0.15;
planner_params.goal_bias_after_goal = 0.04;
planner_params.min_clearance = common.min_clearance_m;
planner_params.max_flight_alt = 500;
planner_params.alpha = 1.0;
planner_params.beta = 0.0;
planner_params.gamma = common.preferred_agl_cost_weight;
planner_params.preferred_agl = cfg.preferred_agl;
planner_params.edge_sample_spacing = common.edge_check_spacing_m;
planner_params.shadow_bias = 0.85;
planner_params.use_parallel_rewire = logical(getfield_with_default(cfg, 'rrt_use_parallel_rewire', false));
planner_params.animate = false;
planner_params.plot_interval = 500;
planner_params.height_query_mode = 'lut';
planner_params.radar_hard_constraint = true;
planner_params.radar_visibility_threshold = cfg.visibility_threshold;
planner_params.max_planning_time = common.time_budget_s;

min_alt = terrain_min + planner_params.min_clearance;
max_alt = terrain_max + planner_params.max_flight_alt;
planner_params.bounds = [tm.bounds(1), tm.bounds(2), tm.bounds(3), tm.bounds(4), -max_alt, -min_alt];

result = init_planner_result('Informed RRT*');
t_total = tic;
try
    [path_ned, info] = rrt_star_radar(start_ned, goal_ned, tm, threat, planner_params);
    result = finalize_ned_path_result(result, path_ned, threat, cfg.common_stealth_alpha, cfg.visibility_threshold, cfg.common_validation_spacing);
    result.success = logical(info.success);
    result.error_message = '';
    result.preprocess_time_s = 0.0;
    result.heuristic_time_s = 0.0;
    result.search_time_s = read_info_field(info, 'planning_time');
    result.total_time_s = toc(t_total);
    result.first_solution_time_s = read_info_field(info, 'first_solution_time');
    result.iterations = read_info_field(info, 'iterations');
    result.tree_size = read_info_field(info, 'tree_size');
    result.stop_reason = infer_rrt_stop_reason(info);
    result.native_path_cost = read_info_field(info, 'path_cost');
    result.terminated_by_time_limit = read_info_field(info, 'terminated_by_time_limit');
    result.native_summary = extract_native_summary(info, {'first_solution_iter', 'path_risk'});
    if ~result.success
        result.error_message = derive_failure_message(info, result);
    end
catch ME
    result.total_time_s = toc(t_total);
    result.error_message = ME.message;
end
end

function result = init_planner_result(label)
result = struct();
result.label = label;
result.success = false;
result.error_message = 'Planner did not run.';
result.path_world = [];
result.path_ned = [];
result.num_waypoints = 0;
result.preprocess_time_s = NaN;
result.heuristic_time_s = NaN;
result.search_time_s = NaN;
result.total_time_s = NaN;
result.first_solution_time_s = NaN;
result.iterations = NaN;
result.tree_size = NaN;
result.stop_reason = '';
result.path_length_m = NaN;
result.risk_max = NaN;
result.risk_avg = NaN;
result.risk_integrated = NaN;
result.common_stealth_cost = NaN;
result.strict_zero_risk = NaN;
result.dense_hidden_ok = false;
result.dense_max_risk = NaN;
result.dense_visible_fraction = NaN;
result.dense_violation_count = NaN;
result.dense_sample_count = NaN;
result.native_path_cost = NaN;
result.warmstart_cost = NaN;
result.volume_ratio = NaN;
result.fsm_used = NaN;
result.fmm_used = NaN;
result.corridor_time_s = NaN;
result.corridor_goal_connected = NaN;
result.corridor_skeleton_points = NaN;
result.samples_generated = NaN;
result.samples_accepted = NaN;
result.connection_radius_last = NaN;
result.first_solution_batch = NaN;
result.terminated_by_time_limit = NaN;
result.native_summary = struct();
end

function result = finalize_world_path_result(result, path_world, threat, alpha_cost, visibility_threshold, validation_spacing)
result.path_world = path_world;
if isempty(path_world)
    result.path_ned = [];
    result.num_waypoints = 0;
    return;
end

path_ned = [path_world(1, :); path_world(2, :); -path_world(3, :)];
result = finalize_ned_path_result(result, path_ned, threat, alpha_cost, visibility_threshold, validation_spacing);
end

function result = finalize_ned_path_result(result, path_ned, threat, alpha_cost, visibility_threshold, validation_spacing)
result.path_ned = path_ned;
if isempty(path_ned)
    result.path_world = [];
    result.num_waypoints = 0;
    return;
end

result.path_world = [path_ned(1, :); path_ned(2, :); -path_ned(3, :)];
result.num_waypoints = size(path_ned, 2);
result.path_length_m = compute_path_length(path_ned);
[result.risk_max, result.risk_avg, result.risk_integrated] = threat.evaluate_path(path_ned, 1);
result.common_stealth_cost = compute_common_stealth_cost(path_ned, threat, alpha_cost);
[result.dense_hidden_ok, result.dense_max_risk, result.dense_visible_fraction, ...
    result.dense_violation_count, result.dense_sample_count] = ...
    validate_path_zero_risk(path_ned, threat, visibility_threshold, validation_spacing);
end

function comparison = build_comparison_struct(planner_benchmark)
keys = {'ait_star', 'bit_star', 'rrt_star'};
labels = {'AIT*', 'BIT*', 'Informed RRT*'};
colors = [planner_benchmark.cfg.colors.ait_star; ...
          planner_benchmark.cfg.colors.bit_star; ...
          planner_benchmark.cfg.colors.rrt_star];

n = numel(keys);
success = false(1, n);
total_time_s = nan(1, n);
preprocess_time_s = nan(1, n);
search_time_s = nan(1, n);
path_length_m = nan(1, n);
risk_integrated = nan(1, n);
common_stealth_cost = nan(1, n);
dense_hidden_ok = false(1, n);
dense_max_risk = nan(1, n);

for i = 1:n
    item = planner_benchmark.planners.(keys{i});
    success(i) = logical(item.success);
    total_time_s(i) = item.total_time_s;
    preprocess_time_s(i) = item.preprocess_time_s;
    search_time_s(i) = item.search_time_s;
    path_length_m(i) = item.path_length_m;
    risk_integrated(i) = item.risk_integrated;
    common_stealth_cost(i) = item.common_stealth_cost;
    dense_hidden_ok(i) = logical(item.dense_hidden_ok);
    dense_max_risk(i) = item.dense_max_risk;
end

comparison = struct();
comparison.keys = keys;
comparison.labels = labels;
comparison.colors = colors;
comparison.success = success;
comparison.total_time_s = total_time_s;
comparison.preprocess_time_s = preprocess_time_s;
comparison.search_time_s = search_time_s;
comparison.path_length_m = path_length_m;
comparison.risk_integrated = risk_integrated;
comparison.common_stealth_cost = common_stealth_cost;
comparison.dense_hidden_ok = dense_hidden_ok;
comparison.dense_max_risk = dense_max_risk;
comparison.time_stack_s = [preprocess_time_s(:), search_time_s(:)];
end

function plot_benchmark_metrics(planner_benchmark)
cmp = planner_benchmark.comparison;

fig = figure('Name', 'Planner Benchmark Metrics', ...
    'Position', [90, 70, 1400, 820], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Unified Radar Planner Benchmark');

ax1 = nexttile(tl, 1);
plot_scalar_metric(ax1, cmp.total_time_s, cmp.labels, cmp.colors, 'Total Planner Time [s]', 'Total Time');

ax2 = nexttile(tl, 2);
plot_stacked_time_metric(ax2, cmp.time_stack_s, cmp.labels, cmp.colors, cmp.success);

ax3 = nexttile(tl, 3);
plot_scalar_metric(ax3, cmp.path_length_m, cmp.labels, cmp.colors, 'Path Length [m]', 'Path Length');

ax4 = nexttile(tl, 4);
plot_scalar_metric(ax4, cmp.risk_integrated, cmp.labels, cmp.colors, 'Integrated Risk', 'Integrated Risk');

ax5 = nexttile(tl, 5);
plot_scalar_metric(ax5, cmp.common_stealth_cost, cmp.labels, cmp.colors, 'Common Stealth Cost', 'Shared Cost');

ax6 = nexttile(tl, 6);
axis(ax6, 'off');
text(ax6, 0.02, 0.92, sprintf('DEM source: %s', planner_benchmark.environment.dem_source), 'FontSize', 11);
text(ax6, 0.02, 0.80, sprintf('Threat resolution: [%.1f, %.1f] m', ...
    planner_benchmark.environment.threat_resolution(1), planner_benchmark.environment.threat_resolution(2)), 'FontSize', 11);
text(ax6, 0.02, 0.68, sprintf('Setup time: %.2f s', planner_benchmark.shared_times.total_setup_s), 'FontSize', 11);
text(ax6, 0.02, 0.56, sprintf('Start risk: %.3f', planner_benchmark.environment.start_risk), 'FontSize', 11);
text(ax6, 0.02, 0.44, sprintf('Goal risk: %.3f', planner_benchmark.environment.goal_risk), 'FontSize', 11);
for i = 1:numel(cmp.labels)
    planner = planner_benchmark.planners.(cmp.keys{i});
    y = 0.30 - 0.12 * (i - 1);
    if planner.success
        status_text = sprintf('%s: success', cmp.labels{i});
    else
        status_text = sprintf('%s: failed', cmp.labels{i});
    end
    text(ax6, 0.02, y, status_text, 'FontSize', 11, 'Color', cmp.colors(i, :));
end
end

function fig = plot_overlay_figure(planner_benchmark, tm, threat)
cmp = planner_benchmark.comparison;
fig = figure('Name', 'Planner Trajectory Overlay', ...
    'Position', [100, 90, 1360, 920], 'Color', 'w');
ax = axes('Parent', fig);
plot_terrain_surface(ax, tm, []);
hold(ax, 'on');
cb = colorbar(ax);
cb.Label.String = 'Elevation [m]';

for r = 1:numel(threat.radars)
    rr = threat.radars{r};
    draw_radar_sphere(ax, rr, planner_benchmark.cfg.colors.radar);
    plot3(ax, rr.position(1), rr.position(2), rr.position(3), 'kd', ...
        'MarkerFaceColor', planner_benchmark.cfg.colors.radar, 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
end

for i = 1:numel(cmp.keys)
    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success && ~isempty(planner.path_world)
        plot3(ax, planner.path_world(1, :), planner.path_world(2, :), planner.path_world(3, :), ...
            '-', 'Color', cmp.colors(i, :), 'LineWidth', 2.4, 'DisplayName', cmp.labels{i});
    end
end

plot3(ax, planner_benchmark.environment.start_world(1), planner_benchmark.environment.start_world(2), planner_benchmark.environment.start_world(3), ...
    'go', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
plot3(ax, planner_benchmark.environment.goal_world(1), planner_benchmark.environment.goal_world(2), planner_benchmark.environment.goal_world(3), ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 10, 'DisplayName', 'Goal');

style_3d_axes(ax, tm, planner_benchmark);
title(ax, sprintf('%s Radar-Masking Scenario', ...
    getfield_with_default(planner_benchmark.cfg, 'scenario_name', 'Representative scenario')));
legend(ax, 'Location', 'northeastoutside');
end

function plot_planner_panels(planner_benchmark, tm, threat)
cmp = planner_benchmark.comparison;
fig = figure('Name', 'Per-Planner 3D Views', ...
    'Position', [110, 100, 1500, 560], 'Color', 'w');
tl = tiledlayout(fig, 1, numel(cmp.keys), 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(cmp.keys)
    ax = nexttile(tl, i);
    plot_terrain_surface(ax, tm, planner_benchmark.cfg.colors.terrain);
    hold(ax, 'on');

    for r = 1:numel(threat.radars)
        rr = threat.radars{r};
        draw_radar_sphere(ax, rr, planner_benchmark.cfg.colors.radar);
        plot3(ax, rr.position(1), rr.position(2), rr.position(3), 'kd', ...
            'MarkerFaceColor', planner_benchmark.cfg.colors.radar, 'MarkerSize', 7);
    end

    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success && ~isempty(planner.path_world)
        plot3(ax, planner.path_world(1, :), planner.path_world(2, :), planner.path_world(3, :), ...
            '-', 'Color', cmp.colors(i, :), 'LineWidth', 2.5);
    else
        text(ax, mean(tm.bounds(1:2)), mean(tm.bounds(3:4)), max(tm.Z(:)) + 60, ...
            'FAILED', 'Color', [0.7, 0.1, 0.1], 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    plot3(ax, planner_benchmark.environment.start_world(1), planner_benchmark.environment.start_world(2), planner_benchmark.environment.start_world(3), ...
        'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot3(ax, planner_benchmark.environment.goal_world(1), planner_benchmark.environment.goal_world(2), planner_benchmark.environment.goal_world(3), ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

    style_3d_axes(ax, tm, planner_benchmark);
    title(ax, cmp.labels{i});
end
end

function output_files = export_report_assets(planner_benchmark, tm)
output_files = struct();
if ~getfield_with_default(planner_benchmark.cfg, 'report_export', false)
    return;
end

out_dir = planner_benchmark.cfg.report_output_dir;
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

overlay_fig = plot_overlay_figure( ...
    planner_benchmark.representative_scenario.representative_trial, ...
    tm, ...
    planner_benchmark.representative_scenario.representative_trial.runtime.threat);

overlay_path = fullfile(out_dir, 'planner_benchmark_paths.png');
combined_radar_map_path = fullfile(out_dir, 'planner_benchmark_radar_selection_map.png');
summary_path = fullfile(out_dir, 'planner_benchmark_summary.txt');
scenario_csv_path = fullfile(out_dir, 'planner_benchmark_scenarios.csv');
radar_csv_path = fullfile(out_dir, 'planner_benchmark_radars.csv');
radar_map_paths = export_radar_position_maps(planner_benchmark, tm, out_dir);
combined_radar_map_fig = export_combined_radar_position_map(planner_benchmark, tm);
mat_path = fullfile(out_dir, 'planner_benchmark_results.mat');

exportgraphics(overlay_fig, overlay_path, 'Resolution', 240);
exportgraphics(combined_radar_map_fig, combined_radar_map_path, 'Resolution', 240);
write_benchmark_summary_text(planner_benchmark, summary_path);
write_scenario_csv(planner_benchmark, scenario_csv_path);
write_radar_csv(planner_benchmark, radar_csv_path);
benchmark_export = struct( ...
    'cfg', planner_benchmark.cfg, ...
    'num_scenarios', planner_benchmark.num_scenarios, ...
    'scenarios', planner_benchmark.scenarios, ...
    'overlay_figure', overlay_path, ...
    'combined_radar_map', combined_radar_map_path, ...
    'scenario_csv', scenario_csv_path, ...
    'radar_csv', radar_csv_path, ...
    'radar_maps', {radar_map_paths}, ...
    'summary_text', summary_path);
save(mat_path, 'benchmark_export');

output_files.overlay_figure = overlay_path;
output_files.combined_radar_map = combined_radar_map_path;
output_files.scenario_csv = scenario_csv_path;
output_files.radar_csv = radar_csv_path;
output_files.radar_maps = radar_map_paths;
output_files.summary_text = summary_path;
output_files.results_mat = mat_path;
output_files.report_overlay_rel = 'figures/planner_benchmark_paths.png';
output_files.report_combined_radar_map_rel = 'figures/planner_benchmark_radar_selection_map.png';

close(overlay_fig);
close(combined_radar_map_fig);
end

function write_benchmark_summary_text(planner_benchmark, summary_path)
fid = fopen(summary_path, 'w');
if fid < 0
    warning('Could not write benchmark summary file: %s', summary_path);
    return;
end

cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, 'Planner Benchmark Summary\n');
fprintf(fid, 'Scenarios: %d\n\n', planner_benchmark.num_scenarios);

for s = 1:planner_benchmark.num_scenarios
    scenario = planner_benchmark.scenarios(s);
    cmp = scenario.comparison;
    fprintf(fid, '%s\n', scenario.scenario_name);
    fprintf(fid, '  Seeds: %d\n', scenario.num_trials);
    fprintf(fid, '  Representative seed: %d\n', getfield_with_default(scenario.representative_trial, 'trial_seed', NaN));
    fprintf(fid, '  DEM source: %s\n', scenario.representative_trial.environment.dem_source);
    fprintf(fid, '  Binary masking threshold: %.3f\n', scenario.conditions.visibility_threshold);
    if ~isempty(getfield_with_default(scenario, 'scenario_notes', ''))
        fprintf(fid, '  Notes: %s\n', scenario.scenario_notes);
    end
    for r = 1:numel(scenario.representative_trial.environment.radar_sites)
        radar_site_info = scenario.representative_trial.environment.radar_sites(r);
        fprintf(fid, '  %s position [N,E,H]: [%.1f, %.1f, %.1f] m | range: %.1f m\n', ...
            radar_site_info.name, radar_site_info.position(1), radar_site_info.position(2), ...
            radar_site_info.position(3), radar_site_info.range_m);
    end
    for i = 1:numel(cmp.keys)
        fprintf(fid, '  %s\n', cmp.labels{i});
        fprintf(fid, '    Success rate: %.2f\n', cmp.success_rate(i));
        fprintf(fid, '    Total time [s]: %.6f +/- %.6f\n', cmp.total_time_mean_s(i), cmp.total_time_std_s(i));
        fprintf(fid, '    Search time [s]: %.6f +/- %.6f\n', cmp.search_time_mean_s(i), cmp.search_time_std_s(i));
        fprintf(fid, '    Path length [m]: %.6f +/- %.6f\n', cmp.path_length_mean_m(i), cmp.path_length_std_m(i));
        fprintf(fid, '    Integrated risk: %.6f +/- %.6f\n', cmp.risk_integrated_mean(i), cmp.risk_integrated_std(i));
    end
    fprintf(fid, '\n');
end
end

function write_scenario_csv(planner_benchmark, scenario_csv_path)
fid = fopen(scenario_csv_path, 'w');
if fid < 0
    warning('Could not write benchmark CSV file: %s', scenario_csv_path);
    return;
end

cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, 'scenario,planner,success_rate,time_mean_s,time_std_s,path_length_mean_m,path_length_std_m,risk_mean,risk_std\n');
for s = 1:planner_benchmark.num_scenarios
    scenario = planner_benchmark.scenarios(s);
    cmp = scenario.comparison;
    for i = 1:numel(cmp.keys)
        fprintf(fid, '%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            scenario.scenario_slug, cmp.keys{i}, cmp.success_rate(i), ...
            cmp.total_time_mean_s(i), cmp.total_time_std_s(i), ...
            cmp.path_length_mean_m(i), cmp.path_length_std_m(i), ...
            cmp.risk_integrated_mean(i), cmp.risk_integrated_std(i));
    end
end
end

function write_radar_csv(planner_benchmark, radar_csv_path)
fid = fopen(radar_csv_path, 'w');
if fid < 0
    warning('Could not write benchmark radar CSV file: %s', radar_csv_path);
    return;
end

cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, 'scenario,radar_name,north_m,east_m,ground_height_m,radar_height_m,range_m\n');
for s = 1:planner_benchmark.num_scenarios
    scenario = planner_benchmark.scenarios(s);
    radar_sites = scenario.representative_trial.environment.radar_sites;
    radar_ne = scenario.representative_trial.environment.radar_positions_ne;
    for r = 1:numel(radar_sites)
        radar_pos = radar_sites(r).position;
        ground_h = radar_pos(3) - planner_benchmark.cfg.radar_alt_offset;
        fprintf(fid, '%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            scenario.scenario_slug, radar_sites(r).name, radar_ne(1, r), radar_ne(2, r), ...
            ground_h, radar_pos(3), radar_sites(r).range_m);
    end
end
end

function radar_map_paths = export_radar_position_maps(planner_benchmark, tm, out_dir)
radar_map_paths = cell(1, planner_benchmark.num_scenarios);
scenario_colors = [0.86, 0.23, 0.19; ...
                   0.95, 0.65, 0.14; ...
                   0.18, 0.46, 0.76];

for s = 1:planner_benchmark.num_scenarios
    scenario = planner_benchmark.scenarios(s);
    fig = figure('Name', sprintf('%s Radar Position', scenario.scenario_name), ...
        'Position', [120, 120, 980, 760], 'Color', 'w', 'Visible', 'off');
    ax = axes('Parent', fig);
    imagesc(ax, tm.N_vec, tm.E_vec, tm.Z);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, 'tight');
    hold(ax, 'on');
    [N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
    contour(ax, N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
    colormap(ax, benchmark_terrain_colormap());
    cb = colorbar(ax);
    cb.Label.String = 'Elevation [m]';

    color_idx = min(s, size(scenario_colors, 1));
    scenario_color = scenario_colors(color_idx, :);
    env = scenario.representative_trial.environment;
    start_xy = env.start_world(1:2);
    goal_xy = env.goal_world(1:2);
    radar_sites = scenario.representative_trial.environment.radar_sites;
    for r = 1:numel(radar_sites)
        radar_xy = radar_sites(r).position(1:2);
        plot_range_circle(ax, radar_xy, radar_sites(r).range_m, scenario_color, 0.035, 0.18, 1.1);
        plot(ax, radar_xy(1), radar_xy(2), 'd', ...
            'MarkerSize', 13, 'MarkerFaceColor', planner_benchmark.cfg.colors.radar, ...
            'MarkerEdgeColor', 'k');
        text(ax, radar_xy(1) + 45, radar_xy(2) + 45, ...
            sprintf('Radar %d', r), 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10);
    end
    plot(ax, start_xy(1), start_xy(2), 'o', 'MarkerSize', 11, ...
        'MarkerFaceColor', [0.10, 0.85, 0.20], 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
    plot(ax, goal_xy(1), goal_xy(2), 'o', 'MarkerSize', 11, ...
        'MarkerFaceColor', [0.95, 0.10, 0.10], 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
    text(ax, start_xy(1) + 45, start_xy(2) - 55, 'Start', 'Color', 'k', ...
        'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 1.2);
    text(ax, goal_xy(1) + 45, goal_xy(2) - 55, 'Goal', 'Color', 'k', ...
        'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 1.2);

    xlabel(ax, 'North [m]');
    ylabel(ax, 'East [m]');
    title(ax, sprintf('%s: Radar, Start, Goal, and Range Footprint', scenario.scenario_name));

    out_path = fullfile(out_dir, sprintf('%s_radar_map.png', scenario.scenario_slug));
    exportgraphics(fig, out_path, 'Resolution', 240);
    radar_map_paths{s} = out_path;
    close(fig);
end
end

function fig = export_combined_radar_position_map(planner_benchmark, tm)
scenario_colors = [0.86, 0.23, 0.19; ...
                   0.95, 0.65, 0.14; ...
                   0.18, 0.46, 0.76];

fig = figure('Name', 'Benchmark Radar Selection Map', ...
    'Position', [120, 120, 1020, 780], 'Color', 'w', 'Visible', 'off');
ax = axes('Parent', fig);
imagesc(ax, tm.N_vec, tm.E_vec, tm.Z);
set(ax, 'YDir', 'normal');
axis(ax, 'equal');
axis(ax, 'tight');
hold(ax, 'on');
[N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
contour(ax, N_grid, E_grid, tm.Z, 18, 'k-', 'LineWidth', 0.5);
colormap(ax, benchmark_terrain_colormap());
cb = colorbar(ax);
cb.Label.String = 'Elevation [m]';

legend_handles = gobjects(1, planner_benchmark.num_scenarios + 2);
legend_labels = cell(1, planner_benchmark.num_scenarios + 2);
for s = 1:planner_benchmark.num_scenarios
    scenario = planner_benchmark.scenarios(s);
    env = scenario.representative_trial.environment;
    radar_site = scenario.representative_trial.environment.radar_sites(1);
    radar_pos = radar_site.position;
    start_xy = env.start_world(1:2);
    goal_xy = env.goal_world(1:2);
    color_idx = min(s, size(scenario_colors, 1));
    scenario_color = scenario_colors(color_idx, :);
    plot_range_circle(ax, radar_pos(1:2), radar_site.range_m, scenario_color, 0.1, 0.25, 1.0);
    legend_handles(s) = plot(ax, radar_pos(1), radar_pos(2), 'd', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', scenario_color, ...
        'MarkerEdgeColor', 'k', ...
        'LineStyle', 'none');
    legend_labels{s} = sprintf('%s radar', scenario.scenario_name);
    plot(ax, start_xy(1), start_xy(2), '^', ...
        'MarkerSize', 9, 'MarkerFaceColor', scenario_color, ...
        'MarkerEdgeColor', 'k', 'LineStyle', 'none');
    plot(ax, goal_xy(1), goal_xy(2), 'o', ...
        'MarkerSize', 9, 'MarkerFaceColor', scenario_color, ...
        'MarkerEdgeColor', 'k', 'LineStyle', 'none');
    text(ax, radar_pos(1) + 55, radar_pos(2) + 45, scenario.scenario_name, ...
        'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', 'w', 'Margin', 1.5);
    text(ax, start_xy(1) + 45, start_xy(2) - 55, sprintf('S%d start', s), ...
        'Color', 'k', 'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 1.0);
    text(ax, goal_xy(1) + 45, goal_xy(2) - 55, sprintf('S%d goal', s), ...
        'Color', 'k', 'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 1.0);
end

legend_handles(end-1) = plot(ax, nan, nan, '^', 'MarkerSize', 9, ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
legend_labels{end-1} = 'Start';
legend_handles(end) = plot(ax, nan, nan, 'o', 'MarkerSize', 9, ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
legend_labels{end} = 'Goal';

xlabel(ax, 'North [m]');
ylabel(ax, 'East [m]');
title(ax, 'Single DEM Map with Radar Positions, Start/Goal Points, and Radar Ranges');
legend(ax, legend_handles, legend_labels, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');
end

function show_native_report_figures(planner_benchmark, tm)
scenario = planner_benchmark.representative_scenario;

fig_metrics = plot_scenario_comparison_figure(planner_benchmark);
set(fig_metrics, 'Visible', 'on');

fig_overlay = plot_overlay_figure( ...
    scenario.representative_trial, ...
    tm, ...
    scenario.representative_trial.runtime.threat);
set(fig_overlay, 'Visible', 'on');

fig_radar = export_combined_radar_position_map(planner_benchmark, tm);
set(fig_radar, 'Visible', 'on');

drawnow;
end

function fig = plot_scenario_comparison_figure(planner_benchmark)
scenario_names = planner_benchmark.scenario_names;
n_scenarios = numel(scenario_names);
planner_labels = planner_benchmark.scenarios(1).comparison.labels;
planner_colors = planner_benchmark.scenarios(1).comparison.colors;
n_planners = numel(planner_labels);

time_data = nan(n_scenarios, n_planners);
length_data = nan(n_scenarios, n_planners);
risk_data = nan(n_scenarios, n_planners);

for s = 1:n_scenarios
    cmp = planner_benchmark.scenarios(s).comparison;
    time_data(s, :) = cmp.total_time_mean_s(:)';
    length_data(s, :) = cmp.path_length_mean_m(:)';
    risk_data(s, :) = cmp.risk_integrated_mean(:)';
end

fig = figure('Name', 'Planner Benchmark Scenario Comparison', ...
    'Position', [110, 80, 1380, 860], 'Color', 'w', 'Visible', 'off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Scenario-Level Planner Benchmark');

ax1 = nexttile(tl, 1);
plot_grouped_metric(ax1, time_data, scenario_names, planner_labels, planner_colors, 'Total Time [s]', 'Computation Time');

ax2 = nexttile(tl, 2);
plot_grouped_metric(ax2, length_data, scenario_names, planner_labels, planner_colors, 'Path Length [m]', 'Path Length');

ax3 = nexttile(tl, 3);
plot_grouped_metric(ax3, risk_data, scenario_names, planner_labels, planner_colors, 'Integrated Risk', 'Integrated Risk');

ax4 = nexttile(tl, 4);
axis(ax4, 'off');
text(ax4, 0.02, 0.92, sprintf('Representative overlay scenario: %s', scenario_names{planner_benchmark.representative_scenario_index}), 'FontSize', 11);
text(ax4, 0.02, 0.78, sprintf('DEM source: %s', planner_benchmark.representative_scenario.representative_trial.environment.dem_source), 'FontSize', 11);
text(ax4, 0.02, 0.64, sprintf('Binary visibility threshold: %.2f', planner_benchmark.representative_scenario.conditions.visibility_threshold), 'FontSize', 11);
text(ax4, 0.02, 0.50, sprintf('Runs per scenario: %d', planner_benchmark.representative_scenario.num_trials), 'FontSize', 11);
text(ax4, 0.02, 0.36, 'All current scenarios return zero integrated radar risk.', 'FontSize', 11);
legend(ax1, planner_labels, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
end

function plot_grouped_metric(ax, values, group_labels, planner_labels, planner_colors, y_label, ttl)
axes(ax);
h = bar(ax, values, 'grouped');
for i = 1:numel(h)
    h(i).FaceColor = planner_colors(i, :);
end
set(ax, 'XTick', 1:numel(group_labels), 'XTickLabel', group_labels);
ylabel(ax, y_label);
title(ax, ttl);
grid(ax, 'on');
xtickangle(ax, 0);
end

function plot_range_circle(ax, center_xy, radius_m, color_rgb, face_alpha, edge_alpha, line_width)
theta = linspace(0, 2 * pi, 240);
x = center_xy(1) + radius_m * cos(theta);
y = center_xy(2) + radius_m * sin(theta);
patch(ax, x, y, color_rgb, ...
    'FaceAlpha', face_alpha, ...
    'EdgeColor', color_rgb, ...
    'EdgeAlpha', edge_alpha, ...
    'LineWidth', line_width);
end

function plot_mean_std_metric(ax, means, stds, labels, colors, y_label, ttl, success_rate)
axes(ax);
h = bar(ax, means, 0.65, 'FaceColor', 'flat');
for i = 1:numel(means)
    h.CData(i, :) = colors(i, :);
end
hold(ax, 'on');
errorbar(ax, 1:numel(means), means, stds, '.k', 'LineWidth', 1.2, 'CapSize', 12);
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels);
ylabel(ax, y_label);
title(ax, ttl);
grid(ax, 'on');

yl = ylim(ax);
offset = 0.04 * max(yl(2) - yl(1), 1);
for i = 1:numel(means)
    if isfinite(means(i))
        text(ax, i, means(i) + stds(i) + offset, ...
            sprintf('%.2f +/- %.2f\nsucc=%.0f%%', means(i), stds(i), 100 * success_rate(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    else
        text(ax, i, yl(1) + 0.85 * (yl(2) - yl(1)), 'FAILED', ...
            'HorizontalAlignment', 'center', 'Color', [0.7, 0.1, 0.1], ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end
end

function plot_scalar_metric(ax, values, labels, colors, y_label, ttl)
axes(ax);
bar_data = values;
h = bar(ax, bar_data, 0.65, 'FaceColor', 'flat');
for i = 1:numel(values)
    h.CData(i, :) = colors(i, :);
end
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels);
ylabel(ax, y_label);
title(ax, ttl);
grid(ax, 'on');

yl = ylim(ax);
for i = 1:numel(values)
    if isfinite(values(i))
        text(ax, i, values(i) + 0.03 * max(yl(2) - yl(1), 1), sprintf('%.2f', values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    else
        text(ax, i, yl(1) + 0.85 * (yl(2) - yl(1)), 'FAILED', ...
            'HorizontalAlignment', 'center', 'Color', [0.7, 0.1, 0.1], 'FontSize', 10, 'FontWeight', 'bold');
    end
end
end

function plot_stacked_time_metric(ax, time_stack_s, labels, colors, success)
axes(ax);
h = bar(ax, time_stack_s, 'stacked', 'BarWidth', 0.65);
h(1).FaceColor = [0.78, 0.83, 0.89];
h(2).FaceColor = [0.27, 0.48, 0.69];
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels);
ylabel(ax, 'Time [s]');
title(ax, 'Preprocess vs Search');
legend(ax, {'Preprocess', 'Search'}, 'Location', 'best');
grid(ax, 'on');

totals = sum(time_stack_s, 2, 'omitnan');
yl = ylim(ax);
for i = 1:numel(labels)
    if success(i) && isfinite(totals(i))
        text(ax, i, totals(i) + 0.03 * max(yl(2) - yl(1), 1), sprintf('%.2f', totals(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    else
        text(ax, i, yl(1) + 0.85 * (yl(2) - yl(1)), 'FAILED', ...
            'HorizontalAlignment', 'center', 'Color', [0.7, 0.1, 0.1], 'FontSize', 10, 'FontWeight', 'bold');
    end
end
end

function plot_terrain_surface(ax, tm, face_color)
[N_grid, E_grid] = meshgrid(tm.N_vec, tm.E_vec);
if isempty(face_color)
    surf(ax, N_grid, E_grid, tm.Z, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp', ...
        'CData', tm.Z, ...
        'FaceAlpha', 1.0);
    contour3(ax, N_grid, E_grid, tm.Z, 15, 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    colormap(ax, benchmark_terrain_colormap());
    light(ax, 'Position', [-1 -1 1], 'Style', 'infinite');
    lighting(ax, 'gouraud');
    material(ax, 'dull');
else
    surf(ax, N_grid, E_grid, tm.Z, ...
        'EdgeColor', 'none', ...
        'FaceColor', face_color, ...
        'FaceAlpha', 0.92);
    colormap(ax, turbo);
end
end

function style_3d_axes(ax, tm, planner_benchmark)
all_alt = [planner_benchmark.environment.start_world(3), planner_benchmark.environment.goal_world(3)];
cmp = planner_benchmark.comparison;
for i = 1:numel(cmp.keys)
    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success && ~isempty(planner.path_world)
        all_alt = [all_alt, planner.path_world(3, :)]; %#ok<AGROW>
    end
end

terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));
alt_max = max([terrain_max + 100, all_alt + 50]);

xlim(ax, [tm.bounds(1), tm.bounds(2)]);
ylim(ax, [tm.bounds(3), tm.bounds(4)]);
zlim(ax, [terrain_min - 20, alt_max]);
xlabel(ax, 'North [m]');
ylabel(ax, 'East [m]');
zlabel(ax, 'Height [m]');
view(ax, 38, 30);
grid(ax, 'on');
axis(ax, 'tight');
end

function draw_radar_sphere(ax, radar_obj, color_rgb)
[sx, sy, sz] = sphere(28);
surf(ax, ...
    sx * radar_obj.R_max + radar_obj.position(1), ...
    sy * radar_obj.R_max + radar_obj.position(2), ...
    sz * radar_obj.R_max + radar_obj.position(3), ...
    'FaceColor', [1.0, 0.35, 0.15], ...
    'FaceAlpha', 0.05, ...
    'EdgeColor', [1.0, 0.45, 0.18], ...
    'EdgeAlpha', 0.10, ...
    'HandleVisibility', 'off');
end

function cmap = benchmark_terrain_colormap()
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

function print_benchmark_summary(planner_benchmark)
cmp = planner_benchmark.comparison;
fprintf('\n=== Benchmark Summary ===\n');
print_shared_conditions_summary(planner_benchmark.conditions);
fprintf('Shared setup: %.2f s (DEM %.2f s, threat %.2f s)\n', ...
    planner_benchmark.shared_times.total_setup_s, ...
    planner_benchmark.shared_times.dem_load_s, ...
    planner_benchmark.shared_times.threat_compute_s);

for i = 1:numel(cmp.keys)
    planner = planner_benchmark.planners.(cmp.keys{i});
    if planner.success
        fprintf(['%-10s | total=%.2f s | preprocess=%.2f s | search=%.2f s | ', ...
                 'length=%.1f m | hidden=%d | max_risk=%.3f | stealth=%.3f\n'], ...
            cmp.labels{i}, planner.total_time_s, planner.preprocess_time_s, planner.search_time_s, ...
            planner.path_length_m, planner.dense_hidden_ok, planner.dense_max_risk, planner.common_stealth_cost);
    else
        fprintf('%-10s | FAILED | total=%.2f s | preprocess=%.2f s | search=%.2f s | %s\n', ...
            cmp.labels{i}, planner.total_time_s, planner.preprocess_time_s, planner.search_time_s, planner.error_message);
    end
end

fprintf('planner_benchmark exported to workspace.\n');
end

function print_shared_conditions_summary(conditions)
if isempty(conditions)
    return;
end

if conditions.preferred_agl_enabled
    objective_text = sprintf('length + AGL penalty (gamma=%.3f, preferred AGL=%.1f m)', ...
        conditions.preferred_agl_cost_weight, conditions.preferred_agl);
else
    objective_text = 'length only';
end

fprintf(['Shared conditions: objective=%s | strict_zero_risk=%d | vis_thresh=%.3f | ', ...
         'clearance=%.1f m | budget=%.1f s | validation=%.1f m | tuning=%s | ', ...
         'AIT parallel: reverse=%d edge=%d\n'], ...
    objective_text, conditions.strict_zero_risk, conditions.visibility_threshold, ...
    conditions.min_clearance_m, conditions.time_budget_s, conditions.validation_spacing_m, ...
    conditions.planner_tuning_policy, conditions.ait_use_parallel_reverse, conditions.ait_use_parallel_edge_checks);
end

function [NE, h, risk] = choose_endpoint(threat, tm, target_NE, agl, visibility_threshold)
span_n = tm.bounds(2) - tm.bounds(1);
span_e = tm.bounds(4) - tm.bounds(3);
max_offset_n = max(700, 0.38 * span_n);
max_offset_e = max(700, 0.38 * span_e);
step_n = max(120, round(max_offset_n / 10));
step_e = max(120, round(max_offset_e / 10));
offset_n = -max_offset_n:step_n:max_offset_n;
offset_e = -max_offset_e:step_e:max_offset_e;
offsets = zeros(numel(offset_n) * numel(offset_e), 2);
cursor = 1;
for dn = offset_n
    for de = offset_e
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
    visible_penalty = 1e6 * double(cand_risk >= visibility_threshold);
    score = visible_penalty + 1e4 * cand_risk + norm(offsets(i, :));
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

function cost_val = compute_common_stealth_cost(path_ned, threat, alpha)
cost_val = NaN;
if isempty(path_ned) || size(path_ned, 2) < 2
    return;
end

cost_val = 0.0;
for k = 2:size(path_ned, 2)
    p0 = path_ned(:, k - 1);
    p1 = path_ned(:, k);
    seg_len = norm(p1 - p0);
    if seg_len < 1e-9
        continue;
    end
    n_samples = max(2, ceil(seg_len / 30));
    lambda = linspace(0, 1, n_samples);
    seg_cost = 0.0;
    for s = 1:n_samples
        ps = (1 - lambda(s)) * p0 + lambda(s) * p1;
        risk_s = threat.get_risk(ps(1), ps(2), -ps(3));
        seg_cost = seg_cost + exp(alpha * risk_s);
    end
    cost_val = cost_val + seg_len * seg_cost / n_samples;
end
end

function [hidden_ok, max_risk, visible_fraction, violation_count, sample_count] = ...
        validate_path_zero_risk(path_ned, threat, visibility_threshold, sample_spacing)
hidden_ok = false;
max_risk = NaN;
visible_fraction = NaN;
violation_count = NaN;
sample_count = NaN;

if isempty(path_ned) || size(path_ned, 2) < 2
    return;
end

path_dense = resample_path_by_spacing(path_ned, sample_spacing);
if isempty(path_dense)
    return;
end

risks = threat.get_risk(path_dense(1, :), path_dense(2, :), -path_dense(3, :));
risks = reshape(risks, 1, []);
sample_count = numel(risks);
if sample_count == 0
    return;
end

max_risk = max(risks);
violation_count = sum(risks >= visibility_threshold);
visible_fraction = violation_count / sample_count;
hidden_ok = violation_count == 0;
end

function path_dense = resample_path_by_spacing(path_xyz, sample_spacing)
if isempty(path_xyz) || size(path_xyz, 2) < 2
    path_dense = zeros(3, 0);
    return;
end

sample_spacing = max(sample_spacing, 1);
path_dense = path_xyz(:, 1);
for k = 2:size(path_xyz, 2)
    p0 = path_xyz(:, k - 1);
    p1 = path_xyz(:, k);
    seg_len = norm(p1 - p0);
    if seg_len < 1e-9
        continue;
    end
    n_seg = max(2, ceil(seg_len / sample_spacing) + 1);
    t = linspace(0, 1, n_seg);
    seg_pts = p0 * (1 - t) + p1 * t;
    path_dense = [path_dense, seg_pts(:, 2:end)]; %#ok<AGROW>
end
end

function len = compute_path_length(path_xyz)
if isempty(path_xyz) || size(path_xyz, 2) < 2
    len = NaN;
    return;
end
len = sum(vecnorm(diff(path_xyz, 1, 2), 2, 1));
end

function n_samples = compute_edge_check_samples(max_edge_length, sample_spacing)
sample_spacing = max(sample_spacing, 1);
max_edge_length = max(max_edge_length, sample_spacing);
n_samples = max(9, ceil(max_edge_length / sample_spacing) + 1);
end

function x_ned = world_to_ned(x_world)
x_ned = [x_world(1); x_world(2); -x_world(3)];
end

function val = read_info_field(info, field_name)
val = NaN;
if isstruct(info) && isfield(info, field_name) && ~isempty(info.(field_name))
    val = info.(field_name);
end
end

function txt = read_info_text(info, field_name)
txt = '';
if isstruct(info) && isfield(info, field_name) && ~isempty(info.(field_name))
    raw = info.(field_name);
    if isstring(raw)
        txt = char(raw);
    elseif ischar(raw)
        txt = raw;
    end
end
end

function stop_reason = infer_rrt_stop_reason(info)
stop_reason = '';
if isfield(info, 'terminated_by_time_limit') && logical(info.terminated_by_time_limit)
    stop_reason = 'time_limit';
elseif isfield(info, 'success') && logical(info.success)
    stop_reason = 'goal_reached';
else
    stop_reason = 'max_iter_or_failure';
end
end

function msg = derive_failure_message(info, result)
msg = read_info_text(info, 'stop_reason');
if isempty(msg) && isfield(result, 'corridor_goal_connected') && islogical(result.corridor_goal_connected) && ~result.corridor_goal_connected
    msg = 'Strict corridor disconnected under the binary masking threshold.';
end
if isempty(msg) && isfield(info, 'terminated_by_time_limit') && logical(info.terminated_by_time_limit)
    msg = 'Planner hit the time limit before finding a path.';
end
if isempty(msg)
    msg = 'Planner finished without returning a feasible path.';
end
end

function summary = extract_native_summary(info, field_names)
summary = struct();
for i = 1:numel(field_names)
    name = field_names{i};
    if isfield(info, name) && ~isempty(info.(name))
        summary.(name) = info.(name);
    end
end
end

function val = getfield_with_default(s, field_name, default_val)
val = default_val;
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    val = s.(field_name);
end
end

function val = finite_or_nan(val)
if ~isfinite(val)
    val = NaN;
end
end

function val = finite_or_zero(val)
if ~isfinite(val)
    val = 0;
end
end
