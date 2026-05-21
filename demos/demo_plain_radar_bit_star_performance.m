%% DEMO_PLAIN_RADAR_BIT_STAR_PERFORMANCE
% Performance demo for the corridor-free radar-aware BIT* planner.

clearvars -except plain_bit_demo_cfg; clc; close all;
clear plain_radar_bit_star_planner;
rehash;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(this_dir), 'project'));
setup_project_paths();

fprintf('=== Plain Radar BIT* Performance Demo ===\n');
fprintf('Using plain_radar_bit_star_planner: %s\n', which('plain_radar_bit_star_planner'));

cfg = default_config();
if exist('plain_bit_demo_cfg', 'var') && isstruct(plain_bit_demo_cfg)
    cfg = merge_structs(cfg, plain_bit_demo_cfg);
end

rng(cfg.rng_seed, 'twister');

%% 1) Terrain and continuous risk grid
t_setup = tic;
[terrain_data, dem_source] = load_demo_terrain(cfg);
tm = terrain_map(terrain_data);
los = los_checker(tm);

N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
diag_len = hypot(N_span, E_span);
terrain_min = min(tm.Z(:));
terrain_max = max(tm.Z(:));

threat_params = struct();
threat_params.bounds = tm.bounds;
threat_params.resolution = [max(cfg.threat_horiz_res, diag_len / cfg.threat_cells_per_diag), cfg.threat_vert_res];
threat_params.alt_range = [terrain_min + cfg.min_clearance, terrain_max + cfg.max_alt_above_terrain];

threat = threat_map(tm, los, cfg.target_rcs_m2, threat_params);
radars = create_demo_radars(tm, cfg, diag_len);
for r = 1:numel(radars)
    threat.add_radar(radars{r});
end

fprintf('Computing continuous threat map (%s)...\n', cfg.threat_method);
t_threat = tic;
threat.compute_map(cfg.threat_method, struct('use_parallel', cfg.use_parallel_threat, 'show_progress', false));
t_threat_s = toc(t_threat);

env = struct();
env.risk_grid = threat.risk_grid;
env.N_vec = threat.N_vec;
env.E_vec = threat.E_vec;
env.alt_vec = threat.alt_vec;
env.bounds_world = [threat.N_vec(1), threat.N_vec(end), threat.E_vec(1), threat.E_vec(end), ...
                    threat.alt_vec(1), threat.alt_vec(end)];
env.terrain_map = tm;

risk_query = make_risk_query(threat);
setup_time_s = toc(t_setup);
fprintf('Setup complete: DEM=%s | threat=%.2fs | total=%.2fs\n', dem_source, t_threat_s, setup_time_s);

%% 2) Pre-sample comparable start/goal trials
trials = make_trial_endpoints(tm, threat, cfg, diag_len);

%% 3) Run planner across radar weights
num_weights = numel(cfg.weights);
num_trials = cfg.num_trials;
results = repmat(empty_result(), num_weights, num_trials);
representative_paths = cell(num_weights, 1);

for w_idx = 1:num_weights
    weight = cfg.weights(w_idx);
    fprintf('\n=== Weight %.2f (%d/%d) ===\n', weight, w_idx, num_weights);

    for trial_idx = 1:num_trials
        rng(cfg.rng_seed + 1000 * w_idx + trial_idx, 'twister');

        x_start = trials(trial_idx).start_world;
        x_goal = trials(trial_idx).goal_world;
        params = make_planner_params(cfg, diag_len, weight);

        t_plan = tic;
        [path_world, path_cost, info] = plain_radar_bit_star_planner(x_start, x_goal, env, params);
        total_plan_s = toc(t_plan);

        result = empty_result();
        result.weight = weight;
        result.trial = trial_idx;
        result.success = ~isempty(path_world) && info.success;
        result.planning_time = read_info(info, 'planning_time', total_plan_s);
        result.first_solution_time = read_info(info, 'first_solution_time', NaN);
        result.iterations = read_info(info, 'iterations', NaN);
        result.tree_size = read_info(info, 'tree_size', NaN);
        result.samples_generated = read_info(info, 'samples_generated', NaN);
        result.samples_accepted = read_info(info, 'samples_accepted', NaN);
        result.path_cost = path_cost;
        result.stop_reason = read_text_info(info, 'stop_reason', '');

        if result.success
            metrics = validate_and_score_path(path_world, tm, risk_query, env.bounds_world, cfg);
            result.path_length_m = metrics.path_length_m;
            result.mean_risk = metrics.mean_risk;
            result.max_risk = metrics.max_risk;
            result.min_clearance_m = metrics.min_clearance_m;
            result.valid_dense = metrics.valid_dense;
            if isempty(representative_paths{w_idx}) && metrics.valid_dense
                representative_paths{w_idx} = path_world;
            end
        end

        results(w_idx, trial_idx) = result;
        fprintf('  trial %02d/%02d | success=%d | time=%.2fs | len=%.1f | mean risk=%.3f | %s\n', ...
            trial_idx, num_trials, result.success, result.planning_time, ...
            result.path_length_m, result.mean_risk, result.stop_reason);
    end
end

summary = summarize_results(results, cfg.weights);
print_summary(summary);

plain_bit_benchmark = struct();
plain_bit_benchmark.cfg = cfg;
plain_bit_benchmark.dem_source = dem_source;
plain_bit_benchmark.setup_time_s = setup_time_s;
plain_bit_benchmark.threat_time_s = t_threat_s;
plain_bit_benchmark.results = results;
plain_bit_benchmark.summary = summary;
plain_bit_benchmark.trials = trials;
plain_bit_benchmark.representative_paths = representative_paths;
plain_bit_benchmark.runtime = struct('terrain', tm, 'threat', threat, 'radars', {radars});

assignin('base', 'plain_bit_benchmark', plain_bit_benchmark);

%% 4) Figures
if cfg.show_figures
    plot_terrain_paths(tm, radars, representative_paths, cfg);
    plot_risk_paths(threat, radars, representative_paths, trials, cfg);
    plot_summary_bars(summary, cfg);
    plot_length_risk_scatter(results, cfg);
end

fprintf('\nplain_bit_benchmark exported to workspace.\n');

%% ---------------- Local helpers ----------------
function cfg = default_config()
cfg = struct();
cfg.dem_file = 'DEM/agri.tif';
cfg.dem_target_resolution = 60;
cfg.dem_fill_nodata = 'nearest';
cfg.dem_crop_half_size = 1800;
cfg.rng_seed = 17;
cfg.num_trials = 20;
cfg.weights = [0, 2.5, 5, 10];
cfg.min_clearance = 10;
cfg.endpoint_agl = 120;
cfg.endpoint_margin_frac = 0.10;
cfg.min_start_goal_dist_frac = 0.55;
cfg.target_rcs_m2 = 0.1;
cfg.threat_method = 'prob';
cfg.threat_horiz_res = 70;
cfg.threat_vert_res = 45;
cfg.threat_cells_per_diag = 100;
cfg.max_alt_above_terrain = 380;
cfg.use_parallel_threat = false;
cfg.radar_alt_offset = 25;
cfg.radar_range = 2400;
cfg.radar_power = 120e3;
cfg.radar_fractions = [0.22, 0.50, 0.78; ...
                       0.26, 0.62, 0.30];
cfg.max_batches = 100;
cfg.batch_size = 320;
cfg.max_nodes = 18000;
cfg.post_solution_batches = 10;
cfg.max_time = 20;
cfg.edge_sample_spacing = 35;
cfg.max_edge_samples = 40;
cfg.goal_radius_frac = 0.018;
cfg.min_connection_radius_frac = 0.012;
cfg.max_connection_radius_frac = 0.075;
cfg.validation_spacing = 20;
cfg.show_figures = true;
cfg.colors = lines(numel(cfg.weights));
end

function cfg = merge_structs(cfg, override)
fields = fieldnames(override);
for i = 1:numel(fields)
    cfg.(fields{i}) = override.(fields{i});
end
if ~isfield(cfg, 'colors') || size(cfg.colors, 1) < numel(cfg.weights)
    cfg.colors = lines(numel(cfg.weights));
end
end

function [terrain_data, dem_source] = load_demo_terrain(cfg)
terrain_data = [];
dem_source = 'synthetic_fallback';

try
    if exist(cfg.dem_file, 'file')
        dem_params = struct();
        dem_params.target_resolution = cfg.dem_target_resolution;
        dem_params.fill_nodata = cfg.dem_fill_nodata;
        if cfg.dem_crop_half_size > 0
            h = cfg.dem_crop_half_size;
            dem_params.crop_bounds = [-h, h, -h, h];
        end
        terrain_data = dem_loader(cfg.dem_file, dem_params);
        dem_source = cfg.dem_file;
    end
catch ME
    warning('demo_plain_radar_bit_star_performance:DEMLoadFailed', ...
        'DEM load failed (%s). Using synthetic fallback.', ME.message);
    terrain_data = [];
end

if isempty(terrain_data)
    tparams = struct();
    tparams.type = 'ridge';
    tparams.bounds = [-1800, 1800, -1800, 1800];
    tparams.resolution = cfg.dem_target_resolution;
    tparams.base_height = 80;
    tparams.amplitude = 220;
    tparams.wavelength = 700;
    tparams.seed = cfg.rng_seed;
    terrain_data = terrain_generator(tparams);
end
end

function radars = create_demo_radars(tm, cfg, diag_len)
fractions = cfg.radar_fractions;
num_radars = size(fractions, 2);
radars = cell(1, num_radars);
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
for r = 1:num_radars
    N = tm.bounds(1) + fractions(1, r) * N_span;
    E = tm.bounds(3) + fractions(2, r) * E_span;
    h = tm.get_height(N, E);
    params = struct();
    params.P_t = cfg.radar_power;
    params.R_max = min(cfg.radar_range, 0.85 * diag_len);
    params.elevation_limits = [-8, 88];
    radars{r} = radar_site([N; E; h + cfg.radar_alt_offset], sprintf('Plain-Radar-%d', r), 'tracking', params);
end
end

function trials = make_trial_endpoints(tm, threat, cfg, diag_len)
trials = repmat(struct('start_world', [], 'goal_world', [], 'start_risk', NaN, 'goal_risk', NaN), 1, cfg.num_trials);
for i = 1:cfg.num_trials
    found = false;
    for attempt = 1:200
        s = sample_endpoint(tm, cfg);
        g = sample_endpoint(tm, cfg);
        if isempty(s) || isempty(g)
            continue;
        end
        if norm(s(1:2) - g(1:2)) < cfg.min_start_goal_dist_frac * diag_len
            continue;
        end
        trials(i).start_world = s;
        trials(i).goal_world = g;
        trials(i).start_risk = threat.get_risk(s(1), s(2), s(3));
        trials(i).goal_risk = threat.get_risk(g(1), g(2), g(3));
        found = true;
        break;
    end
    if ~found
        error('Could not sample a valid start/goal pair for trial %d.', i);
    end
end
end

function p = sample_endpoint(tm, cfg)
p = [];
N_span = tm.bounds(2) - tm.bounds(1);
E_span = tm.bounds(4) - tm.bounds(3);
margin_n = cfg.endpoint_margin_frac * N_span;
margin_e = cfg.endpoint_margin_frac * E_span;
for attempt = 1:80
    N = tm.bounds(1) + margin_n + rand() * max(N_span - 2 * margin_n, 1);
    E = tm.bounds(3) + margin_e + rand() * max(E_span - 2 * margin_e, 1);
    h = tm.get_height(N, E);
    if isfinite(h)
        p = [N; E; h + cfg.endpoint_agl];
        return;
    end
end
end

function params = make_planner_params(cfg, diag_len, weight)
params = struct();
params.stealth_weight = weight;
params.point_clearance = cfg.min_clearance;
params.edge_sample_spacing = cfg.edge_sample_spacing;
params.max_edge_samples = cfg.max_edge_samples;
params.max_batches = cfg.max_batches;
params.batch_size = cfg.batch_size;
params.max_nodes = cfg.max_nodes;
params.post_solution_batches = cfg.post_solution_batches;
params.max_time = cfg.max_time;
params.goal_radius = max(25, cfg.goal_radius_frac * diag_len);
params.min_connection_radius = max(25, cfg.min_connection_radius_frac * diag_len);
params.max_connection_radius = max(90, cfg.max_connection_radius_frac * diag_len);
params.eta_radius = 1.8;
end

function risk_query = make_risk_query(threat)
F = griddedInterpolant({threat.E_vec, threat.N_vec, threat.alt_vec}, threat.risk_grid, 'linear', 'nearest');
risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);
end

function metrics = validate_and_score_path(path_world, tm, risk_query, bounds, cfg)
dense = resample_path_by_spacing(path_world, cfg.validation_spacing);
in_bounds = dense(1, :) >= bounds(1) & dense(1, :) <= bounds(2) & ...
            dense(2, :) >= bounds(3) & dense(2, :) <= bounds(4) & ...
            dense(3, :) >= bounds(5) & dense(3, :) <= bounds(6);
h = tm.get_height(dense(1, :), dense(2, :));
clearance = dense(3, :) - h;
risk = risk_query(dense);
risk = max(0, min(1, risk));

metrics = struct();
metrics.path_length_m = compute_path_length(path_world);
metrics.mean_risk = mean(risk(isfinite(risk)));
metrics.max_risk = max(risk(isfinite(risk)));
metrics.min_clearance_m = min(clearance(isfinite(clearance)));
metrics.valid_dense = all(in_bounds) && all(isfinite(risk)) && ...
    all(isfinite(clearance)) && metrics.min_clearance_m >= cfg.min_clearance;
end

function dense = resample_path_by_spacing(path_world, spacing)
if isempty(path_world) || size(path_world, 2) < 2
    dense = zeros(3, 0);
    return;
end
dense = path_world(:, 1);
for k = 1:size(path_world, 2)-1
    p0 = path_world(:, k);
    p1 = path_world(:, k+1);
    seg_len = norm(p1 - p0);
    n = max(2, ceil(seg_len / spacing) + 1);
    t = linspace(0, 1, n);
    seg = p0 .* (1 - t) + p1 .* t;
    dense = [dense, seg(:, 2:end)]; %#ok<AGROW>
end
end

function len = compute_path_length(path_world)
if isempty(path_world) || size(path_world, 2) < 2
    len = NaN;
else
    len = sum(vecnorm(diff(path_world, 1, 2), 2, 1));
end
end

function result = empty_result()
result = struct();
result.weight = NaN;
result.trial = NaN;
result.success = false;
result.planning_time = NaN;
result.first_solution_time = NaN;
result.iterations = NaN;
result.tree_size = NaN;
result.samples_generated = NaN;
result.samples_accepted = NaN;
result.path_cost = NaN;
result.path_length_m = NaN;
result.mean_risk = NaN;
result.max_risk = NaN;
result.min_clearance_m = NaN;
result.valid_dense = false;
result.stop_reason = '';
end

function summary = summarize_results(results, weights)
num_weights = numel(weights);
summary = repmat(struct(), 1, num_weights);
for w = 1:num_weights
    r = results(w, :);
    success = [r.success] & [r.valid_dense];
    accepted_ratio = [r.samples_accepted] ./ max([r.samples_generated], 1);
    summary(w).weight = weights(w);
    summary(w).success_rate = mean(success);
    summary(w).mean_planning_time = mean([r(success).planning_time], 'omitnan');
    summary(w).median_planning_time = median([r(success).planning_time], 'omitnan');
    summary(w).mean_first_solution_time = mean([r(success).first_solution_time], 'omitnan');
    summary(w).mean_path_length_m = mean([r(success).path_length_m], 'omitnan');
    summary(w).mean_risk = mean([r(success).mean_risk], 'omitnan');
    summary(w).max_risk = mean([r(success).max_risk], 'omitnan');
    summary(w).mean_accepted_ratio = mean(accepted_ratio, 'omitnan');
    summary(w).mean_tree_size = mean([r.tree_size], 'omitnan');
end
end

function print_summary(summary)
fprintf('\n=== Plain BIT* Summary ===\n');
for i = 1:numel(summary)
    s = summary(i);
    fprintf(['W=%5.2f | success=%5.1f%% | mean time=%.2fs | median time=%.2fs | ', ...
        'first=%.2fs | len=%.1fm | mean risk=%.3f | max risk=%.3f | accepted=%.2f%% | tree=%.0f\n'], ...
        s.weight, 100 * s.success_rate, s.mean_planning_time, s.median_planning_time, ...
        s.mean_first_solution_time, s.mean_path_length_m, s.mean_risk, s.max_risk, ...
        100 * s.mean_accepted_ratio, s.mean_tree_size);
end
end

function plot_terrain_paths(tm, radars, representative_paths, cfg)
fig = figure('Name', 'Plain Radar BIT*: Terrain Paths', 'Position', [80, 80, 1180, 820], 'Color', 'w');
tm.plot(fig);
hold on;
for r = 1:numel(radars)
    rd = radars{r};
    plot3(rd.position(1), rd.position(2), rd.position(3), 'kd', ...
        'MarkerFaceColor', [1.0, 0.55, 0.1], 'MarkerSize', 9);
end
for w = 1:numel(representative_paths)
    path = representative_paths{w};
    if ~isempty(path)
        plot3(path(1, :), path(2, :), path(3, :), '-', ...
            'Color', cfg.colors(w, :), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('W=%.1f', cfg.weights(w)));
    end
end
title('Representative Plain Radar BIT* Paths on Terrain');
legend('Location', 'northeastoutside');
view(38, 30);
grid on;
end

function plot_risk_paths(threat, radars, representative_paths, trials, cfg)
figure('Name', 'Plain Radar BIT*: Top-Down Risk', 'Position', [110, 100, 1050, 760], 'Color', 'w');
max_risk = max(threat.risk_grid, [], 3);
imagesc(threat.N_vec, threat.E_vec, max_risk);
set(gca, 'YDir', 'normal');
axis equal tight;
hold on;
colormap(turbo);
colorbar;
for r = 1:numel(radars)
    rd = radars{r};
    plot(rd.position(1), rd.position(2), 'kd', 'MarkerFaceColor', [1.0, 0.55, 0.1], 'MarkerSize', 8);
end
for w = 1:numel(representative_paths)
    path = representative_paths{w};
    if ~isempty(path)
        plot(path(1, :), path(2, :), '-', 'Color', cfg.colors(w, :), 'LineWidth', 2.0);
    end
end
plot(trials(1).start_world(1), trials(1).start_world(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(trials(1).goal_world(1), trials(1).goal_world(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('North [m]');
ylabel('East [m]');
title('Max Continuous Radar Risk with Representative Paths');
end

function plot_summary_bars(summary, ~)
weights = [summary.weight];
fig = figure('Name', 'Plain Radar BIT*: Summary Metrics', 'Position', [120, 120, 1180, 720], 'Color', 'w');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(weights, 100 * [summary.success_rate], 'FaceColor', [0.2, 0.45, 0.75]);
ylabel('Success [%]');
xlabel('Radar weight W');
grid on;

nexttile;
bar(weights, [summary.mean_planning_time], 'FaceColor', [0.25, 0.65, 0.45]);
ylabel('Mean planning time [s]');
xlabel('Radar weight W');
grid on;

nexttile;
bar(weights, [summary.mean_path_length_m], 'FaceColor', [0.75, 0.45, 0.25]);
ylabel('Mean path length [m]');
xlabel('Radar weight W');
grid on;

nexttile;
bar(weights, [summary.mean_risk], 'FaceColor', [0.75, 0.25, 0.35]);
ylabel('Mean path risk');
xlabel('Radar weight W');
grid on;
end

function plot_length_risk_scatter(results, cfg)
figure('Name', 'Plain Radar BIT*: Length vs Risk', 'Position', [140, 140, 900, 650], 'Color', 'w');
hold on;
for w = 1:numel(cfg.weights)
    r = results(w, :);
    success = [r.success] & [r.valid_dense];
    if any(success)
        scatter([r(success).path_length_m], [r(success).mean_risk], 46, ...
            cfg.colors(w, :), 'filled', 'DisplayName', sprintf('W=%.1f', cfg.weights(w)));
    end
end
xlabel('Path length [m]');
ylabel('Mean radar risk');
title('Path Length vs Radar Exposure');
legend('Location', 'best');
grid on;
end

function val = read_info(info, field_name, default_val)
if isstruct(info) && isfield(info, field_name) && ~isempty(info.(field_name))
    val = info.(field_name);
else
    val = default_val;
end
end

function txt = read_text_info(info, field_name, default_val)
txt = default_val;
if isstruct(info) && isfield(info, field_name) && ~isempty(info.(field_name))
    raw = info.(field_name);
    if isstring(raw)
        txt = char(raw);
    elseif ischar(raw)
        txt = raw;
    end
end
end
