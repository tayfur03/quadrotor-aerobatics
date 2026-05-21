function manifest = generate_radar_uav_dataset(cfg)
%GENERATE_RADAR_UAV_DATASET Generate AIT*-planned radar-aware UAV training data.
%
% Coordinates exported to HDF5:
%   x -> East [m], centered at DEM midpoint
%   y -> North [m], centered at DEM midpoint
%   z -> height above local terrain (AGL) [m]
% Planner state is [x; y; z_ASL]. The AIT* planner receives an AGL-indexed
% risk grid through build_direct_planner_context via env.risk_altitude_mode.

if nargin < 1
    cfg = struct();
end

repo_root = initialize_project();
cfg = default_config(cfg, repo_root);
cfg = configure_casadi(cfg);
cfg = configure_planner_backend(cfg);

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

parallel_status = maybe_start_pool(cfg);
rng(cfg.base_seed, 'twister');

map_specs = build_map_specs(cfg);
manifest = repmat(struct('file', '', 'seed', NaN, 'source', '', ...
    'synthetic', false, 'num_radars', NaN, 'num_trajectories', NaN, ...
    'mean_rcs', NaN, 'elapsed_s', NaN, 'skipped', false, 'skip_reason', ''), 1, numel(map_specs));

fprintf('=== AIT* Radar UAV Dataset Generator ===\n');
fprintf('Maps=%d | K=%d | N=%d | parallel=%d\n', numel(map_specs), cfg.K, cfg.N_samples, cfg.use_parallel);
fprintf('Parallel status: %s\n', parallel_status.message);
fprintf('Output: %s\n', cfg.output_dir);
fprintf('Resume: start_map_idx=%d | end_map_idx=%s | skip_existing=%d | partial_maps=%d\n', ...
    cfg.start_map_idx, format_optional_idx(cfg.end_map_idx), cfg.skip_existing, cfg.write_partial_maps);

for map_idx = 1:numel(map_specs)
    map_number = map_idx - 1;
    if map_number < cfg.start_map_idx || (~isempty(cfg.end_map_idx) && map_number > cfg.end_map_idx)
        manifest(map_idx).seed = cfg.base_seed + 10000 * map_number;
        manifest(map_idx).source = map_specs(map_idx).source;
        manifest(map_idx).synthetic = map_specs(map_idx).synthetic;
        manifest(map_idx).skipped = true;
        manifest(map_idx).skip_reason = 'Outside requested resume range.';
        continue;
    end

    map_seed = cfg.base_seed + 10000 * map_number;
    rng(map_seed, 'twister');
    map_tic = tic;
    out_file = fullfile(cfg.output_dir, sprintf('map_%04d.h5', map_number));

    fprintf('\n=== map_%04d | seed=%d | %s ===\n', map_number, map_seed, map_specs(map_idx).source);
    if cfg.skip_existing && exist(out_file, 'file') == 2
        existing_k = read_existing_trajectory_count(out_file);
        fprintf('map_%04d.h5: existing file found, skipping (trajectories=%d)\n', map_number, existing_k);
        manifest(map_idx).file = out_file;
        manifest(map_idx).seed = map_seed;
        manifest(map_idx).source = map_specs(map_idx).source;
        manifest(map_idx).synthetic = map_specs(map_idx).synthetic;
        manifest(map_idx).num_radars = NaN;
        manifest(map_idx).num_trajectories = existing_k;
        manifest(map_idx).mean_rcs = NaN;
        manifest(map_idx).elapsed_s = 0;
        manifest(map_idx).skipped = false;
        manifest(map_idx).skip_reason = 'Existing output preserved.';
        save_generation_checkpoint(cfg, manifest);
        continue;
    end

    try
        map = load_or_create_map(map_specs(map_idx), cfg, map_seed);
        radars = scatter_radars_enu(map, cfg);
        sandwich = build_sandwich_maps(map, radars, cfg);
        [cost_field, z_layers_agl] = build_radar_cost_field_agl(map, radars, sandwich, cfg);
        mex_cleanup = initialize_mex_planner_for_map(map, cost_field, z_layers_agl, cfg); %#ok<NASGU>
        [waypoints, starts, goals, rcs_costs, planner_times, trajectory_dts, durations, radar_metrics, target_weights, visible_counts, visibility_budget_used, feasibility_stage, num_valid] = generate_ait_trajectories(map, radars, sandwich, cost_field, z_layers_agl, cfg, map_seed);
    catch ME
        elapsed_s = toc(map_tic);
        fprintf('map_%04d.h5: skipped after %.1fs (%s)\n', map_number, elapsed_s, ME.message);
        manifest(map_idx).file = '';
        manifest(map_idx).seed = map_seed;
        manifest(map_idx).source = map_specs(map_idx).source;
        manifest(map_idx).synthetic = map_specs(map_idx).synthetic;
        manifest(map_idx).num_radars = NaN;
        manifest(map_idx).num_trajectories = 0;
        manifest(map_idx).mean_rcs = NaN;
        manifest(map_idx).elapsed_s = elapsed_s;
        manifest(map_idx).skipped = true;
        manifest(map_idx).skip_reason = ME.message;
        save_generation_checkpoint(cfg, manifest);
        continue;
    end

    if num_valid < cfg.K
        elapsed_s = toc(map_tic);
        wrote_partial = false;
        if num_valid > 0 && cfg.write_partial_maps
            bundle = make_bundle(map, radars, sandwich, cost_field, z_layers_agl, ...
                waypoints, starts, goals, rcs_costs, planner_times, trajectory_dts, durations, radar_metrics, ...
                target_weights, visible_counts, visibility_budget_used, feasibility_stage, num_valid, cfg);
            bundle.partial = true;
            validate_bundle(bundle, cfg);
            write_bundle_h5(out_file, bundle, cfg);
            wrote_partial = true;
        end
        fprintf('map_%04d.h5: %s, generated %d/%d trajectories after %d attempts (%.1fs)\n', ...
            map_number, ternary(wrote_partial, 'partial saved', 'skipped'), ...
            num_valid, cfg.K, cfg.max_path_attempts, elapsed_s);
        manifest(map_idx).file = ternary(wrote_partial, out_file, '');
        manifest(map_idx).seed = map_seed;
        manifest(map_idx).source = map_specs(map_idx).source;
        manifest(map_idx).synthetic = map_specs(map_idx).synthetic;
        manifest(map_idx).num_radars = size(radars, 1);
        manifest(map_idx).num_trajectories = num_valid;
        if num_valid > 0
            manifest(map_idx).mean_rcs = mean(rcs_costs(1:num_valid), 'omitnan');
        else
            manifest(map_idx).mean_rcs = NaN;
        end
        manifest(map_idx).elapsed_s = elapsed_s;
        manifest(map_idx).skipped = true;
        manifest(map_idx).skip_reason = sprintf('Only %d/%d valid trajectories after %d attempts.', ...
            num_valid, cfg.K, cfg.max_path_attempts);
        save_generation_checkpoint(cfg, manifest);
        continue;
    end

    bundle = make_bundle(map, radars, sandwich, cost_field, z_layers_agl, ...
        waypoints, starts, goals, rcs_costs, planner_times, trajectory_dts, durations, radar_metrics, ...
        target_weights, visible_counts, visibility_budget_used, feasibility_stage, num_valid, cfg);

    validate_bundle(bundle, cfg);
    write_bundle_h5(out_file, bundle, cfg);

    elapsed_s = toc(map_tic);
    mean_rcs = mean(rcs_costs, 'omitnan');
    fprintf('map_%04d.h5: K=%d, mean_rcs=%.3f, ok (%.1fs)\n', map_idx - 1, cfg.K, mean_rcs, elapsed_s);

    manifest(map_idx).file = out_file;
    manifest(map_idx).seed = map_seed;
    manifest(map_idx).source = map_specs(map_idx).source;
    manifest(map_idx).synthetic = map_specs(map_idx).synthetic;
    manifest(map_idx).num_radars = size(radars, 1);
    manifest(map_idx).num_trajectories = num_valid;
    manifest(map_idx).mean_rcs = mean_rcs;
    manifest(map_idx).elapsed_s = elapsed_s;
    manifest(map_idx).skipped = false;
    manifest(map_idx).skip_reason = '';
    save_generation_checkpoint(cfg, manifest);
end
save_generation_checkpoint(cfg, manifest);
end

function repo_root = initialize_project()
this_file = mfilename('fullpath');
repo_root = fileparts(this_file);
addpath(fullfile(repo_root, 'project'));
if exist('setup_project_paths', 'file') == 2
    setup_project_paths();
else
    addpath(fullfile(repo_root, 'terrain'));
    addpath(fullfile(repo_root, 'radar'));
    addpath(fullfile(repo_root, 'motion_planner'));
    addpath(fullfile(repo_root, 'flight'));
end
end

function bundle = make_bundle(map, radars, sandwich, cost_field, z_layers_agl, waypoints, starts, goals, rcs_costs, planner_times, trajectory_dts, durations, radar_metrics, target_weights, visible_counts, visibility_budget_used, feasibility_stage, num_valid, cfg)
idx = 1:num_valid;
bundle = struct();
bundle.elevation = single(map.elevation);
bundle.min_visible_alt = single(sandwich.min_visible_alt);
bundle.thickness = single(sandwich.thickness);
bundle.scene_input = single(cat(3, bundle.elevation, bundle.min_visible_alt, bundle.thickness));
bundle.resolution_m = single(map.resolution_m);
bundle.origin_xy = single(map.origin_xy);
bundle.radar_positions = single(radars);
bundle.cost_field = single(cost_field);
bundle.z_layers = single(z_layers_agl(:));
bundle.cost_resolution_m = single(map.resolution_m);
bundle.waypoints = single(waypoints(idx, :, :));
bundle.starts = single(starts(idx, :));
bundle.goals = single(goals(idx, :));
bundle.rcs_costs = single(rcs_costs(idx));
bundle.planner_times = single(planner_times(idx));
bundle.trajectory_dts = single(trajectory_dts(idx));
bundle.durations = single(durations(idx));
bundle.scenario_ids = int32(radar_metrics.scenario_ids(idx));
bundle.straight_rcs_costs = single(radar_metrics.straight_rcs_costs(idx));
bundle.straight_visible_ratios = single(radar_metrics.straight_visible_ratios(idx));
bundle.planned_visible_ratios = single(radar_metrics.planned_visible_ratios(idx));
bundle.radar_improvement_ratios = single(radar_metrics.radar_improvement_ratios(idx));
bundle.target_weights = single(target_weights(idx));
bundle.visible_counts = int32(visible_counts(idx));
bundle.visibility_budget_used = int32(visibility_budget_used(idx));
bundle.feasibility_stage = int32(feasibility_stage(idx));
region = build_region_bundle(map, sandwich, cost_field, z_layers_agl, bundle.waypoints, cfg);
bundle.region_heuristic = single(region.heuristic);
bundle.region_safe_mask = single(region.safe_mask);
bundle.region_threat_distance = single(region.threat_distance);
bundle.num_trajectories = int32(num_valid);
bundle.requested_K = int32(cfg.K);
bundle.partial = num_valid < cfg.K;
bundle.dt = single(cfg.dt);
bundle.dt_mode = cfg.dt_mode;
if cfg.export_time_channel
    bundle.waypoint_channels = 'x,y,z_agl,t';
else
    bundle.waypoint_channels = 'x,y,z_agl';
end
bundle.v_max = single(cfg.v_max);
bundle.a_max = single(cfg.a_max);
bundle.start_goal_mode = cfg.start_goal_mode;
bundle.enforce_start_goal_reachability = cfg.enforce_start_goal_reachability;
bundle.solver_backend = cfg.local_solver_resolved;
bundle.casadi_version = cfg.casadi_version;
bundle.planner_backend = cfg.planner_backend;
bundle.generated_trajectory = single(cfg.generate_trajectory);
bundle.used_local_optimizer = single(cfg.generate_trajectory && cfg.use_local_optimizer && cfg.casadi_available);
if ~cfg.generate_trajectory
    bundle.path_source = 'ait_star_raw_resampled';
    bundle.planner_label = sprintf('%s_raw_region', cfg.planner_backend);
elseif cfg.use_local_optimizer && cfg.casadi_available
    bundle.path_source = sprintf('minsnap_casadi_%s', cfg.local_solver_resolved);
    bundle.planner_label = sprintf('%s_minsnap_casadi_%s', cfg.planner_backend, cfg.local_solver_resolved);
else
    bundle.path_source = 'minsnap';
    bundle.planner_label = sprintf('%s_minsnap', cfg.planner_backend);
end
if cfg.export_bspline
    bundle.bspline_control_points = single(fit_bspline_control_points(bundle.waypoints(:, :, 1:3), cfg.bspline_control_points));
end
end

function region = build_region_bundle(map, sandwich, cost_field, z_layers_agl, waypoints, cfg)
H = cfg.dem_size(1);
W = cfg.dem_size(2);
safe_mask = build_region_safe_mask(map, sandwich, cost_field, z_layers_agl, cfg);
threat_distance = build_region_threat_distance_m(safe_mask, double(map.resolution_m));
heuristic = zeros(H, W, 'single');

if ~use_region_label_mode(cfg) || isempty(waypoints)
    region = struct('heuristic', heuristic, 'safe_mask', single(safe_mask), ...
        'threat_distance', single(threat_distance));
    return;
end

if strcmpi(string(cfg.region_corridor_mode), "threat_distance_adaptive")
    heuristic = build_adaptive_corridor_heuristic(map, waypoints, safe_mask, threat_distance, cfg);
    region = struct('heuristic', heuristic, 'safe_mask', single(safe_mask), ...
        'threat_distance', single(threat_distance));
    return;
end

positive_seed = false(H, W);
visible_seed = false(H, W);
K = size(waypoints, 1);
for k = 1:K
    wp = double(squeeze(waypoints(k, :, 1:3)));
    [hidden, visible] = waypoint_visibility_masks(map, sandwich, cost_field, z_layers_agl, wp, cfg);
    [rows, cols] = xy_to_grid_indices(map, wp(:, 1), wp(:, 2));
    valid = rows >= 1 & rows <= H & cols >= 1 & cols <= W & isfinite(wp(:, 1)) & isfinite(wp(:, 2));
    if any(valid & hidden)
        positive_seed(sub2ind([H, W], rows(valid & hidden), cols(valid & hidden))) = true;
    end
    if any(valid & visible)
        visible_seed(sub2ind([H, W], rows(valid & visible), cols(valid & visible))) = true;
    end
end

if any(positive_seed(:))
    dist_positive = distance_to_seed_m(positive_seed, double(map.resolution_m));
    sigma = max(double(cfg.region_corridor_sigma_m), double(map.resolution_m));
    heuristic = single(exp(-(dist_positive.^2) ./ (2 * sigma^2)));

    if any(visible_seed(:))
        dist_visible = distance_to_seed_m(visible_seed, double(map.resolution_m));
        sigma_v = max(double(cfg.region_visible_exclusion_sigma_m), double(map.resolution_m));
        visible_suppression = 1 - exp(-(dist_visible.^2) ./ (2 * sigma_v^2));
        heuristic = heuristic .* single(visible_suppression);
    end
    heuristic = heuristic .* single(safe_mask);
    max_h = max(heuristic(:));
    if max_h > 0
        heuristic = heuristic ./ max_h;
    end
end

heuristic = single(min(max(heuristic, 0), 1));
region = struct('heuristic', heuristic, 'safe_mask', single(safe_mask), ...
    'threat_distance', single(threat_distance));
end

function threat_distance_m = build_region_threat_distance_m(safe_mask, resolution_m)
unsafe_mask = ~logical(safe_mask);
if all(~unsafe_mask(:))
    threat_distance_m = inf(size(safe_mask));
elseif exist('bwdist', 'file') == 2
    threat_distance_m = bwdist(unsafe_mask) * resolution_m;
else
    threat_distance_m = distance_to_seed_m(unsafe_mask, resolution_m);
end
finite_vals = threat_distance_m(isfinite(threat_distance_m));
if isempty(finite_vals)
    threat_distance_m(:) = 0;
else
    fill_val = max(finite_vals(:));
    threat_distance_m(~isfinite(threat_distance_m)) = fill_val;
end
threat_distance_m = single(max(threat_distance_m, 0));
end

function heuristic = build_adaptive_corridor_heuristic(map, waypoints, safe_mask, threat_distance_m, cfg)
H = cfg.dem_size(1);
W = cfg.dem_size(2);
heuristic = zeros(H, W, 'single');
[X_grid, Y_grid] = meshgrid(double(map.xs), double(map.ys));
sigma_min = max(double(cfg.region_sigma_min_m), double(map.resolution_m));
sigma_max = max(sigma_min, double(cfg.region_sigma_max_m));
sigma_gain = max(0, double(cfg.region_sigma_threat_gain));

for k = 1:size(waypoints, 1)
    wp = double(squeeze(waypoints(k, :, 1:3)));
    if isempty(wp)
        continue;
    end
    [rows, cols] = xy_to_grid_indices(map, wp(:, 1), wp(:, 2));
    valid = rows >= 1 & rows <= H & cols >= 1 & cols <= W & ...
        isfinite(wp(:, 1)) & isfinite(wp(:, 2));
    valid_idx = find(valid(:)');
    for ii = valid_idx
        d_threat = double(threat_distance_m(rows(ii), cols(ii)));
        sigma_i = min(max(sigma_min + sigma_gain * d_threat, sigma_min), sigma_max);
        dist2 = (X_grid - wp(ii, 1)).^2 + (Y_grid - wp(ii, 2)).^2;
        contrib = single(exp(-dist2 ./ (2 * sigma_i^2)));
        heuristic = max(heuristic, contrib);
    end
end

heuristic = heuristic .* single(safe_mask);
max_h = max(heuristic(:));
if max_h > 0
    heuristic = heuristic ./ max_h;
end
heuristic = single(min(max(heuristic, 0), 1));
end

function safe_mask = build_region_safe_mask(map, sandwich, cost_field, z_layers_agl, cfg)
cost_hidden = min(double(cost_field), [], 3) <= double(cfg.visible_cost_threshold);
safe_mask = cost_hidden;
if cfg.region_use_sandwich_visibility && isfield(sandwich, 'min_visible_alt') && ~isempty(sandwich.min_visible_alt)
    terrain = double(map.elevation);
    sandwich_hidden_any = false(size(terrain));
    z_candidates = unique([double(cfg.region_safe_agl), double(z_layers_agl(:)')]);
    for i = 1:numel(z_candidates)
        z_asl = terrain + z_candidates(i);
        visible = z_asl - double(sandwich.min_visible_alt) + double(cfg.visibility_margin_m) > 0;
        sandwich_hidden_any = sandwich_hidden_any | ~visible;
    end
    safe_mask = safe_mask & sandwich_hidden_any;
end
safe_mask = single(safe_mask);
end

function [hidden, visible] = waypoint_visibility_masks(map, sandwich, cost_field, z_layers_agl, wp_agl, cfg)
costs = query_cost_field_agl(map, cost_field, z_layers_agl, wp_agl(:, 1), wp_agl(:, 2), wp_agl(:, 3));
cost_visible = costs > cfg.visible_cost_threshold;
if cfg.region_use_sandwich_visibility && isfield(sandwich, 'min_visible_alt') && ~isempty(sandwich.min_visible_alt)
    terrain = terrain_height(map, wp_agl(:, 1), wp_agl(:, 2));
    z_asl = terrain(:) + wp_agl(:, 3);
    min_visible_asl = interp_regular_2d(double(sandwich.min_visible_alt), ...
        double(map.xs), double(map.ys), wp_agl(:, 1), wp_agl(:, 2));
    sandwich_visible = z_asl(:) - min_visible_asl(:) + double(cfg.visibility_margin_m) > 0;
else
    sandwich_visible = false(size(cost_visible));
end
visible = cost_visible(:) | sandwich_visible(:);
hidden = ~visible;
hidden(~isfinite(wp_agl(:, 1)) | ~isfinite(wp_agl(:, 2)) | ~isfinite(wp_agl(:, 3))) = false;
visible(~isfinite(wp_agl(:, 1)) | ~isfinite(wp_agl(:, 2)) | ~isfinite(wp_agl(:, 3))) = false;
end

function [rows, cols] = xy_to_grid_indices(map, x, y)
dx = max(double(map.xs(2) - map.xs(1)), eps);
dy = max(double(map.ys(2) - map.ys(1)), eps);
cols = round(1 + (double(x(:)) - double(map.xs(1))) ./ dx);
rows = round(1 + (double(y(:)) - double(map.ys(1))) ./ dy);
end

function dist_m = distance_to_seed_m(seed, resolution_m)
if exist('bwdist', 'file') == 2
    dist_m = bwdist(seed) * resolution_m;
    return;
end
[H, W] = size(seed);
[rr, cc] = find(seed);
if isempty(rr)
    dist_m = inf(H, W);
    return;
end
[R, C] = ndgrid(1:H, 1:W);
dist2 = inf(H, W);
for i = 1:numel(rr)
    dist2 = min(dist2, (R - rr(i)).^2 + (C - cc(i)).^2);
end
dist_m = sqrt(dist2) * resolution_m;
end

function ctrl = fit_bspline_control_points(waypoints, n_ctrl)
K = size(waypoints, 1);
ctrl = zeros(K, n_ctrl, 3, 'single');
s = linspace(0, 1, size(waypoints, 2));
sq = linspace(0, 1, n_ctrl);
for k = 1:K
    for d = 1:3
        ctrl(k, :, d) = single(interp1(s, squeeze(waypoints(k, :, d)), sq, 'pchip'));
    end
end
end

function s = format_optional_idx(idx)
if isempty(idx)
    s = 'end';
else
    s = sprintf('%d', idx);
end
end

function n = read_existing_trajectory_count(filename)
n = NaN;
try
    info = h5info(filename, '/trajectories/waypoints');
    if ~isempty(info.Dataspace.Size)
        n = info.Dataspace.Size(1);
    end
catch
end
end

function save_generation_checkpoint(cfg, manifest)
if ~cfg.checkpoint_manifest
    return;
end
try
    checkpoint_file = cfg.manifest_file;
    save(checkpoint_file, 'manifest', 'cfg');
catch ME
    warning('generate_radar_uav_dataset:CheckpointFailed', ...
        'Could not save checkpoint manifest (%s).', ME.message);
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function cfg = configure_casadi(cfg)
cfg.casadi_available = false;
cfg.casadi_version = '';
cfg.local_solver_resolved = 'none';
if ~cfg.use_local_optimizer
    return;
end
if ~isempty(cfg.casadi_path) && exist(cfg.casadi_path, 'dir') == 7
    addpath(cfg.casadi_path);
end
try
    import casadi.*
    cfg.casadi_available = exist('casadi.Opti', 'class') == 8 || exist('casadi.MX', 'class') == 8;
    if cfg.casadi_available
        cfg.casadi_version = char(CasadiMeta.version());
        cfg.local_solver_resolved = resolve_casadi_solver(cfg.local_solver);
    end
catch ME
    warning('generate_radar_uav_dataset:CasadiUnavailable', ...
        'CasADi local optimizer disabled (%s).', ME.message);
    cfg.casadi_available = false;
    cfg.local_solver_resolved = 'none';
end
if cfg.use_local_optimizer && cfg.casadi_available
    fprintf('CasADi %s enabled, local_solver=%s\n', cfg.casadi_version, cfg.local_solver_resolved);
elseif cfg.use_local_optimizer
    fprintf('CasADi unavailable; using AIT*/minimum-snap fallback.\n');
end
end

function solver_name = resolve_casadi_solver(requested_solver)
requested_solver = lower(char(requested_solver));
if strcmp(requested_solver, 'auto')
    candidates = {'fatrop', 'ipopt'};
else
    candidates = {requested_solver};
end
solver_name = 'none';
for i = 1:numel(candidates)
    candidate = candidates{i};
    try
        import casadi.*
        x = MX.sym('x');
        nlp = struct('x', x, 'f', x^2);
        nlpsol('solver_probe', candidate, nlp);
        solver_name = candidate;
        return;
    catch
    end
end
if strcmp(requested_solver, 'auto')
    solver_name = 'ipopt';
else
    solver_name = requested_solver;
end
end

function cfg = configure_planner_backend(cfg)
cfg.planner_backend = lower(char(cfg.planner_backend));
if strcmp(cfg.planner_backend, 'ait_star_mex')
    if ~isempty(cfg.mex_path) && exist(cfg.mex_path, 'dir') == 7
        addpath(cfg.mex_path);
    end
    if isempty(which('ait_star_mex'))
        warning('generate_radar_uav_dataset:MexUnavailable', ...
            'ait_star_mex not found. Falling back to ait_star_matlab backend.');
        cfg.planner_backend = 'ait_star_matlab';
    end
elseif ~strcmp(cfg.planner_backend, 'ait_star_matlab')
    warning('generate_radar_uav_dataset:UnknownPlannerBackend', ...
        'Unknown planner_backend=%s. Falling back to ait_star_matlab.', cfg.planner_backend);
    cfg.planner_backend = 'ait_star_matlab';
end
end

function cleanup_obj = initialize_mex_planner_for_map(map, cost_field, z_layers_agl, cfg)
cleanup_obj = [];
if ~strcmpi(string(cfg.planner_backend), "ait_star_mex")
    return;
end
clear ait_star_mex;
ait_star_mex('init_map', double(map.elevation), double(map.resolution_m), ...
    double([map.xs(1), map.ys(1)]), double(cfg.min_clearance), ...
    double(cfg.mex_horizontal_clearance), 'dataset_map');
ait_star_mex('init_radar3d', single(cost_field), double(map.xs(:)'), double(map.ys(:)'), ...
    double(z_layers_agl(:)'), 'agl', double(cfg.planner_weight));
cleanup_obj = onCleanup(@() release_mex_planner_quietly());
end

function release_mex_planner_quietly()
try
    if ~isempty(which('ait_star_mex'))
        ait_star_mex('release');
    end
catch
end
end

function cfg = default_config(cfg, repo_root)
cfg.repo_root = get_cfg(cfg, 'repo_root', repo_root);
cfg.dem_dir = get_cfg(cfg, 'dem_dir', fullfile(repo_root, 'DEM'));
cfg.output_dir = get_cfg(cfg, 'output_dir', fullfile(repo_root, 'exports'));
cfg.base_seed = get_cfg(cfg, 'base_seed', 20260424);
cfg.total_maps = get_cfg(cfg, 'total_maps', 500);
cfg.map_limit = get_cfg(cfg, 'map_limit', []);
cfg.start_map_idx = get_cfg(cfg, 'start_map_idx', 0);
cfg.end_map_idx = get_cfg(cfg, 'end_map_idx', []);
cfg.skip_existing = get_cfg(cfg, 'skip_existing', true);
cfg.write_partial_maps = get_cfg(cfg, 'write_partial_maps', true);
cfg.checkpoint_manifest = get_cfg(cfg, 'checkpoint_manifest', true);
cfg.manifest_file = get_cfg(cfg, 'manifest_file', fullfile(cfg.output_dir, 'generation_manifest.mat'));
cfg.K = get_cfg(cfg, 'K', 20);
cfg.N_samples = get_cfg(cfg, 'N_samples', 64);
cfg.dt = get_cfg(cfg, 'dt', 0.05);
cfg.dt_mode = get_cfg(cfg, 'dt_mode', 'adaptive');
cfg.nominal_speed_mps = get_cfg(cfg, 'nominal_speed_mps', 28);
cfg.min_dt = get_cfg(cfg, 'min_dt', 0.05);
cfg.max_dt = get_cfg(cfg, 'max_dt', 1.25);
cfg.dt_safety_factor = get_cfg(cfg, 'dt_safety_factor', 1.20);
cfg.export_time_channel = get_cfg(cfg, 'export_time_channel', true);
cfg.dem_size = get_cfg(cfg, 'dem_size', [256, 256]);
cfg.map_resolution_m = get_cfg(cfg, 'map_resolution_m', 5.0);
cfg.real_crop_mode = get_cfg(cfg, 'real_crop_mode', 'random');
cfg.synthetic_types = get_cfg(cfg, 'synthetic_types', {'ridge', 'valley', 'mountain', 'hills', 'random'});
cfg.synthetic_resolution_m = get_cfg(cfg, 'synthetic_resolution_m', cfg.map_resolution_m);
cfg.synthetic_amplitude = get_cfg(cfg, 'synthetic_amplitude', 85);

cfg.num_radars_range = get_cfg(cfg, 'num_radars_range', [1, 2]);
cfg.radar_height_agl_range = get_cfg(cfg, 'radar_height_agl_range', [8, 40]);
cfg.z_layers_agl = get_cfg(cfg, 'z_layers_agl', linspace(8, 180, 14));
cfg.ray_samples = get_cfg(cfg, 'ray_samples', 28);
cfg.ray_chunk_size = get_cfg(cfg, 'ray_chunk_size', 512);
cfg.ray_pair_batch_size = get_cfg(cfg, 'ray_pair_batch_size', 250000);
cfg.use_sandwich_mode = get_cfg(cfg, 'use_sandwich_mode', true);
cfg.sandwich_z_agl_layers = get_cfg(cfg, 'sandwich_z_agl_layers', linspace(5, 220, 24));
cfg.sandwich_default_z_agl = get_cfg(cfg, 'sandwich_default_z_agl', max(cfg.sandwich_z_agl_layers));
cfg.sandwich_stride = get_cfg(cfg, 'sandwich_stride', max(1, round(max(cfg.dem_size) / 64)));
cfg.sigmoid_k = get_cfg(cfg, 'sigmoid_k', 0.12);
cfg.radar_cost_gain = get_cfg(cfg, 'radar_cost_gain', 2.5e8);
cfg.radar_cost_normalizer = get_cfg(cfg, 'radar_cost_normalizer', 1.0);
cfg.radar_Pt = get_cfg(cfg, 'radar_Pt', 1000);
cfg.radar_G = get_cfg(cfg, 'radar_G', 1e3);
cfg.radar_lambda = get_cfg(cfg, 'radar_lambda', 0.03);
cfg.radar_sigma = get_cfg(cfg, 'radar_sigma', 0.01);
cfg.radar_kB = get_cfg(cfg, 'radar_kB', 1.38e-23);
cfg.radar_T0 = get_cfg(cfg, 'radar_T0', 290);
cfg.radar_B = get_cfg(cfg, 'radar_B', 1e6);
cfg.radar_F = get_cfg(cfg, 'radar_F', 10^(3/10));
cfg.radar_snr_constant = get_cfg(cfg, 'radar_snr_constant', ...
    (cfg.radar_Pt * cfg.radar_G^2 * cfg.radar_lambda^2 * cfg.radar_sigma) / ...
    (((4*pi)^3) * cfg.radar_kB * cfg.radar_T0 * cfg.radar_B * cfg.radar_F));
cfg.radar_effective_range_m = get_cfg(cfg, 'radar_effective_range_m', NaN);
cfg.radar_effective_range_frac = get_cfg(cfg, 'radar_effective_range_frac', 0.55);
cfg.radar_range_softness_m = get_cfg(cfg, 'radar_range_softness_m', NaN);
cfg.radar_range_softness_frac = get_cfg(cfg, 'radar_range_softness_frac', 0.08);

cfg.start_goal_margin_m = get_cfg(cfg, 'start_goal_margin_m', 30);
cfg.start_goal_agl_range = get_cfg(cfg, 'start_goal_agl_range', [20, 120]);
cfg.start_goal_mode = get_cfg(cfg, 'start_goal_mode', 'different_quadrants');
cfg.start_goal_corner_jitter_frac = get_cfg(cfg, 'start_goal_corner_jitter_frac', 0.12);
cfg.scenario_sampling_mode = get_cfg(cfg, 'scenario_sampling_mode', 'radar_relevant');
cfg.scenario_mix = get_cfg(cfg, 'scenario_mix', [0.15, 0.50, 0.35]); % easy, medium, hard
cfg.enforce_scenario_quota = get_cfg(cfg, 'enforce_scenario_quota', true);
cfg.visible_cost_threshold = get_cfg(cfg, 'visible_cost_threshold', 0.05);
cfg.easy_straight_visible_range = get_cfg(cfg, 'easy_straight_visible_range', [0.00, 0.25]);
cfg.medium_straight_visible_range = get_cfg(cfg, 'medium_straight_visible_range', [0.25, 0.70]);
cfg.hard_straight_visible_range = get_cfg(cfg, 'hard_straight_visible_range', [0.60, 1.00]);
cfg.max_planned_visible_ratio = get_cfg(cfg, 'max_planned_visible_ratio', 0.35);
cfg.max_planned_mean_radar_cost = get_cfg(cfg, 'max_planned_mean_radar_cost', 0.20);
cfg.min_radar_improvement_ratio = get_cfg(cfg, 'min_radar_improvement_ratio', 0.30);
cfg.enforce_visibility_acceptance = get_cfg(cfg, 'enforce_visibility_acceptance', true);
cfg.visibility_margin_m = get_cfg(cfg, 'visibility_margin_m', 0.1);
cfg.max_planned_visibility_violations = get_cfg(cfg, 'max_planned_visibility_violations', [12, 8, 5]); % easy, medium, hard
cfg.label_mode = get_cfg(cfg, 'label_mode', 'trajectory_and_region');
if strcmpi(string(cfg.label_mode), "trajectory_and_region")
    cfg.region_label_type = get_cfg(cfg, 'region_label_type', '2d_corridor');
else
    cfg.region_label_type = get_cfg(cfg, 'region_label_type', 'none');
end
cfg.radar_constraint_mode = get_cfg(cfg, 'radar_constraint_mode', 'soft_cost');
cfg.visibility_budget_schedule = get_cfg(cfg, 'visibility_budget_schedule', [0, 2, 5, 10, 20]);
cfg.relaxed_target_weight = get_cfg(cfg, 'relaxed_target_weight', 0.35);
cfg.min_violation_target_weight = get_cfg(cfg, 'min_violation_target_weight', 0.15);
cfg.region_corridor_sigma_m = get_cfg(cfg, 'region_corridor_sigma_m', 250);
cfg.region_visible_exclusion_sigma_m = get_cfg(cfg, 'region_visible_exclusion_sigma_m', 150);
cfg.region_corridor_mode = get_cfg(cfg, 'region_corridor_mode', 'threat_distance_adaptive');
cfg.region_sigma_min_m = get_cfg(cfg, 'region_sigma_min_m', 60);
cfg.region_sigma_max_m = get_cfg(cfg, 'region_sigma_max_m', 300);
cfg.region_sigma_threat_gain = get_cfg(cfg, 'region_sigma_threat_gain', 0.5);
cfg.region_safe_agl = get_cfg(cfg, 'region_safe_agl', median(double(cfg.z_layers_agl(:))));
cfg.region_use_sandwich_visibility = get_cfg(cfg, 'region_use_sandwich_visibility', true);
cfg.radar_relevant_sampling_attempts = get_cfg(cfg, 'radar_relevant_sampling_attempts', 500);
cfg.enforce_start_goal_reachability = get_cfg(cfg, 'enforce_start_goal_reachability', true);
cfg.min_start_goal_dist_frac = get_cfg(cfg, 'min_start_goal_dist_frac', 0.35);
cfg.v_max = get_cfg(cfg, 'v_max', 35);
cfg.a_max = get_cfg(cfg, 'a_max', 12);
cfg.reachability_distance_factor = get_cfg(cfg, 'reachability_distance_factor', 0.70);
if strcmpi(string(cfg.dt_mode), "adaptive")
    default_horizon_s = max((cfg.N_samples - 1) * cfg.max_dt, cfg.max_dt);
else
    default_horizon_s = max((cfg.N_samples - 1) * cfg.dt, cfg.dt);
end
cfg.max_start_goal_dist_m = get_cfg(cfg, 'max_start_goal_dist_m', cfg.reachability_distance_factor * cfg.v_max * default_horizon_s);
cfg.min_clearance = get_cfg(cfg, 'min_clearance', 8);
cfg.max_clearance = get_cfg(cfg, 'max_clearance', 180);
cfg.validation_min_agl = get_cfg(cfg, 'validation_min_agl', 5);
cfg.max_path_attempts = get_cfg(cfg, 'max_path_attempts', max(12 * cfg.K, 300));
cfg.trajectory_batch_size = get_cfg(cfg, 'trajectory_batch_size', 8);
cfg.trajectory_attempt_oversample_factor = get_cfg(cfg, 'trajectory_attempt_oversample_factor', 1);
cfg.log_rejection_summary = get_cfg(cfg, 'log_rejection_summary', true);
cfg.log_rejection_examples = get_cfg(cfg, 'log_rejection_examples', 3);

cfg.planner_weight = get_cfg(cfg, 'planner_weight', 5.0);
cfg.planner_altitude_weight = get_cfg(cfg, 'planner_altitude_weight', 0.0);
cfg.planner_preferred_agl = get_cfg(cfg, 'planner_preferred_agl', median(double(cfg.start_goal_agl_range(:))));
cfg.planner_max_time = get_cfg(cfg, 'planner_max_time', 4);
cfg.planner_max_batches = get_cfg(cfg, 'planner_max_batches', 18);
cfg.planner_batch_size = get_cfg(cfg, 'planner_batch_size', 160);
cfg.planner_max_nodes = get_cfg(cfg, 'planner_max_nodes', 3500);
cfg.planner_post_solution_batches = get_cfg(cfg, 'planner_post_solution_batches', 2);
cfg.scenario_planner_time_scale = get_cfg(cfg, 'scenario_planner_time_scale', [1.0, 1.8, 3.0]); % easy, medium, hard
cfg.scenario_planner_batch_scale = get_cfg(cfg, 'scenario_planner_batch_scale', [1.0, 1.6, 2.6]);
cfg.scenario_planner_node_scale = get_cfg(cfg, 'scenario_planner_node_scale', [1.0, 2.0, 4.0]);
cfg.scenario_planner_post_solution_scale = get_cfg(cfg, 'scenario_planner_post_solution_scale', [1.0, 1.5, 2.5]);
cfg.edge_check_samples = get_cfg(cfg, 'edge_check_samples', 9);
cfg.goal_sample_fraction = get_cfg(cfg, 'goal_sample_fraction', 0.22);
cfg.line_sample_fraction = get_cfg(cfg, 'line_sample_fraction', 0.30);
cfg.goal_probe_count = get_cfg(cfg, 'goal_probe_count', 10);
cfg.reverse_k_neighbors = get_cfg(cfg, 'reverse_k_neighbors', 22);
cfg.forward_k_neighbors = get_cfg(cfg, 'forward_k_neighbors', 16);
cfg.reverse_batch_repair_limit = get_cfg(cfg, 'reverse_batch_repair_limit', 128);
cfg.forward_expand_chunk = get_cfg(cfg, 'forward_expand_chunk', 32);
cfg.use_parallel_ait_reverse = get_cfg(cfg, 'use_parallel_ait_reverse', true);
cfg.use_parallel_ait_edge_checks = get_cfg(cfg, 'use_parallel_ait_edge_checks', true);
cfg.graph_rebuild_interval = get_cfg(cfg, 'graph_rebuild_interval', 3);
cfg.reverse_full_rebuild_interval = get_cfg(cfg, 'reverse_full_rebuild_interval', cfg.graph_rebuild_interval);

cfg.use_optimized_time_allocation = get_cfg(cfg, 'use_optimized_time_allocation', true);
cfg.smooth_control_points = get_cfg(cfg, 'smooth_control_points', 7);
cfg.smooth_v_max_range = get_cfg(cfg, 'smooth_v_max_range', [15, 25]);
cfg.smooth_a_max_range = get_cfg(cfg, 'smooth_a_max_range', [6, 15]);
cfg.smooth_cruise_speed_range = get_cfg(cfg, 'smooth_cruise_speed_range', [10, 18]);
cfg.smooth_aggressiveness_range = get_cfg(cfg, 'smooth_aggressiveness_range', [3, 5]);
cfg.smooth_v_safety_range = get_cfg(cfg, 'smooth_v_safety_range', [1.15, 1.35]);
cfg.smooth_samples_per_seg = get_cfg(cfg, 'smooth_samples_per_seg', 50);
cfg.smooth_vel_bc_mode = get_cfg(cfg, 'smooth_vel_bc_mode', 'free');
cfg.smooth_verbose = get_cfg(cfg, 'smooth_verbose', false);

cfg.use_local_optimizer = get_cfg(cfg, 'use_local_optimizer', false);
cfg.use_local_optimizer_template = get_cfg(cfg, 'use_local_optimizer_template', true);
cfg.max_radars = get_cfg(cfg, 'max_radars', 5);
cfg.casadi_path = get_cfg(cfg, 'casadi_path', 'D:\casadi-3.7.2-windows64-matlab2018b');
cfg.local_solver = get_cfg(cfg, 'local_solver', 'fatrop');
cfg.local_optimizer_max_iter = get_cfg(cfg, 'local_optimizer_max_iter', 80);
cfg.local_optimizer_time_limit = get_cfg(cfg, 'local_optimizer_time_limit', 0.2);
cfg.local_tracking_weight = get_cfg(cfg, 'local_tracking_weight', 0.04);
cfg.local_smooth_weight = get_cfg(cfg, 'local_smooth_weight', 1e-3);
cfg.local_control_weight = get_cfg(cfg, 'local_control_weight', 1e-4);
cfg.local_radar_weight = get_cfg(cfg, 'local_radar_weight', 1.0);
cfg.local_use_mex = get_cfg(cfg, 'local_use_mex', false);
cfg.export_bspline = get_cfg(cfg, 'export_bspline', false);
cfg.bspline_control_points = get_cfg(cfg, 'bspline_control_points', 16);

cfg.planner_backend = get_cfg(cfg, 'planner_backend', 'ait_star_mex');
cfg.mex_path = get_cfg(cfg, 'mex_path', fullfile(repo_root, 'ait_star', 'build', 'Release'));
cfg.mex_horizontal_clearance = get_cfg(cfg, 'mex_horizontal_clearance', 50);
cfg.generate_trajectory = get_cfg(cfg, 'generate_trajectory', false);
if ~cfg.generate_trajectory
    cfg.use_local_optimizer = false;
end

cfg.use_parallel = get_cfg(cfg, 'use_parallel', true);
cfg.parallel_pool = get_cfg(cfg, 'parallel_pool', 'threads');
cfg.use_parallel_cost_field = get_cfg(cfg, 'use_parallel_cost_field', true);
cfg.use_parallel_trajectories = get_cfg(cfg, 'use_parallel_trajectories', true);
cfg.use_gpu = get_cfg(cfg, 'use_gpu', false);
cfg.use_gpu_cost_field = get_cfg(cfg, 'use_gpu_cost_field', cfg.use_gpu);
cfg.gpu_device_index = get_cfg(cfg, 'gpu_device_index', 1);
cfg.gpu_gather_each_layer = get_cfg(cfg, 'gpu_gather_each_layer', true);
end

function map_specs = build_map_specs(cfg)
dem_files = dir(fullfile(cfg.dem_dir, '*.tif'));
num_real = numel(dem_files);
num_maps = max(cfg.total_maps, num_real);
map_specs = repmat(struct('source', '', 'path', '', 'synthetic', false, 'synthetic_type', ''), 1, num_maps);
for i = 1:num_real
    map_specs(i).source = dem_files(i).name;
    map_specs(i).path = fullfile(dem_files(i).folder, dem_files(i).name);
    map_specs(i).synthetic = false;
end
for i = (num_real + 1):num_maps
    type = cfg.synthetic_types{1 + mod(i - num_real - 1, numel(cfg.synthetic_types))};
    map_specs(i).source = sprintf('synthetic_%s_%02d', type, i - num_real);
    map_specs(i).synthetic = true;
    map_specs(i).synthetic_type = type;
end
if ~isempty(cfg.map_limit)
    map_specs = map_specs(1:min(numel(map_specs), cfg.map_limit));
end
end

function map = load_or_create_map(spec, cfg, seed)
if spec.synthetic
    map = create_synthetic_map(spec.synthetic_type, cfg, seed);
else
    map = load_geotiff_map(spec.path, cfg, seed);
end
end

function map = load_geotiff_map(path, cfg, seed)
try
    map = load_geotiff_map_with_dem_loader(path, cfg, seed);
    return;
catch ME
    warning('generate_radar_uav_dataset:DEMLoaderFallback', ...
        'dem_loader failed for %s (%s). Falling back to direct readgeoraster.', path, ME.message);
end

[Z_raw, R] = readgeoraster(path);
if size(Z_raw, 3) > 1
    Z_raw = Z_raw(:, :, 1);
end
Z = double(Z_raw);
Z(Z < -1000 | Z > 10000) = NaN;
Z = fill_nearest_nan(Z);

[res_x, res_y] = estimate_resolution_m(R, size(Z));
res_native = mean([res_x, res_y], 'omitnan');
if ~isfinite(res_native) || res_native <= 0
    res_native = 30;
end

target_h = cfg.dem_size(1);
target_w = cfg.dem_size(2);
Z_crop = crop_and_resample_dem(Z, res_native, cfg, seed);
resolution_m = single(cfg.map_resolution_m);

xs = ((0:target_w-1) - (target_w - 1) / 2) * double(resolution_m);
ys = ((0:target_h-1) - (target_h - 1) / 2) * double(resolution_m);

map = struct();
map.xs = xs;
map.ys = ys;
map.elevation = single(Z_crop);
map.resolution_m = single(resolution_m);
map.origin_xy = single([xs(1), ys(1)]);
map.source = path;
end

function map = load_geotiff_map_with_dem_loader(path, cfg, seed)
dem_params = struct();
dem_params.fill_nodata = 'nearest';
terrain_data = dem_loader(path, dem_params);
map = terrain_data_to_centered_map(terrain_data, cfg, path, seed);
end

function map = terrain_data_to_centered_map(terrain_data, cfg, source_path, seed)
E_vec = double(terrain_data.E_vec(:)');
N_vec = double(terrain_data.N_vec(:)');
Z = double(terrain_data.Z);

if isequal(size(Z), [numel(E_vec), numel(N_vec)])
    Z_ne = Z.';
elseif isequal(size(Z), [numel(N_vec), numel(E_vec)])
    Z_ne = Z;
else
    error('generate_radar_uav_dataset:DEMSizeMismatch', ...
        'DEM Z size [%d %d] does not match E/N vectors [%d %d].', ...
        size(Z, 1), size(Z, 2), numel(E_vec), numel(N_vec));
end

target_h = cfg.dem_size(1);
target_w = cfg.dem_size(2);
res_e = mean(abs(diff(E_vec)), 'omitnan');
res_n = mean(abs(diff(N_vec)), 'omitnan');
res_native = mean([res_e, res_n], 'omitnan');
if ~isfinite(res_native) || res_native <= 0
    res_native = double(terrain_data.resolution);
end
Z_resized = crop_and_resample_dem(fill_nearest_nan(Z_ne), res_native, cfg, seed);
resolution_m = single(cfg.map_resolution_m);

xs = ((0:target_w-1) - (target_w - 1) / 2) * double(resolution_m);
ys = ((0:target_h-1) - (target_h - 1) / 2) * double(resolution_m);

map = struct();
map.xs = xs;
map.ys = ys;
map.elevation = single(Z_resized);
map.resolution_m = single(resolution_m);
map.origin_xy = single([xs(1), ys(1)]);
map.source = source_path;
end

function [res_x, res_y] = estimate_resolution_m(R, raster_size)
res_x = NaN;
res_y = NaN;
if isprop(R, 'CellExtentInWorldX') && isprop(R, 'CellExtentInWorldY')
    res_x = abs(R.CellExtentInWorldX);
    res_y = abs(R.CellExtentInWorldY);
elseif isprop(R, 'SampleSpacingInWorldX') && isprop(R, 'SampleSpacingInWorldY')
    res_x = abs(R.SampleSpacingInWorldX);
    res_y = abs(R.SampleSpacingInWorldY);
end

if isprop(R, 'LatitudeLimits') && isprop(R, 'LongitudeLimits')
    lat_limits = R.LatitudeLimits;
    lon_limits = R.LongitudeLimits;
    center_lat = mean(lat_limits);
    meters_per_deg_lat = 1000 * deg2km_local(1);
    meters_per_deg_lon = meters_per_deg_lat * cosd(center_lat);
    res_y = abs(diff(lat_limits)) / max(raster_size(1), 1) * meters_per_deg_lat;
    res_x = abs(diff(lon_limits)) / max(raster_size(2), 1) * meters_per_deg_lon;
elseif isprop(R, 'XWorldLimits') && isprop(R, 'YWorldLimits')
    res_x = abs(diff(R.XWorldLimits)) / max(raster_size(2), 1);
    res_y = abs(diff(R.YWorldLimits)) / max(raster_size(1), 1);
end
end

function km = deg2km_local(deg)
km = 111.32 * deg;
end

function map = create_synthetic_map(type, cfg, seed)
rng(seed, 'twister');
H = cfg.dem_size(1);
W = cfg.dem_size(2);
res = cfg.synthetic_resolution_m;
xs = ((0:W-1) - (W - 1) / 2) * res;
ys = ((0:H-1) - (H - 1) / 2) * res;
[X, Y] = meshgrid(xs, ys);
amp = cfg.synthetic_amplitude;
base = 120 + 80 * rand();
noise = fractal_value_noise(H, W, 6, 0.52);
noise = (noise - mean(noise(:))) ./ max(std(noise(:)), eps);
ridge_noise = 1 - abs(max(min(noise / 2.5, 1), -1));
ridge_noise = ridge_noise - mean(ridge_noise(:));
map_span = max(range(xs), range(ys));

switch lower(type)
    case 'ridge'
        shape = exp(-(Y .^ 2) / (2 * (0.18 * max(abs(ys)))^2)) + 0.35 * ridge_noise;
    case 'valley'
        shape = 0.9 * abs(Y) / max(abs(ys)) - 0.55 * exp(-(Y .^ 2) / (2 * (0.16 * max(abs(ys)))^2)) + 0.35 * noise;
    case 'mountain'
        shape = exp(-(X.^2 + Y.^2) / (2 * (0.24 * max(abs(xs)))^2)) + 0.45 * noise;
    case 'hills'
        shape = 0.45 * sin(2*pi*X / max(map_span, 1)) .* cos(2*pi*Y / max(map_span, 1)) + 0.55 * noise;
    otherwise
        shape = 0.65 * noise + 0.35 * ridge_noise;
end
Z = base + amp * shape;
Z = smoothdata(smoothdata(Z, 1, 'movmean', 5), 2, 'movmean', 5);
Z = max(Z, 0);

map = struct();
map.xs = xs;
map.ys = ys;
map.elevation = single(Z);
map.resolution_m = single(res);
map.origin_xy = single([xs(1), ys(1)]);
map.source = ['synthetic_', type];
end

function Z = fill_nearest_nan(Z)
mask = isfinite(Z);
if all(mask(:))
    return;
end
if ~any(mask(:))
    Z(:) = 0;
    return;
end
if exist('bwdist', 'file') == 2
    [~, idx] = bwdist(mask);
    Z(~mask) = Z(idx(~mask));
else
    Z = fillmissing(Z, 'nearest', 1);
    Z = fillmissing(Z, 'nearest', 2);
    Z(~isfinite(Z)) = mean(Z(isfinite(Z)), 'omitnan');
end
end

function Zout = resize_bilinear(Z, target_size)
if exist('imresize', 'file') == 2
    Zout = imresize(Z, target_size, 'bilinear');
    return;
end
[h, w] = size(Z);
[X, Y] = meshgrid(1:w, 1:h);
xq = linspace(1, w, target_size(2));
yq = linspace(1, h, target_size(1));
[Xq, Yq] = meshgrid(xq, yq);
Zout = interp2(X, Y, Z, Xq, Yq, 'linear');
Zout = fill_nearest_nan(Zout);
end

function Zout = crop_and_resample_dem(Z, native_resolution_m, cfg, seed)
target_h = cfg.dem_size(1);
target_w = cfg.dem_size(2);
target_res = double(cfg.map_resolution_m);
native_resolution_m = double(native_resolution_m);
if ~isfinite(native_resolution_m) || native_resolution_m <= 0
    native_resolution_m = target_res;
end
crop_h = min(size(Z, 1), max(2, round(((target_h - 1) * target_res) / native_resolution_m) + 1));
crop_w = min(size(Z, 2), max(2, round(((target_w - 1) * target_res) / native_resolution_m) + 1));
if strcmpi(cfg.real_crop_mode, 'center') || size(Z, 1) == crop_h
    r0 = floor((size(Z, 1) - crop_h) / 2) + 1;
else
    rng(seed + 17, 'twister');
    r0 = randi(size(Z, 1) - crop_h + 1);
end
if strcmpi(cfg.real_crop_mode, 'center') || size(Z, 2) == crop_w
    c0 = floor((size(Z, 2) - crop_w) / 2) + 1;
else
    rng(seed + 31, 'twister');
    c0 = randi(size(Z, 2) - crop_w + 1);
end
Z_crop = Z(r0:(r0 + crop_h - 1), c0:(c0 + crop_w - 1));
Zout = resize_makima(Z_crop, [target_h, target_w]);
end

function Zout = resize_makima(Z, target_size)
[h, w] = size(Z);
[X, Y] = meshgrid(1:w, 1:h);
xq = linspace(1, w, target_size(2));
yq = linspace(1, h, target_size(1));
[Xq, Yq] = meshgrid(xq, yq);
try
    Zout = interp2(X, Y, Z, Xq, Yq, 'makima');
catch
    Zout = interp2(X, Y, Z, Xq, Yq, 'linear');
end
Zout = fill_nearest_nan(Zout);
end

function N = fractal_value_noise(H, W, octaves, persistence)
N = zeros(H, W);
amp = 1;
amp_sum = 0;
for octave = 1:octaves
    gh = max(3, ceil(H / (2 ^ (octaves - octave + 1))));
    gw = max(3, ceil(W / (2 ^ (octaves - octave + 1))));
    coarse = randn(gh, gw);
    layer = resize_makima(coarse, [H, W]);
    N = N + amp * layer;
    amp_sum = amp_sum + amp;
    amp = amp * persistence;
end
N = N ./ max(amp_sum, eps);
end

function radars = scatter_radars_enu(map, cfg)
R = randi(cfg.num_radars_range);
radars = zeros(R, 3);
x_min = map.xs(1); x_max = map.xs(end);
y_min = map.ys(1); y_max = map.ys(end);
for r = 1:R
    x = x_min + rand() * (x_max - x_min);
    y = y_min + rand() * (y_max - y_min);
    terrain = terrain_height(map, x, y);
    z = terrain + rand_range(cfg.radar_height_agl_range);
    radars(r, :) = [x, y, z];
end
end

function sandwich = build_sandwich_maps(map, radars, cfg)
H = numel(map.ys);
W = numel(map.xs);
R = size(radars, 1);
sandwich = struct();
sandwich.z_min_by_radar = repmat(single(cfg.sandwich_default_z_agl), H, W, R);
sandwich.min_visible_alt = single(map.elevation) + single(cfg.sandwich_default_z_agl);
sandwich.thickness = repmat(single(cfg.sandwich_default_z_agl), H, W);
if ~cfg.use_sandwich_mode || R == 0
    return;
end

mesh = build_numeric_terrain_mesh(map);
cost_ctx = build_cost_context(map);
z_layers = double(cfg.sandwich_z_agl_layers(:)');
stride = max(1, round(cfg.sandwich_stride));
rows = unique([1:stride:H, H]);
cols = unique([1:stride:W, W]);
[Cq, Rq] = meshgrid(cols, rows);
query_idx = sub2ind([H, W], Rq(:).', Cq(:).');
Hq = numel(rows);
Wq = numel(cols);
fprintf('Computing sandwich min-visible-alt maps: %d x %d query grid -> %d x %d full, %d radars x %d layers\n', ...
    Hq, Wq, H, W, R, numel(z_layers));
t_sandwich = tic;
for r = 1:R
    radar = radars(r, :).';
    z_min = inf(1, numel(query_idx));
    unresolved = true(1, numel(query_idx));
    for k = 1:numel(z_layers)
        if ~any(unresolved)
            break;
        end
        idx = find(unresolved);
        map_idx = query_idx(idx);
        pts = [cost_ctx.x(map_idx); cost_ctx.y(map_idx); cost_ctx.terrain(map_idx) + z_layers(k)];
        visible = ~ray_tracing(mesh, radar, pts, cfg.ray_chunk_size, cfg.ray_pair_batch_size);
        if any(visible)
            hit_idx = idx(visible);
            z_min(hit_idx) = z_layers(k);
            unresolved(hit_idx) = false;
        end
    end
    z_min(~isfinite(z_min)) = cfg.sandwich_default_z_agl;
    z_min_full = resize_makima(reshape(z_min, Hq, Wq), [H, W]);
    sandwich.z_min_by_radar(:, :, r) = single(z_min_full);
end
sandwich.thickness = min(sandwich.z_min_by_radar, [], 3);
sandwich.min_visible_alt = single(map.elevation) + sandwich.thickness;
fprintf('  sandwich: done in %.1fs\n', toc(t_sandwich));
end

function [cost_field, z_layers_agl] = build_radar_cost_field_agl(map, radars, sandwich, cfg)
z_layers_agl = single(cfg.z_layers_agl(:)');
H = numel(map.ys);
W = numel(map.xs);
ZL = numel(z_layers_agl);
cost_field = zeros(H, W, ZL, 'single');
mesh = build_numeric_terrain_mesh(map);
cost_ctx = build_cost_context(map);
cfg_cost = cfg;
diag_len = hypot(range(map.xs), range(map.ys));
if ~isfinite(cfg_cost.radar_effective_range_m) || cfg_cost.radar_effective_range_m <= 0
    cfg_cost.radar_effective_range_m = cfg_cost.radar_effective_range_frac * diag_len;
end
if ~isfinite(cfg_cost.radar_range_softness_m) || cfg_cost.radar_range_softness_m <= 0
    cfg_cost.radar_range_softness_m = max(cfg_cost.radar_range_softness_frac * diag_len, double(map.resolution_m));
end
fprintf('Computing AGL radar cost field with vectorized Moller-Trumbore: %d x %d x %d\n', H, W, ZL);
fprintf('  radar coverage: effective_range=%.1fm softness=%.1fm\n', ...
    cfg_cost.radar_effective_range_m, cfg_cost.radar_range_softness_m);
cost_tic = tic;

gpu_status = maybe_select_gpu(cfg_cost);
use_gpu_cost = gpu_status.enabled && cfg_cost.use_sandwich_mode && isfield(sandwich, 'z_min_by_radar');
if gpu_status.enabled && ~use_gpu_cost
    warning('generate_radar_uav_dataset:GpuCostFieldFallback', ...
        'GPU cost field is only used for sandwich analytical mode. Falling back to CPU cost field.');
end

if use_gpu_cost
    fprintf('  cost field: GPU over %d altitude layers on %s\n', ZL, gpu_status.name);
    gpu_ctx = build_gpu_cost_context(cost_ctx, sandwich);
    if cfg_cost.gpu_gather_each_layer
        for k = 1:ZL
            cost_field(:, :, k) = compute_cost_layer_gpu(gpu_ctx, radars, double(z_layers_agl(k)), cfg_cost, true);
        end
    else
        gpu_cost_field = gpuArray.zeros(H, W, ZL, 'single');
        for k = 1:ZL
            gpu_cost_field(:, :, k) = compute_cost_layer_gpu(gpu_ctx, radars, double(z_layers_agl(k)), cfg_cost, false);
        end
        cost_field = gather(gpu_cost_field);
    end
else
    use_par_cost = cfg.use_parallel && cfg.use_parallel_cost_field && can_use_parfor() && ~isempty(gcp('nocreate'));
    if use_par_cost
        pool = gcp('nocreate');
        fprintf('  cost field: parfor over %d altitude layers on %d workers\n', ZL, pool.NumWorkers);
        parfor k = 1:ZL
            cost_field(:, :, k) = compute_cost_layer(cost_ctx, mesh, radars, sandwich, double(z_layers_agl(k)), cfg_cost);
        end
    else
        fprintf('  cost field: serial loop over %d altitude layers\n', ZL);
        for k = 1:ZL
            cost_field(:, :, k) = compute_cost_layer(cost_ctx, mesh, radars, sandwich, double(z_layers_agl(k)), cfg_cost);
        end
    end
end
fprintf('  cost field: done in %.1fs\n', toc(cost_tic));
end

function ctx = build_cost_context(map)
[X, Y] = meshgrid(double(map.xs), double(map.ys));
ctx.H = numel(map.ys);
ctx.W = numel(map.xs);
ctx.x = X(:).';
ctx.y = Y(:).';
ctx.terrain = double(map.elevation(:)).';
end

function gpu_ctx = build_gpu_cost_context(cost_ctx, sandwich)
gpu_ctx = struct();
gpu_ctx.H = cost_ctx.H;
gpu_ctx.W = cost_ctx.W;
gpu_ctx.x = gpuArray(cost_ctx.x);
gpu_ctx.y = gpuArray(cost_ctx.y);
gpu_ctx.terrain = gpuArray(cost_ctx.terrain);
gpu_ctx.z_min_by_radar = gpuArray(single(sandwich.z_min_by_radar));
end

function layer = compute_cost_layer(cost_ctx, mesh, radars, sandwich, z_agl, cfg)
layer = zeros(cost_ctx.H, cost_ctx.W, 'single');
z_asl = cost_ctx.terrain + z_agl;
pts = [cost_ctx.x; cost_ctx.y; z_asl];
for r = 1:size(radars, 1)
    radar = radars(r, :).';
    if cfg.use_sandwich_mode && isfield(sandwich, 'z_min_by_radar')
        z_min_agl = reshape(double(sandwich.z_min_by_radar(:, :, r)), 1, []);
        visible = z_agl >= z_min_agl;
        snr_cost = analytical_sigmoid_radar_cost(radar, pts, z_asl, cost_ctx.terrain + z_min_agl, cfg);
    else
        visible = ~ray_tracing(mesh, radar, pts, cfg.ray_chunk_size, cfg.ray_pair_batch_size);
        snr_cost = radar_snr_cost(radar, pts, cfg);
    end
    vals = zeros(1, numel(visible), 'single');
    vals(visible) = snr_cost(visible);
    layer = max(layer, reshape(vals, cost_ctx.H, cost_ctx.W));
end
end

function layer = compute_cost_layer_gpu(gpu_ctx, radars, z_agl, cfg, gather_output)
if nargin < 5
    gather_output = true;
end
layer_gpu = gpuArray.zeros(gpu_ctx.H, gpu_ctx.W, 'single');
z_asl = gpu_ctx.terrain + z_agl;
for r = 1:size(radars, 1)
    radar = gpuArray(double(radars(r, :).'));
    z_min_agl = reshape(double(gpu_ctx.z_min_by_radar(:, :, r)), 1, []);
    visible = z_agl >= z_min_agl;

    dx = gpu_ctx.x - radar(1);
    dy = gpu_ctx.y - radar(2);
    dz = z_asl - radar(3);
    d2 = max(dx.^2 + dy.^2 + dz.^2, 1);
    d = sqrt(d2);

    z_min_asl = gpu_ctx.terrain + z_min_agl;
    sigmoid = 1 ./ (1 + exp(-cfg.sigmoid_k .* (z_asl - z_min_asl)));
    if isfield(cfg, 'radar_effective_range_m') && isfinite(cfg.radar_effective_range_m) && cfg.radar_effective_range_m > 0
        softness = max(cfg.radar_range_softness_m, 1);
        range_gate = 1 ./ (1 + exp((d - cfg.radar_effective_range_m) ./ softness));
    else
        range_gate = 1;
    end
    raw = range_gate .* sigmoid .* cfg.radar_cost_gain ./ (d2 .^ 2);
    snr_cost = single(min(max(raw ./ max(cfg.radar_cost_normalizer, eps), 0), 1));

    vals = gpuArray.zeros(1, numel(visible), 'single');
    vals(visible) = snr_cost(visible);
    layer_gpu = max(layer_gpu, reshape(vals, gpu_ctx.H, gpu_ctx.W));
end
if gather_output
    layer = gather(layer_gpu);
else
    layer = layer_gpu;
end
end

function cost = analytical_sigmoid_radar_cost(radar, pts, z_asl, z_min_asl, cfg)
d = vecnorm(pts - radar, 2, 1);
d = max(d, 1);
sigmoid = 1 ./ (1 + exp(-cfg.sigmoid_k .* (z_asl - z_min_asl)));
if isfield(cfg, 'radar_effective_range_m') && isfinite(cfg.radar_effective_range_m) && cfg.radar_effective_range_m > 0
    softness = max(cfg.radar_range_softness_m, 1);
    range_gate = 1 ./ (1 + exp((d - cfg.radar_effective_range_m) ./ softness));
else
    range_gate = 1;
end
raw = range_gate .* sigmoid .* cfg.radar_cost_gain ./ (d .^ 4);
cost = single(min(max(raw ./ max(cfg.radar_cost_normalizer, eps), 0), 1));
end

function cost = radar_snr_cost(radar, pts, cfg)
d = vecnorm(pts - radar, 2, 1);
d = max(d, 1);
snr = cfg.radar_snr_constant ./ (d.^4);
snr_db = 10 * log10(max(snr, realmin));
cost = single(min(max((snr_db + 50) ./ 100, 0), 1));
end

function cost = radar_snr_cost_gpu(radar, pts_gpu, cfg)
radar_gpu = gpuArray(double(radar(:)));
dx = pts_gpu(1, :) - radar_gpu(1);
dy = pts_gpu(2, :) - radar_gpu(2);
dz = pts_gpu(3, :) - radar_gpu(3);
d2 = max(dx.^2 + dy.^2 + dz.^2, 1);
snr = cfg.radar_snr_constant ./ (d2 .^ 2);
snr_db = 10 * log10(max(snr, realmin));
cost = single(min(max((snr_db + 50) ./ 100, 0), 1));
end

function [waypoints, starts, goals, rcs_costs, planner_times, trajectory_dts, durations, radar_metrics, target_weights, visible_counts, visibility_budget_used, feasibility_stage, success_count] = generate_ait_trajectories(map, radars, sandwich, cost_field, z_layers_agl, cfg, map_seed)
waypoint_dim = 3 + double(cfg.export_time_channel);
waypoints = nan(cfg.K, cfg.N_samples, waypoint_dim, 'single');
starts = nan(cfg.K, 3, 'single');
goals = nan(cfg.K, 3, 'single');
rcs_costs = nan(cfg.K, 1, 'single');
planner_times = nan(cfg.K, 1, 'single');
trajectory_dts = nan(cfg.K, 1, 'single');
durations = nan(cfg.K, 1, 'single');
radar_metrics = empty_radar_metrics(cfg.K);
target_weights = nan(cfg.K, 1, 'single');
visible_counts = zeros(cfg.K, 1, 'int32');
visibility_budget_used = zeros(cfg.K, 1, 'int32');
feasibility_stage = zeros(cfg.K, 1, 'int32');

env = make_ait_env(map, cost_field, z_layers_agl, cfg);
params = make_ait_params(map, cfg);
use_par_traj = cfg.use_parallel && cfg.use_parallel_trajectories && ...
    ~strcmpi(string(cfg.planner_backend), "ait_star_mex") && ...
    can_use_parfor() && ~isempty(gcp('nocreate'));
if use_par_traj
    pool = gcp('nocreate');
    fprintf('  trajectories: parfor attempt batches on %d workers, batch_size=%d\n', pool.NumWorkers, cfg.trajectory_batch_size);
else
    fprintf('  trajectories: serial attempt batches, batch_size=%d\n', cfg.trajectory_batch_size);
end
success_count = 0;
attempt_cursor = 0;
use_scenario_quota = strcmpi(string(cfg.scenario_sampling_mode), "radar_relevant") && cfg.enforce_scenario_quota;
scenario_targets = scenario_targets_for_k(cfg.K, cfg.scenario_mix);
scenario_counts = zeros(3, 1);
while success_count < cfg.K && attempt_cursor < cfg.max_path_attempts
    remaining_successes = cfg.K - success_count;
    target_attempts_this_batch = max(remaining_successes * cfg.trajectory_attempt_oversample_factor, 1);
    batch_n = min([cfg.trajectory_batch_size, target_attempts_this_batch, cfg.max_path_attempts - attempt_cursor]);
    attempt_ids = attempt_cursor + (1:batch_n);
    attempt_cursor = attempt_cursor + batch_n;
    forced_scenario_ids = zeros(1, batch_n, 'int32');
    if use_scenario_quota
        planned_counts = scenario_counts;
        for j = 1:batch_n
            forced_scenario_ids(j) = next_needed_scenario_id(planned_counts, scenario_targets);
            if forced_scenario_ids(j) >= 1 && forced_scenario_ids(j) <= 3
                planned_counts(forced_scenario_ids(j)) = planned_counts(forced_scenario_ids(j)) + 1;
            end
        end
    end
    batch_tic = tic;
    fprintf('  attempts %d-%d started | success=%d/%d\n', attempt_ids(1), attempt_ids(end), success_count, cfg.K);

    if use_par_traj
        batch_results = cell(1, batch_n);
        parfor j = 1:batch_n
            batch_results{j} = run_single_attempt(map, radars, sandwich, cost_field, z_layers_agl, env, params, cfg, map_seed, attempt_ids(j), forced_scenario_ids(j));
        end
    else
        batch_results = cell(1, batch_n);
        for j = 1:batch_n
            batch_results{j} = run_single_attempt(map, radars, sandwich, cost_field, z_layers_agl, env, params, cfg, map_seed, attempt_ids(j), forced_scenario_ids(j));
        end
    end

    accepted_before = success_count;
    reject_reasons = strings(1, batch_n);
    for j = 1:batch_n
        result = batch_results{j};
        if ~result.success
            reject_reasons(j) = string(result.reject_reason);
            continue;
        end
        if use_scenario_quota && result.scenario_id >= 1 && result.scenario_id <= 3
            if scenario_counts(result.scenario_id) >= scenario_targets(result.scenario_id)
                reject_reasons(j) = "quota_full";
                continue;
            end
        end
        success_count = success_count + 1;
        waypoints(success_count, :, :) = single(reshape(result.waypoints, [1, cfg.N_samples, waypoint_dim]));
        starts(success_count, :) = single(result.start);
        goals(success_count, :) = single(result.goal);
        rcs_costs(success_count) = single(result.mean_cost);
        planner_times(success_count) = single(result.planner_time);
        trajectory_dts(success_count) = single(result.dt);
        durations(success_count) = single(result.duration);
        radar_metrics.scenario_ids(success_count) = int32(result.scenario_id);
        radar_metrics.straight_rcs_costs(success_count) = single(result.straight_mean_cost);
        radar_metrics.straight_visible_ratios(success_count) = single(result.straight_visible_ratio);
        radar_metrics.planned_visible_ratios(success_count) = single(result.planned_visible_ratio);
        radar_metrics.radar_improvement_ratios(success_count) = single(result.radar_improvement_ratio);
        target_weights(success_count) = single(result.target_weight);
        visible_counts(success_count) = int32(result.visible_count);
        visibility_budget_used(success_count) = int32(result.visibility_budget_used);
        feasibility_stage(success_count) = int32(result.feasibility_stage);
        if use_scenario_quota && result.scenario_id >= 1 && result.scenario_id <= 3
            scenario_counts(result.scenario_id) = scenario_counts(result.scenario_id) + 1;
        end
        fprintf('  traj %03d/%03d | attempt=%d | scen=%s | stage=%d | vis=%d/%d | w=%.2f | plan=%.2fs | traj_dt=%.3fs | T=%.1fs | mean_rcs=%.3f | straight_vis=%.2f planned_vis=%.2f improve=%.2f\n', ...
            success_count, cfg.K, result.attempt_id, scenario_name(result.scenario_id), ...
            result.feasibility_stage, result.visible_count, result.visibility_budget_used, result.target_weight, ...
            result.planner_time, result.dt, result.duration, ...
            result.mean_cost, result.straight_visible_ratio, result.planned_visible_ratio, result.radar_improvement_ratio);
        if success_count >= cfg.K
            break;
        end
    end
    fprintf('  attempts %d-%d done | accepted=%d | success=%d/%d | batch_time=%.1fs\n', ...
        attempt_ids(1), attempt_ids(end), success_count - accepted_before, success_count, cfg.K, toc(batch_tic));
    if cfg.log_rejection_summary
        print_rejection_summary(batch_results, reject_reasons, cfg);
    end
end

if success_count < cfg.K
    fprintf('  max attempts reached: generated %d/%d valid trajectories after %d attempts\n', ...
        success_count, cfg.K, attempt_cursor);
end
end

function print_rejection_summary(batch_results, reject_reasons, cfg)
reject_reasons = reject_reasons(strlength(reject_reasons) > 0);
if isempty(reject_reasons)
    return;
end
[uniq, ~, ic] = unique(reject_reasons, 'stable');
counts = accumarray(ic(:), 1);
parts = strings(1, numel(uniq));
for i = 1:numel(uniq)
    parts(i) = sprintf('%s=%d', uniq(i), counts(i));
end
fprintf('    rejects: %s\n', strjoin(parts, ', '));

n_examples = max(0, round(cfg.log_rejection_examples));
if n_examples <= 0
    return;
end
printed = 0;
for j = 1:numel(batch_results)
    result = batch_results{j};
    if result.success || strlength(string(result.reject_reason)) == 0
        continue;
    end
    fprintf('      ex attempt=%d reason=%s scen=%s plan=%.2fs straight_vis=%.2f planned_vis=%.2f improve=%.2f err=%s\n', ...
        result.attempt_id, string(result.reject_reason), scenario_name(result.scenario_id), ...
        result.planner_time, result.straight_visible_ratio, result.planned_visible_ratio, ...
        result.radar_improvement_ratio, string(result.error));
    printed = printed + 1;
    if printed >= n_examples
        break;
    end
end
end

function metrics = empty_radar_metrics(K)
metrics = struct();
metrics.scenario_ids = zeros(K, 1, 'int32');
metrics.straight_rcs_costs = nan(K, 1, 'single');
metrics.straight_visible_ratios = nan(K, 1, 'single');
metrics.planned_visible_ratios = nan(K, 1, 'single');
metrics.radar_improvement_ratios = nan(K, 1, 'single');
end

function result = run_single_attempt(map, radars, sandwich, cost_field, z_layers_agl, env, params, cfg, map_seed, attempt_id, forced_scenario_id)
result = empty_attempt_result(attempt_id);
rng(map_seed + attempt_id, 'twister');
if nargin < 11
    forced_scenario_id = int32(0);
end

[start_agl, goal_agl, scenario] = sample_start_goal_agl(map, radars, cost_field, z_layers_agl, cfg, forced_scenario_id);
result.scenario_id = scenario.id;
result.straight_mean_cost = scenario.straight_mean_cost;
result.straight_visible_ratio = scenario.straight_visible_ratio;
if isempty(start_agl)
    result.reject_reason = 'start_goal_sampling';
    return;
end
if ~isfinite(scenario.straight_mean_cost) || ~isfinite(scenario.straight_visible_ratio)
    straight = straight_line_agl(map, start_agl, goal_agl, cfg);
    straight_costs = query_cost_field_agl(map, cost_field, z_layers_agl, ...
        straight(:, 1), straight(:, 2), straight(:, 3));
    straight_stats = radar_path_stats(straight_costs, cfg);
    scenario.straight_mean_cost = single(straight_stats.mean_cost);
    scenario.straight_visible_ratio = single(straight_stats.visible_ratio);
    result.straight_mean_cost = scenario.straight_mean_cost;
    result.straight_visible_ratio = scenario.straight_visible_ratio;
end
% AIT* collision checks operate in absolute altitude (ASL): z_ASL is
% compared with terrain(x,y) to enforce AGL clearance. Dataset export stays
% terrain-relative, so accepted paths are converted back to z_AGL below.
start_asl = agl_to_asl(map, start_agl);
goal_asl = agl_to_asl(map, goal_agl);
params_attempt = apply_scenario_planner_budget(params, scenario.id, cfg);
region_mode = use_region_label_mode(cfg);
[path_asl, info, planner_time, planner_error] = run_planner_with_constraint_fallback(start_asl, goal_asl, env, params_attempt, cfg);
result.planner_time = single(planner_time);

if ~isempty(planner_error)
    result.error = planner_error;
    result.reject_reason = 'planner_error';
    return;
end
if isempty(path_asl) || ~isfield(info, 'success') || ~info.success
    result.reject_reason = 'planner_failed';
    if isfield(info, 'stop_reason')
        result.error = char(string(info.stop_reason));
    end
    return;
end

try
    traj_dt = choose_trajectory_dt(path_asl, cfg);
    if cfg.generate_trajectory
        wp_xyz = refine_or_smooth_path_agl(map, radars, sandwich, path_asl, start_agl, goal_agl, cfg, traj_dt);
    else
        wp_xyz = raw_path_resample_agl(map, path_asl, start_agl, goal_agl, cfg);
    end
    wp = append_time_channel(wp_xyz, traj_dt, cfg);
catch ME
    result.error = ME.message;
    result.reject_reason = 'smoothing_failed';
    return;
end

z_agl = wp(:, 3);
if any(~isfinite(wp(:))) || min(z_agl) < cfg.min_clearance
    result.reject_reason = 'invalid_waypoints';
    return;
end
costs = query_cost_field_agl(map, cost_field, z_layers_agl, wp(:, 1), wp(:, 2), wp(:, 3));
if any(~isfinite(costs))
    result.reject_reason = 'nonfinite_cost';
    return;
end
planned_stats = radar_path_stats(costs, cfg);
visibility_stats = radar_visibility_stats_agl(map, sandwich, wp(:, 1:3), cfg);
[feasibility_stage, visibility_budget_used, target_weight] = classify_trajectory_feasibility(visibility_stats.visible_count, cfg);
if ~region_mode && ~accept_planned_radar_stats(scenario, planned_stats, cfg)
    result.reject_reason = radar_acceptance_reject_reason(scenario, planned_stats, cfg);
    result.planned_visible_ratio = planned_stats.visible_ratio;
    result.mean_cost = single(planned_stats.mean_cost);
    result.radar_improvement_ratio = max(0, scenario.straight_mean_cost - planned_stats.mean_cost) / max(scenario.straight_mean_cost, eps);
    result.visible_count = int32(visibility_stats.visible_count);
    result.visibility_budget_used = int32(visibility_budget_used);
    result.feasibility_stage = int32(feasibility_stage);
    result.target_weight = single(target_weight);
    return;
end
if ~region_mode && ~accept_visibility_stats_agl(map, sandwich, wp(:, 1:3), scenario.id, cfg)
    result.reject_reason = 'visibility_acceptance';
    result.planned_visible_ratio = planned_stats.visible_ratio;
    result.mean_cost = single(planned_stats.mean_cost);
    result.radar_improvement_ratio = max(0, scenario.straight_mean_cost - planned_stats.mean_cost) / max(scenario.straight_mean_cost, eps);
    result.visible_count = int32(visibility_stats.visible_count);
    result.visibility_budget_used = int32(visibility_budget_used);
    result.feasibility_stage = int32(feasibility_stage);
    result.target_weight = single(target_weight);
    return;
end

result.success = true;
result.waypoints = single(wp);
result.start = single(start_agl(:)');
result.goal = single(goal_agl(:)');
result.mean_cost = single(mean(costs));
result.planner_time = single(planner_time);
result.dt = single(traj_dt);
result.duration = single(traj_dt * (cfg.N_samples - 1));
result.scenario_id = scenario.id;
result.straight_mean_cost = scenario.straight_mean_cost;
result.straight_visible_ratio = scenario.straight_visible_ratio;
result.planned_visible_ratio = planned_stats.visible_ratio;
result.radar_improvement_ratio = max(0, scenario.straight_mean_cost - planned_stats.mean_cost) / max(scenario.straight_mean_cost, eps);
result.visible_count = int32(visibility_stats.visible_count);
result.visibility_budget_used = int32(visibility_budget_used);
result.feasibility_stage = int32(feasibility_stage);
result.target_weight = single(target_weight);
end

function result = empty_attempt_result(attempt_id)
result = struct('success', false, 'attempt_id', attempt_id, 'waypoints', [], ...
    'start', [], 'goal', [], 'mean_cost', NaN, 'planner_time', NaN, ...
    'dt', NaN, 'duration', NaN, 'error', '', 'reject_reason', '', 'scenario_id', int32(0), ...
    'straight_mean_cost', NaN, 'straight_visible_ratio', NaN, ...
    'planned_visible_ratio', NaN, 'radar_improvement_ratio', NaN, ...
    'target_weight', single(1), 'visible_count', int32(0), ...
    'visibility_budget_used', int32(0), 'feasibility_stage', int32(0));
end

function [path_asl, info, planner_time, planner_error] = run_planner_with_constraint_fallback(start_asl, goal_asl, env, params, cfg)
planner_error = '';
planner_time = 0;
path_asl = [];
info = struct('success', false, 'stop_reason', 'not_started');

if strcmpi(string(cfg.planner_backend), "ait_star_mex")
    [path_asl, info, planner_time, planner_error] = run_mex_planner(start_asl, goal_asl, params, cfg);
    return;
end

if strcmpi(string(cfg.radar_constraint_mode), "hard_with_budget_fallback")
    hard_params = params;
    hard_params.stealth_weight = inf;
    hard_params.strict_zero_risk_tol = cfg.visible_cost_threshold;
    t_plan = tic;
    try
        [path_asl, ~, info] = ait_star_planner_direct(start_asl(:), goal_asl(:), env, hard_params);
    catch ME
        planner_error = ME.message;
        planner_time = toc(t_plan);
        return;
    end
    planner_time = toc(t_plan);
    if ~isempty(path_asl) && isfield(info, 'success') && info.success
        info.constraint_stage = 'hard';
        return;
    end

    soft_params = params;
    t_plan = tic;
    try
        [path_asl, ~, soft_info] = ait_star_planner_direct(start_asl(:), goal_asl(:), env, soft_params);
    catch ME
        planner_error = ME.message;
        planner_time = planner_time + toc(t_plan);
        return;
    end
    planner_time = planner_time + toc(t_plan);
    info = soft_info;
    info.constraint_stage = 'budget_fallback';
    return;
end

t_plan = tic;
try
    [path_asl, ~, info] = ait_star_planner_direct(start_asl(:), goal_asl(:), env, params);
catch ME
    planner_error = ME.message;
end
planner_time = toc(t_plan);
end

function [path_asl, info, planner_time, planner_error] = run_mex_planner(start_asl, goal_asl, params, cfg)
planner_error = '';
path_asl = [];
info = struct('success', false, 'stop_reason', 'not_started');
if strcmpi(string(cfg.radar_constraint_mode), "hard_with_budget_fallback")
    params.risk_limit = cfg.visible_cost_threshold;
else
    params.risk_limit = 1.0;
end
t_plan = tic;
try
    res = ait_star_mex('plan', double(start_asl(:)), double(goal_asl(:)), params);
    path_asl = double(res.path);
    info = normalize_mex_plan_info(res);
catch ME
    planner_error = ME.message;
end
planner_time = toc(t_plan);
end

function info = normalize_mex_plan_info(res)
info = struct();
if isfield(res, 'success')
    info.success = logical(res.success);
else
    info.success = ~isempty(res.path);
end
info.stop_reason = get_struct_field_or(res, 'stop_reason', '');
info.cost = get_struct_field_or(res, 'cost', NaN);
info.batches_run = get_struct_field_or(res, 'batches_run', NaN);
info.nodes_in_tree = get_struct_field_or(res, 'nodes_in_tree', NaN);
info.total_nodes = get_struct_field_or(res, 'total_nodes', NaN);
info.first_solution_batch = get_struct_field_or(res, 'first_solution_batch', NaN);
info.first_solution_time_s = get_struct_field_or(res, 'first_solution_time_s', NaN);
info.edge_checks = get_struct_field_or(res, 'edge_checks', NaN);
info.constraint_stage = 'mex';
end

function val = get_struct_field_or(s, name, default_val)
if isstruct(s) && isfield(s, name)
    val = s.(name);
else
    val = default_val;
end
end

function tf = use_region_label_mode(cfg)
tf = strcmpi(string(cfg.label_mode), "trajectory_and_region") && ...
     strcmpi(string(cfg.region_label_type), "2d_corridor");
end

function [stage, budget_used, target_weight] = classify_trajectory_feasibility(visible_count, cfg)
schedule = unique(max(0, round(double(cfg.visibility_budget_schedule(:)'))));
if isempty(schedule)
    schedule = [0, 2, 5, 10, 20];
end
visible_count = max(0, round(double(visible_count)));
idx = find(schedule >= visible_count, 1, 'first');
if isempty(idx)
    stage = int32(2);
    budget_used = int32(-1);
    target_weight = single(cfg.min_violation_target_weight);
elseif idx == 1
    stage = int32(0);
    budget_used = int32(schedule(idx));
    target_weight = single(1.0);
else
    stage = int32(1);
    budget_used = int32(schedule(idx));
    target_weight = single(cfg.relaxed_target_weight);
end
end

function traj_dt = choose_trajectory_dt(path_asl, cfg)
if ~strcmpi(string(cfg.dt_mode), "adaptive")
    traj_dt = cfg.dt;
    return;
end
path_len = polyline_length(path_asl);
duration = cfg.dt_safety_factor * path_len / max(cfg.nominal_speed_mps, eps);
traj_dt = duration / max(cfg.N_samples - 1, 1);
traj_dt = clamp(traj_dt, cfg.min_dt, cfg.max_dt);
end

function len = polyline_length(P)
if isempty(P) || size(P, 2) < 2
    len = 0;
    return;
end
len = sum(vecnorm(diff(P, 1, 2), 2, 1));
end

function wp = append_time_channel(wp_xyz, traj_dt, cfg)
if ~cfg.export_time_channel
    wp = wp_xyz;
    return;
end
t = (0:cfg.N_samples - 1).' * traj_dt;
wp = [wp_xyz, t];
end

function env = make_ait_env(map, cost_field, z_layers_agl, cfg)
env = struct();
env.terrain = struct();
env.terrain.Z = double(map.elevation);
env.terrain.N_vec = map.xs(:)';
env.terrain.E_vec = map.ys(:)';
env.terrain.x_vec = map.xs(:)';
env.terrain.y_vec = map.ys(:)';
env.terrain.dx = double(map.resolution_m);
env.terrain.dy = double(map.resolution_m);
env.terrain.alt_vec = double(z_layers_agl(:)');
env.risk_grid = double(cost_field);
env.N_vec = map.xs(:)';
env.E_vec = map.ys(:)';
env.alt_vec = double(z_layers_agl(:)');
env.risk_altitude_mode = 'agl';
env.risk_extrapolation = 'nearest';
env.visibility_threshold = 1.0;
env.bounds_world = [map.xs(1), map.xs(end), map.ys(1), map.ys(end), ...
    min(map.elevation(:)) + cfg.min_clearance, max(map.elevation(:)) + cfg.max_clearance];
end

function params = make_ait_params(map, cfg)
diag_len = hypot(range(map.xs), range(map.ys));
params = struct();
params.max_batches = cfg.planner_max_batches;
params.batch_size = cfg.planner_batch_size;
params.max_nodes = cfg.planner_max_nodes;
params.eps_global = 0.12;
params.stealth_weight = cfg.planner_weight;
params.gamma = cfg.planner_altitude_weight;
params.preferred_agl = cfg.planner_preferred_agl;
params.strict_zero_risk_tol = 1e-6;
params.edge_check_samples = cfg.edge_check_samples;
params.post_solution_batches = cfg.planner_post_solution_batches;
params.goal_radius = max(40, 0.018 * diag_len);
params.min_connection_radius = max(35, 0.012 * diag_len);
params.max_connection_radius = max(120, 0.08 * diag_len);
params.eta_radius = 1.8;
params.max_time = cfg.planner_max_time;
params.min_clearance = cfg.min_clearance;
params.max_clearance = cfg.max_clearance;
params.reverse_k_neighbors = cfg.reverse_k_neighbors;
params.forward_k_neighbors = cfg.forward_k_neighbors;
params.goal_sample_fraction = cfg.goal_sample_fraction;
params.line_sample_fraction = cfg.line_sample_fraction;
params.goal_probe_count = cfg.goal_probe_count;
params.profile_ait = false;
params.neighbor_precompute_mode = 'batch';
use_nested_trajectory_parallel = cfg.use_parallel && cfg.use_parallel_trajectories;
params.use_parallel_reverse = cfg.use_parallel && cfg.use_parallel_ait_reverse && ~use_nested_trajectory_parallel;
params.use_parallel_edge_checks = cfg.use_parallel && cfg.use_parallel_ait_edge_checks && ~use_nested_trajectory_parallel;
params.parallel_threshold_neighbors = 64;
params.parallel_threshold_edges = 32;
params.reverse_batch_repair_limit = cfg.reverse_batch_repair_limit;
params.forward_expand_chunk = cfg.forward_expand_chunk;
params.kd_rebuild_interval = 1;
params.graph_rebuild_interval = cfg.graph_rebuild_interval;
params.reverse_full_rebuild_interval = cfg.reverse_full_rebuild_interval;
params.reverse_debug = false;
end

function params = apply_scenario_planner_budget(params, scenario_id, cfg)
idx = max(1, min(3, double(scenario_id)));
time_scale = scenario_scale_value(cfg.scenario_planner_time_scale, idx, 1.0);
batch_scale = scenario_scale_value(cfg.scenario_planner_batch_scale, idx, 1.0);
node_scale = scenario_scale_value(cfg.scenario_planner_node_scale, idx, 1.0);
post_scale = scenario_scale_value(cfg.scenario_planner_post_solution_scale, idx, 1.0);

params.max_time = max(params.max_time, cfg.planner_max_time * time_scale);
params.max_batches = max(params.max_batches, ceil(cfg.planner_max_batches * batch_scale));
params.max_nodes = max(params.max_nodes, ceil(cfg.planner_max_nodes * node_scale));
params.post_solution_batches = max(params.post_solution_batches, ceil(cfg.planner_post_solution_batches * post_scale));
end

function val = scenario_scale_value(values, idx, fallback)
values = double(values(:)');
if isempty(values)
    val = fallback;
    return;
end
idx = max(1, min(idx, numel(values)));
val = values(idx);
if ~isfinite(val) || val <= 0
    val = fallback;
end
end

function [start_agl, goal_agl, scenario] = sample_start_goal_agl(map, radars, cost_field, z_layers_agl, cfg, forced_scenario_id)
start_agl = [];
goal_agl = [];
scenario = empty_scenario();
if nargin < 6
    forced_scenario_id = int32(0);
end
if strcmpi(string(cfg.scenario_sampling_mode), "radar_relevant")
    [start_agl, goal_agl, scenario] = sample_radar_relevant_start_goal_agl(map, radars, cost_field, z_layers_agl, cfg, forced_scenario_id);
    return;
end

diag_len = hypot(range(map.xs), range(map.ys));
min_dist = cfg.min_start_goal_dist_frac * diag_len;
max_dist = min(cfg.max_start_goal_dist_m, 0.92 * diag_len);
if cfg.enforce_start_goal_reachability
    min_dist = min(min_dist, 0.65 * max_dist);
else
    max_dist = 0.98 * diag_len;
end

mode = lower(string(cfg.start_goal_mode));
if mode == "diagonal" || mode == "opposite_corners" || mode == "opposite_edges" || mode == "different_quadrants"
    for attempt = 1:120
        [s, g] = sample_separated_endpoints_agl(map, cfg, mode);
        if isempty(s) || isempty(g)
            continue;
        end
        d = norm(s(1:2) - g(1:2));
        if d >= min_dist && (~cfg.enforce_start_goal_reachability || d <= max_dist)
            start_agl = s;
            goal_agl = g;
            scenario = build_scenario(1, NaN, NaN);
            return;
        end
    end
    return;
end

for attempt = 1:120
    s = sample_endpoint_agl(map, cfg);
    g = sample_endpoint_agl(map, cfg);
    if isempty(s) || isempty(g)
        continue;
    end
    d = norm(s(1:2) - g(1:2));
    if d >= min_dist && d <= max_dist
        start_agl = s;
        goal_agl = g;
        scenario = build_scenario(1, NaN, NaN);
        return;
    end
end
end

function [start_agl, goal_agl, scenario] = sample_radar_relevant_start_goal_agl(map, radars, cost_field, z_layers_agl, cfg, forced_scenario_id)
start_agl = [];
goal_agl = [];
scenario = empty_scenario();
if nargin < 6
    forced_scenario_id = int32(0);
end
for attempt = 1:cfg.radar_relevant_sampling_attempts
    if forced_scenario_id >= 1 && forced_scenario_id <= 3
        scenario_id = int32(forced_scenario_id);
    else
        scenario_id = sample_scenario_id(cfg);
    end
    mode = lower(string(cfg.start_goal_mode));
    [s, g] = sample_separated_endpoints_agl(map, cfg, mode);
    if isempty(s) || isempty(g)
        continue;
    end
    if ~start_goal_distance_ok(map, s, g, cfg)
        continue;
    end
    if ~line_is_relevant_to_radars(s, g, radars, map, scenario_id)
        continue;
    end
    straight = straight_line_agl(map, s, g, cfg);
    costs = query_cost_field_agl(map, cost_field, z_layers_agl, straight(:, 1), straight(:, 2), straight(:, 3));
    stats = radar_path_stats(costs, cfg);
    if straight_stats_match_scenario(stats, scenario_id, cfg)
        start_agl = s;
        goal_agl = g;
        scenario = build_scenario(scenario_id, stats.mean_cost, stats.visible_ratio);
        return;
    end
end
end

function ok = start_goal_distance_ok(map, s, g, cfg)
diag_len = hypot(range(map.xs), range(map.ys));
min_dist = cfg.min_start_goal_dist_frac * diag_len;
max_dist = min(cfg.max_start_goal_dist_m, 0.92 * diag_len);
if cfg.enforce_start_goal_reachability
    min_dist = min(min_dist, 0.65 * max_dist);
else
    max_dist = 0.98 * diag_len;
end
d = norm(s(1:2) - g(1:2));
ok = d >= min_dist && d <= max_dist;
end

function id = sample_scenario_id(cfg)
mix = double(cfg.scenario_mix(:)');
if numel(mix) ~= 3 || sum(mix) <= 0
    mix = [0.15, 0.50, 0.35];
end
mix = mix ./ sum(mix);
u = rand();
if u < mix(1)
    id = int32(1); % easy
elseif u < mix(1) + mix(2)
    id = int32(2); % medium
else
    id = int32(3); % hard
end
end

function targets = scenario_targets_for_k(K, mix)
mix = double(mix(:));
if numel(mix) ~= 3 || sum(mix) <= 0
    mix = [0.15; 0.50; 0.35];
end
mix = mix ./ sum(mix);
raw = mix .* double(K);
targets = floor(raw);
remaining = K - sum(targets);
[~, order] = sort(raw - targets, 'descend');
for i = 1:remaining
    targets(order(i)) = targets(order(i)) + 1;
end
targets = targets(:);
end

function id = next_needed_scenario_id(counts, targets)
deficit = targets(:) - counts(:);
if all(deficit <= 0)
    id = int32(0);
    return;
end
ratio = deficit ./ max(targets(:), 1);
ratio(deficit <= 0) = -Inf;
[~, idx] = max(ratio);
id = int32(idx);
end

function name = scenario_name(id)
if id == 1
    name = 'easy';
elseif id == 2
    name = 'medium';
elseif id == 3
    name = 'hard';
else
    name = 'unknown';
end
end

function scenario = empty_scenario()
scenario = build_scenario(int32(0), NaN, NaN);
end

function scenario = build_scenario(id, straight_mean_cost, straight_visible_ratio)
scenario = struct();
scenario.id = int32(id);
scenario.straight_mean_cost = single(straight_mean_cost);
scenario.straight_visible_ratio = single(straight_visible_ratio);
end

function ok = line_is_relevant_to_radars(start_agl, goal_agl, radars, map, scenario_id)
if isempty(radars)
    ok = true;
    return;
end
diag_len = hypot(range(map.xs), range(map.ys));
seg = goal_agl(1:2) - start_agl(1:2);
seg_len = norm(seg);
if seg_len <= 1e-6
    ok = false;
    return;
end
rel = radars(:, 1:2) - start_agl(1:2);
u = min(max((rel * seg(:)) ./ (seg_len^2), 0), 1);
closest = start_agl(1:2) + u .* seg;
dist_to_line = sqrt(sum((radars(:, 1:2) - closest).^2, 2));
min_dist = min(dist_to_line);
if scenario_id == 1
    ok = min_dist <= 0.65 * diag_len;
elseif scenario_id == 2
    ok = min_dist <= 0.35 * diag_len;
else
    ok = min_dist <= 0.20 * diag_len;
end
end

function wp = straight_line_agl(map, start_agl, goal_agl, cfg)
u = linspace(0, 1, cfg.N_samples).';
xy = (1 - u) .* start_agl(1:2) + u .* goal_agl(1:2);
terrain = terrain_height(map, xy(:, 1), xy(:, 2));
start_asl = terrain_height(map, start_agl(1), start_agl(2)) + start_agl(3);
goal_asl = terrain_height(map, goal_agl(1), goal_agl(2)) + goal_agl(3);
z_asl = (1 - u) .* start_asl + u .* goal_asl;
z_agl = z_asl - terrain(:);
z_agl = min(max(z_agl, cfg.min_clearance), cfg.max_clearance);
wp = [xy, z_agl(:)];
end

function wp_xyz = raw_path_resample_agl(map, path_asl, start_agl, goal_agl, cfg)
path_asl = remove_duplicate_points(path_asl);
if isempty(path_asl) || size(path_asl, 2) < 2
    path_asl = [agl_to_asl(map, start_agl), agl_to_asl(map, goal_agl)];
end
path_asl(:, 1) = agl_to_asl(map, start_agl);
path_asl(:, end) = agl_to_asl(map, goal_agl);
path_resampled = resample_path_arclength(path_asl, cfg.N_samples);
path_resampled(1, :) = clamp(path_resampled(1, :), map.xs(1), map.xs(end));
path_resampled(2, :) = clamp(path_resampled(2, :), map.ys(1), map.ys(end));
terrain = terrain_height(map, path_resampled(1, :), path_resampled(2, :));
z_agl = path_resampled(3, :) - terrain(:)';
z_agl = min(max(z_agl, cfg.min_clearance), cfg.max_clearance);
wp_xyz = [path_resampled(1, :).', path_resampled(2, :).', z_agl(:)];
end

function stats = radar_path_stats(costs, cfg)
costs = double(costs(:));
stats = struct();
stats.mean_cost = mean(costs, 'omitnan');
stats.visible_ratio = mean(costs > cfg.visible_cost_threshold, 'omitnan');
end

function ok = straight_stats_match_scenario(stats, scenario_id, cfg)
if scenario_id == 1
    range = cfg.easy_straight_visible_range;
elseif scenario_id == 2
    range = cfg.medium_straight_visible_range;
else
    range = cfg.hard_straight_visible_range;
end
ok = stats.visible_ratio >= range(1) && stats.visible_ratio <= range(2);
end

function ok = accept_planned_radar_stats(scenario, planned_stats, cfg)
if scenario.id == 1
    ok = planned_stats.visible_ratio <= max(cfg.max_planned_visible_ratio, scenario.straight_visible_ratio + 0.10);
    return;
end
improvement = max(0, scenario.straight_mean_cost - planned_stats.mean_cost) / max(scenario.straight_mean_cost, eps);
ok = planned_stats.visible_ratio <= cfg.max_planned_visible_ratio && ...
     planned_stats.mean_cost <= cfg.max_planned_mean_radar_cost && ...
     improvement >= cfg.min_radar_improvement_ratio;
end

function reason = radar_acceptance_reject_reason(scenario, planned_stats, cfg)
if scenario.id == 1
    visible_limit = max(cfg.max_planned_visible_ratio, scenario.straight_visible_ratio + 0.10);
    if planned_stats.visible_ratio > visible_limit
        reason = 'radar_acceptance_visible';
    else
        reason = 'radar_acceptance_unknown';
    end
    return;
end

improvement = max(0, scenario.straight_mean_cost - planned_stats.mean_cost) / max(scenario.straight_mean_cost, eps);
if planned_stats.visible_ratio > cfg.max_planned_visible_ratio
    reason = 'radar_acceptance_visible';
elseif planned_stats.mean_cost > cfg.max_planned_mean_radar_cost
    reason = 'radar_acceptance_mean_cost';
elseif improvement < cfg.min_radar_improvement_ratio
    reason = 'radar_acceptance_improvement';
else
    reason = 'radar_acceptance_unknown';
end
end

function ok = accept_visibility_stats_agl(map, sandwich, wp_agl, scenario_id, cfg)
if ~cfg.enforce_visibility_acceptance
    ok = true;
    return;
end
stats = radar_visibility_stats_agl(map, sandwich, wp_agl, cfg);
max_violations = max_visibility_violations_for_scenario(scenario_id, cfg);
ok = stats.visible_count <= max_violations;
end

function stats = radar_visibility_stats_agl(map, sandwich, wp_agl, cfg)
wp_agl = double(wp_agl);
N = size(wp_agl, 1);
stats = struct('visible_count', 0, 'visible_ratio', 0, 'max_margin_m', NaN);
if N == 0 || ~isfield(sandwich, 'min_visible_alt') || isempty(sandwich.min_visible_alt)
    return;
end

terrain = terrain_height(map, wp_agl(:, 1), wp_agl(:, 2));
z_asl = terrain(:) + wp_agl(:, 3);
min_visible_asl = interp_regular_2d(double(sandwich.min_visible_alt), ...
    double(map.xs), double(map.ys), wp_agl(:, 1), wp_agl(:, 2));
min_visible_asl = min_visible_asl(:);

margin = z_asl(:) - min_visible_asl(:) + double(cfg.visibility_margin_m);
finite_mask = isfinite(min_visible_asl) & isfinite(margin);
visible = finite_mask & (margin > 0);
stats.visible_count = sum(visible(:));
stats.visible_ratio = stats.visible_count / max(N, 1);
if any(finite_mask(:))
    stats.max_margin_m = max(margin(finite_mask));
end
end

function max_violations = max_visibility_violations_for_scenario(scenario_id, cfg)
limits = double(cfg.max_planned_visibility_violations(:)');
if isempty(limits)
    limits = [12, 8, 5];
end
if numel(limits) == 1
    max_violations = limits(1);
    return;
end
idx = max(1, min(3, double(scenario_id)));
idx = min(idx, numel(limits));
max_violations = limits(idx);
end

function [s, g] = sample_separated_endpoints_agl(map, cfg, mode)
s = [];
g = [];
x_margin = min(cfg.start_goal_margin_m, 0.2 * range(map.xs));
y_margin = min(cfg.start_goal_margin_m, 0.2 * range(map.ys));
x_min = map.xs(1) + x_margin;
x_max = map.xs(end) - x_margin;
y_min = map.ys(1) + y_margin;
y_max = map.ys(end) - y_margin;
if x_min >= x_max || y_min >= y_max
    x_min = map.xs(1); x_max = map.xs(end);
    y_min = map.ys(1); y_max = map.ys(end);
end

if mode == "diagonal" || mode == "opposite_corners"
    pairs = [x_min, y_min, x_max, y_max; ...
             x_min, y_max, x_max, y_min; ...
             x_max, y_max, x_min, y_min; ...
             x_max, y_min, x_min, y_max];
elseif mode == "different_quadrants"
    xm = 0.5 * (x_min + x_max);
    ym = 0.5 * (y_min + y_max);
    quadrants = [x_min, xm, y_min, ym; ...
                 xm, x_max, y_min, ym; ...
                 x_min, xm, ym, y_max; ...
                 xm, x_max, ym, y_max];
    q1 = randi(4);
    q2 = randi(3);
    other = setdiff(1:4, q1);
    q2 = other(q2);
    s = sample_endpoint_in_box_agl(map, cfg, quadrants(q1, :));
    g = sample_endpoint_in_box_agl(map, cfg, quadrants(q2, :));
    return;
else
    pairs = [x_min, y_min, x_max, y_min; ...
             x_min, y_max, x_max, y_max; ...
             x_min, y_min, x_min, y_max; ...
             x_max, y_min, x_max, y_max; ...
             x_min, 0.5 * (y_min + y_max), x_max, 0.5 * (y_min + y_max); ...
             0.5 * (x_min + x_max), y_min, 0.5 * (x_min + x_max), y_max];
end

pair = pairs(randi(size(pairs, 1)), :);
jx = cfg.start_goal_corner_jitter_frac * max(x_max - x_min, eps);
jy = cfg.start_goal_corner_jitter_frac * max(y_max - y_min, eps);
sx = clamp(pair(1) + (2 * rand() - 1) * jx, x_min, x_max);
sy = clamp(pair(2) + (2 * rand() - 1) * jy, y_min, y_max);
gx = clamp(pair(3) + (2 * rand() - 1) * jx, x_min, x_max);
gy = clamp(pair(4) + (2 * rand() - 1) * jy, y_min, y_max);
s = endpoint_at_xy_agl(map, sx, sy, cfg);
g = endpoint_at_xy_agl(map, gx, gy, cfg);
end

function p = sample_endpoint_in_box_agl(map, cfg, box)
p = [];
for attempt = 1:40
    x = box(1) + rand() * max(box(2) - box(1), eps);
    y = box(3) + rand() * max(box(4) - box(3), eps);
    p = endpoint_at_xy_agl(map, x, y, cfg);
    if ~isempty(p)
        return;
    end
end
end

function p = sample_endpoint_agl(map, cfg)
p = [];
x_margin = min(cfg.start_goal_margin_m, 0.2 * range(map.xs));
y_margin = min(cfg.start_goal_margin_m, 0.2 * range(map.ys));
x_min = map.xs(1) + x_margin;
x_max = map.xs(end) - x_margin;
y_min = map.ys(1) + y_margin;
y_max = map.ys(end) - y_margin;
if x_min >= x_max || y_min >= y_max
    x_min = map.xs(1); x_max = map.xs(end);
    y_min = map.ys(1); y_max = map.ys(end);
end
for attempt = 1:40
    x = x_min + rand() * (x_max - x_min);
    y = y_min + rand() * (y_max - y_min);
    h = terrain_height(map, x, y);
    if isfinite(h)
        z_agl = rand_range(cfg.start_goal_agl_range);
        p = [x, y, z_agl];
        return;
    end
end
end

function p = endpoint_at_xy_agl(map, x, y, cfg)
p = [];
h = terrain_height(map, x, y);
if isfinite(h)
    z_agl = rand_range(cfg.start_goal_agl_range);
    p = [x, y, z_agl];
end
end

function p_asl = agl_to_asl(map, p_agl)
p_asl = [p_agl(1); p_agl(2); terrain_height(map, p_agl(1), p_agl(2)) + p_agl(3)];
end

function wp = smooth_resample_yaw_agl(map, path_asl, cfg)
path_asl = remove_duplicate_points(path_asl);
arc_path = resample_path_arclength(path_asl, cfg.N_samples);
ctrl = resample_path_arclength(path_asl, cfg.smooth_control_points);
seg_len = vecnorm(diff(ctrl, 1, 2), 2, 1);
seg_len = max(seg_len, 1);
smooth_params = sample_smoothing_params(cfg);
[time_points, velBC, accBC] = allocate_optimized_time_points(ctrl, seg_len, smooth_params, cfg);
try
    warn_near = warning('off', 'MATLAB:nearlySingularMatrix');
    warn_sing = warning('off', 'MATLAB:singularMatrix');
    cleanup_warn = onCleanup(@() restore_warnings(warn_near, warn_sing));
    [pos, ~, ~, ~, ~, ~, ~, ~] = minsnappolytraj(ctrl, time_points, cfg.N_samples, ...
        'VelocityBoundaryCondition', velBC, ...
        'AccelerationBoundaryCondition', accBC);
catch
    pos = arc_path;
end

pos(1, :) = clamp(pos(1, :), map.xs(1), map.xs(end));
pos(2, :) = clamp(pos(2, :), map.ys(1), map.ys(end));
agl = pos(3, :) - terrain_height(map, pos(1, :), pos(2, :));
agl = min(max(agl, cfg.min_clearance), cfg.max_clearance);
wp = zeros(cfg.N_samples, 4);
wp(:, 1:3) = [pos(1, :).', pos(2, :).', agl(:)];
dx = [diff(wp(:, 1)); wp(end, 1) - wp(end-1, 1)];
dy = [diff(wp(:, 2)); wp(end, 2) - wp(end-1, 2)];
yaw = atan2(dy, dx);
wp(:, 4) = wrap_pi(yaw);
end

function wp = refine_or_smooth_path_agl(map, radars, sandwich, path_asl, start_agl, goal_agl, cfg, traj_dt)
fallback = smooth_resample_yaw_agl(map, path_asl, cfg);
fallback_xyz = fallback(:, 1:3);
if ~(cfg.use_local_optimizer && cfg.casadi_available)
    wp = fallback_xyz;
    return;
end
try
    if cfg.use_local_optimizer_template
        wp = refine_path_casadi_template_agl(map, radars, sandwich, path_asl, start_agl, goal_agl, cfg, traj_dt);
    else
        wp = refine_path_casadi_agl(map, radars, sandwich, path_asl, start_agl, goal_agl, cfg, traj_dt);
    end
catch ME
    if cfg.smooth_verbose
        warning('generate_radar_uav_dataset:LocalOptimizerFallback', ...
            'Local optimizer failed; using smoothed AIT* path (%s).', ME.message);
    end
    wp = fallback_xyz;
end
end

function wp_agl = refine_path_casadi_agl(map, radars, sandwich, path_asl, start_agl, goal_agl, cfg, traj_dt)
import casadi.*
N = cfg.N_samples;
dt = traj_dt;
warm_asl = resample_path_arclength(path_asl, N);
warm_asl(:, 1) = agl_to_asl(map, start_agl);
warm_asl(:, end) = agl_to_asl(map, goal_agl);

terrain_warm = terrain_height(map, warm_asl(1, :), warm_asl(2, :));
zmin_by_radar = sample_zmin_asl_for_path(map, sandwich, warm_asl);
opti = Opti();
P = opti.variable(3, N);
V = opti.variable(3, N);
A = opti.variable(3, N - 1);

opti.subject_to(P(:, 1) == warm_asl(:, 1));
opti.subject_to(P(:, end) == warm_asl(:, end));
for k = 1:(N - 1)
    opti.subject_to(P(:, k + 1) == P(:, k) + dt * V(:, k) + 0.5 * dt^2 * A(:, k));
    opti.subject_to(V(:, k + 1) == V(:, k) + dt * A(:, k));
    opti.subject_to(sumsqr(V(:, k)) <= cfg.v_max^2);
    opti.subject_to(sumsqr(A(:, k)) <= cfg.a_max^2);
end
opti.subject_to(sumsqr(V(:, N)) <= cfg.v_max^2);
opti.subject_to(P(1, :) >= map.xs(1));
opti.subject_to(P(1, :) <= map.xs(end));
opti.subject_to(P(2, :) >= map.ys(1));
opti.subject_to(P(2, :) <= map.ys(end));
opti.subject_to(P(3, :) >= terrain_warm + cfg.min_clearance);
opti.subject_to(P(3, :) <= terrain_warm + cfg.max_clearance);

V0 = finite_difference_velocity(warm_asl, dt);
A0 = finite_difference_acceleration(V0, dt);
opti.set_initial(P, warm_asl);
opti.set_initial(V, V0);
opti.set_initial(A, A0);

J = cfg.local_tracking_weight * sumsqr(P - warm_asl) + cfg.local_control_weight * sumsqr(A);
if N > 2
    J = J + cfg.local_smooth_weight * sumsqr(A(:, 2:end) - A(:, 1:end-1));
end
for k = 1:N
    for r = 1:size(radars, 1)
        dr = P(:, k) - radars(r, :).';
        range2 = sumsqr(dr) + 1.0;
        sigmoid = 1 / (1 + exp(-cfg.sigmoid_k * (P(3, k) - zmin_by_radar(k, r))));
        J = J + cfg.local_radar_weight * cfg.radar_cost_gain * sigmoid / (range2^2);
    end
end
opti.minimize(J);

solver_order = local_solver_order(cfg.local_solver_resolved);
last_error = '';
for i = 1:numel(solver_order)
    solver_name = solver_order{i};
    try
        solver_opts = casadi_solver_options(solver_name, cfg);
        if cfg.smooth_verbose
            opti.solver(solver_name, solver_opts);
            sol = opti.solve();
        else
            evalc('opti.solver(solver_name, solver_opts); sol = opti.solve();');
        end
        P_sol = full(sol.value(P));
        P_sol(1, :) = clamp(P_sol(1, :), map.xs(1), map.xs(end));
        P_sol(2, :) = clamp(P_sol(2, :), map.ys(1), map.ys(end));
        agl = P_sol(3, :) - terrain_height(map, P_sol(1, :), P_sol(2, :));
        wp_agl = [P_sol(1, :).', P_sol(2, :).', agl(:)];
if any(~isfinite(wp_agl(:))) || min(wp_agl(:, 3)) < cfg.min_clearance
    error('generate_radar_uav_dataset:InvalidLocalSolution', ...
        'Local optimizer produced invalid terrain clearance.');
end
wp_agl(:, 3) = min(wp_agl(:, 3), cfg.max_clearance);
return;
    catch ME
        last_error = ME.message;
    end
end
error('generate_radar_uav_dataset:LocalOptimizerFailed', '%s', last_error);
end

function wp_agl = refine_path_casadi_template_agl(map, radars, sandwich, path_asl, start_agl, goal_agl, cfg, traj_dt)
N = cfg.N_samples;
Rmax = cfg.max_radars;
warm_asl = resample_path_arclength(path_asl, N);
warm_asl(:, 1) = agl_to_asl(map, start_agl);
warm_asl(:, end) = agl_to_asl(map, goal_agl);
terrain_warm = terrain_height(map, warm_asl(1, :), warm_asl(2, :));
zmin_actual = sample_zmin_asl_for_path(map, sandwich, warm_asl);

radar_pad = zeros(3, Rmax);
zmin_pad = repmat(warm_asl(3, :).', 1, Rmax);
radar_mask = zeros(1, Rmax);
R = min(size(radars, 1), Rmax);
if R > 0
    radar_pad(:, 1:R) = radars(1:R, :).';
    if ~isempty(zmin_actual)
        zmin_pad(:, 1:R) = zmin_actual(:, 1:R);
    end
    radar_mask(1:R) = 1;
end

template = get_local_optimizer_template(cfg, N, Rmax);
opti = template.opti;
opti.set_value(template.warm, warm_asl);
opti.set_value(template.terrain, terrain_warm(:).');
opti.set_value(template.zmin, zmin_pad);
opti.set_value(template.radars, radar_pad);
opti.set_value(template.radar_mask, radar_mask);
opti.set_value(template.start_p, warm_asl(:, 1));
opti.set_value(template.goal_p, warm_asl(:, end));
opti.set_value(template.x_min_p, map.xs(1));
opti.set_value(template.x_max_p, map.xs(end));
opti.set_value(template.y_min_p, map.ys(1));
opti.set_value(template.y_max_p, map.ys(end));
opti.set_value(template.v_max_p, cfg.v_max);
opti.set_value(template.a_max_p, cfg.a_max);
opti.set_value(template.dt_p, traj_dt);
opti.set_initial(template.P, warm_asl);
opti.set_initial(template.V, finite_difference_velocity(warm_asl, traj_dt));
opti.set_initial(template.A, finite_difference_acceleration(finite_difference_velocity(warm_asl, traj_dt), traj_dt));

if cfg.smooth_verbose
    sol = opti.solve();
else
    evalc('sol = opti.solve();');
end
P_sol = full(sol.value(template.P));
P_sol(1, :) = clamp(P_sol(1, :), map.xs(1), map.xs(end));
P_sol(2, :) = clamp(P_sol(2, :), map.ys(1), map.ys(end));
agl = P_sol(3, :) - terrain_height(map, P_sol(1, :), P_sol(2, :));
wp_agl = [P_sol(1, :).', P_sol(2, :).', agl(:)];
if any(~isfinite(wp_agl(:))) || min(wp_agl(:, 3)) < cfg.min_clearance
    error('generate_radar_uav_dataset:InvalidLocalSolution', ...
        'Template local optimizer produced invalid terrain clearance.');
end
wp_agl(:, 3) = min(wp_agl(:, 3), cfg.max_clearance);
end

function template = get_local_optimizer_template(cfg, N, Rmax)
persistent cache
solver_name = cfg.local_solver_resolved;
key = sprintf('N%d_R%d_%s_v%d_a%d_dt%.6g_k%.6g', N, Rmax, solver_name, ...
    round(cfg.v_max), round(cfg.a_max), cfg.dt, cfg.sigmoid_k);
if ~isempty(cache) && isfield(cache, 'key') && strcmp(cache.key, key)
    template = cache.template;
    return;
end

import casadi.*
opti = Opti();
P = opti.variable(3, N);
V = opti.variable(3, N);
A = opti.variable(3, N - 1);
warm = opti.parameter(3, N);
terrain = opti.parameter(1, N);
zmin = opti.parameter(N, Rmax);
radars_p = opti.parameter(3, Rmax);
radar_mask = opti.parameter(1, Rmax);
start_p = opti.parameter(3, 1);
goal_p = opti.parameter(3, 1);
x_min_p = opti.parameter();
x_max_p = opti.parameter();
y_min_p = opti.parameter();
y_max_p = opti.parameter();
v_max_p = opti.parameter();
a_max_p = opti.parameter();
dt_p = opti.parameter();

opti.subject_to(P(:, 1) == start_p);
opti.subject_to(P(:, end) == goal_p);
for k = 1:(N - 1)
    opti.subject_to(P(:, k + 1) == P(:, k) + dt_p * V(:, k) + 0.5 * dt_p^2 * A(:, k));
    opti.subject_to(V(:, k + 1) == V(:, k) + dt_p * A(:, k));
    opti.subject_to(sumsqr(V(:, k)) <= v_max_p^2);
    opti.subject_to(sumsqr(A(:, k)) <= a_max_p^2);
end
opti.subject_to(sumsqr(V(:, N)) <= v_max_p^2);
opti.subject_to(P(1, :) >= x_min_p);
opti.subject_to(P(1, :) <= x_max_p);
opti.subject_to(P(2, :) >= y_min_p);
opti.subject_to(P(2, :) <= y_max_p);
opti.subject_to(P(3, :) >= terrain + cfg.min_clearance);
opti.subject_to(P(3, :) <= terrain + cfg.max_clearance);

J = cfg.local_tracking_weight * sumsqr(P - warm) + cfg.local_control_weight * sumsqr(A);
if N > 2
    J = J + cfg.local_smooth_weight * sumsqr(A(:, 2:end) - A(:, 1:end-1));
end
for k = 1:N
    for r = 1:Rmax
        dr = P(:, k) - radars_p(:, r);
        range2 = sumsqr(dr) + 1.0;
        sigmoid = 1 / (1 + exp(-cfg.sigmoid_k * (P(3, k) - zmin(k, r))));
        J = J + radar_mask(r) * cfg.local_radar_weight * cfg.radar_cost_gain * sigmoid / (range2^2);
    end
end
opti.minimize(J);
opti.solver(solver_name, casadi_solver_options(solver_name, cfg));

template = struct();
template.opti = opti;
template.P = P;
template.V = V;
template.A = A;
template.warm = warm;
template.terrain = terrain;
template.zmin = zmin;
template.radars = radars_p;
template.radar_mask = radar_mask;
template.start_p = start_p;
template.goal_p = goal_p;
template.x_min_p = x_min_p;
template.x_max_p = x_max_p;
template.y_min_p = y_min_p;
template.y_max_p = y_max_p;
template.v_max_p = v_max_p;
template.a_max_p = a_max_p;
template.dt_p = dt_p;
cache = struct('key', key, 'template', template);
end

function order = local_solver_order(primary_solver)
primary_solver = lower(char(primary_solver));
if strcmp(primary_solver, 'fatrop')
    order = {'fatrop', 'ipopt'};
elseif strcmp(primary_solver, 'ipopt')
    order = {'ipopt', 'fatrop'};
else
    order = {'ipopt'};
end
end

function opts = casadi_solver_options(solver_name, cfg)
opts = struct();
solver_name = lower(char(solver_name));
if strcmp(solver_name, 'ipopt')
    opts.ipopt.print_level = 0;
    opts.ipopt.max_iter = cfg.local_optimizer_max_iter;
    opts.ipopt.sb = 'yes';
    opts.print_time = 0;
elseif strcmp(solver_name, 'fatrop')
    opts.print_time = 0;
    opts.fatrop.max_iter = cfg.local_optimizer_max_iter;
else
    opts.print_time = 0;
end
end

function zmin_asl = sample_zmin_asl_for_path(map, sandwich, path_asl)
N = size(path_asl, 2);
if ~isfield(sandwich, 'z_min_by_radar') || isempty(sandwich.z_min_by_radar)
    zmin_asl = repmat(path_asl(3, :).', 1, 0);
    return;
end
R = size(sandwich.z_min_by_radar, 3);
terrain = terrain_height(map, path_asl(1, :), path_asl(2, :));
zmin_asl = zeros(N, R);
for r = 1:R
    zmin_agl = interp_regular_2d(double(sandwich.z_min_by_radar(:, :, r)), ...
        double(map.xs), double(map.ys), path_asl(1, :), path_asl(2, :));
    zmin_asl(:, r) = terrain(:) + zmin_agl(:);
end
end

function V = finite_difference_velocity(P, dt)
V = zeros(size(P));
if size(P, 2) == 1
    return;
end
V(:, 1:end-1) = diff(P, 1, 2) ./ dt;
V(:, end) = V(:, end-1);
end

function A = finite_difference_acceleration(V, dt)
A = diff(V, 1, 2) ./ dt;
if isempty(A)
    A = zeros(size(V, 1), 0);
end
end

function smooth_params = sample_smoothing_params(cfg)
smooth_params = struct();
smooth_params.v_max = rand_range(cfg.smooth_v_max_range);
smooth_params.a_max = rand_range(cfg.smooth_a_max_range);
smooth_params.cruise_speed = rand_range(cfg.smooth_cruise_speed_range);
smooth_params.aggressiveness = rand_range(cfg.smooth_aggressiveness_range);
smooth_params.v_safety = rand_range(cfg.smooth_v_safety_range);
smooth_params.samples_per_seg = cfg.smooth_samples_per_seg;
smooth_params.vel_bc_mode = cfg.smooth_vel_bc_mode;
end

function [time_points, velBC, accBC] = allocate_optimized_time_points(ctrl, seg_len, smooth_params, cfg)
n_wp = size(ctrl, 2);
velBC = build_velocity_boundary_conditions(ctrl, smooth_params);
accBC = nan(3, n_wp);
accBC(:, 1) = [0; 0; 0];
accBC(:, end) = [0; 0; 0];

fallback_times = max(seg_len ./ max(smooth_params.cruise_speed, 1), 0.5);
segment_times = fallback_times;
if cfg.use_optimized_time_allocation && exist('optimize_time_allocation', 'file') == 2
    opt_params = struct();
    opt_params.aggressiveness = smooth_params.aggressiveness;
    opt_params.samples_per_seg = smooth_params.samples_per_seg;
    opt_params.v_safety = smooth_params.v_safety;
    try
        warn_near = warning('off', 'MATLAB:nearlySingularMatrix');
        warn_sing = warning('off', 'MATLAB:singularMatrix');
        cleanup_warn = onCleanup(@() restore_warnings(warn_near, warn_sing));
        if cfg.smooth_verbose
            [segment_times, ~] = optimize_time_allocation(ctrl, smooth_params.v_max, ...
                smooth_params.a_max, velBC, opt_params);
        else
            evalc('[segment_times, ~] = optimize_time_allocation(ctrl, smooth_params.v_max, smooth_params.a_max, velBC, opt_params);');
        end
    catch
        segment_times = fallback_times;
    end
end
segment_times = max(segment_times(:).', 0.1);
time_points = [0, cumsum(segment_times)];
end

function velBC = build_velocity_boundary_conditions(ctrl, smooth_params)
n_wp = size(ctrl, 2);
velBC = nan(3, n_wp);
velBC(:, 1) = [0; 0; 0];
velBC(:, end) = [0; 0; 0];
if strcmpi(smooth_params.vel_bc_mode, 'bisector') && n_wp > 2
    cruise_speed = smooth_params.cruise_speed;
    for i = 2:n_wp-1
        dir_in = ctrl(:, i) - ctrl(:, i-1);
        dir_out = ctrl(:, i+1) - ctrl(:, i);
        len_in = norm(dir_in);
        len_out = norm(dir_out);
        if len_in > 1e-6
            dir_in = dir_in / len_in;
        else
            dir_in = [0; 0; 0];
        end
        if len_out > 1e-6
            dir_out = dir_out / len_out;
        else
            dir_out = [0; 0; 0];
        end
        dir_bisect = dir_in + dir_out;
        len_bisect = norm(dir_bisect);
        if len_bisect > 1e-6
            dir_bisect = dir_bisect / len_bisect;
        else
            dir_bisect = dir_in;
        end
        velBC(:, i) = cruise_speed * dir_bisect;
    end
end
end

function restore_warnings(varargin)
for i = 1:nargin
    warning(varargin{i});
end
end

function P = remove_duplicate_points(P)
if size(P, 2) <= 1
    return;
end
d = vecnorm(diff(P, 1, 2), 2, 1);
P = P(:, [true, d > 1e-6]);
if size(P, 2) < 2
    P = [P, P + [1e-3; 0; 0]];
end
end

function Pq = resample_path_arclength(P, n)
P = remove_duplicate_points(P);
d = [0, cumsum(vecnorm(diff(P, 1, 2), 2, 1))];
if d(end) <= 1e-9
    Pq = repmat(P(:, 1), 1, n);
    return;
end
[d_unique, ia] = unique(d, 'stable');
P_unique = P(:, ia);
s = linspace(0, d_unique(end), n);
Pq = interp1(d_unique(:), P_unique.', s(:), 'pchip').';
end

function c = query_cost_field_agl(map, cost_field, z_layers_agl, x, y, z_agl)
c = interp_regular_3d(double(cost_field), double(map.xs), double(map.ys), double(z_layers_agl(:)'), x, y, z_agl);
c = max(0, min(1, c(:)));
end

function h = terrain_height(map, x, y)
h = interp_regular_2d(double(map.elevation), double(map.xs), double(map.ys), x, y);
end

function v = interp_regular_2d(V, xs, ys, x, y)
out_size = size(x);
x = double(x);
y = double(y);
H = numel(ys);
W = numel(xs);
dx = max(xs(2) - xs(1), eps);
dy = max(ys(2) - ys(1), eps);
fx = 1 + (x(:).' - xs(1)) ./ dx;
fy = 1 + (y(:).' - ys(1)) ./ dy;
fx = min(max(fx, 1), W);
fy = min(max(fy, 1), H);
x0 = min(max(floor(fx), 1), W - 1);
y0 = min(max(floor(fy), 1), H - 1);
ax = fx - x0;
ay = fy - y0;
x1 = x0 + 1;
y1 = y0 + 1;

i00 = sub2ind([H, W], y0, x0);
i10 = sub2ind([H, W], y0, x1);
i01 = sub2ind([H, W], y1, x0);
i11 = sub2ind([H, W], y1, x1);
v = (1 - ax) .* (1 - ay) .* V(i00) + ...
    ax .* (1 - ay) .* V(i10) + ...
    (1 - ax) .* ay .* V(i01) + ...
    ax .* ay .* V(i11);
v = reshape(v, out_size);
end

function v = interp_regular_3d(V, xs, ys, zs, x, y, z)
H = numel(ys);
W = numel(xs);
ZL = numel(zs);
dx = max(xs(2) - xs(1), eps);
dy = max(ys(2) - ys(1), eps);
dz = max(zs(2) - zs(1), eps);
fx = 1 + (double(x(:).') - xs(1)) ./ dx;
fy = 1 + (double(y(:).') - ys(1)) ./ dy;
fz = 1 + (double(z(:).') - zs(1)) ./ dz;
fx = min(max(fx, 1), W);
fy = min(max(fy, 1), H);
fz = min(max(fz, 1), ZL);
x0 = min(max(floor(fx), 1), W - 1);
y0 = min(max(floor(fy), 1), H - 1);
z0 = min(max(floor(fz), 1), ZL - 1);
ax = fx - x0;
ay = fy - y0;
az = fz - z0;
x1 = x0 + 1;
y1 = y0 + 1;
z1 = z0 + 1;

i000 = sub2ind([H, W, ZL], y0, x0, z0);
i100 = sub2ind([H, W, ZL], y0, x1, z0);
i010 = sub2ind([H, W, ZL], y1, x0, z0);
i110 = sub2ind([H, W, ZL], y1, x1, z0);
i001 = sub2ind([H, W, ZL], y0, x0, z1);
i101 = sub2ind([H, W, ZL], y0, x1, z1);
i011 = sub2ind([H, W, ZL], y1, x0, z1);
i111 = sub2ind([H, W, ZL], y1, x1, z1);

c00 = (1 - ax) .* V(i000) + ax .* V(i100);
c10 = (1 - ax) .* V(i010) + ax .* V(i110);
c01 = (1 - ax) .* V(i001) + ax .* V(i101);
c11 = (1 - ax) .* V(i011) + ax .* V(i111);
c0 = (1 - ay) .* c00 + ay .* c10;
c1 = (1 - ay) .* c01 + ay .* c11;
v = (1 - az) .* c0 + az .* c1;
v = v(:);
end

function mesh = build_numeric_terrain_mesh(map)
xs = double(map.xs(:)');
ys = double(map.ys(:)');
Z = double(map.elevation);
H = numel(ys);
W = numel(xs);
[X, Y] = meshgrid(xs, ys);
vertices = [X(:), Y(:), Z(:)];

grid_idx = reshape(uint32(1:(H * W)), H, W);
a = grid_idx(1:end-1, 1:end-1);
b = grid_idx(1:end-1, 2:end);
c = grid_idx(2:end, 2:end);
d = grid_idx(2:end, 1:end-1);
triangles = [a(:), b(:), c(:); a(:), c(:), d(:)];
tri_rows = double(triangles);
tri_x = reshape(vertices(tri_rows, 1), size(tri_rows));
tri_y = reshape(vertices(tri_rows, 2), size(tri_rows));
tri_z = reshape(vertices(tri_rows, 3), size(tri_rows));

mesh = struct();
mesh.vertices = vertices;
mesh.triangles = triangles;
mesh.v0 = vertices(tri_rows(:, 1), :).';
mesh.v1 = vertices(tri_rows(:, 2), :).';
mesh.v2 = vertices(tri_rows(:, 3), :).';
mesh.tri_x_min = min(tri_x, [], 2).';
mesh.tri_x_max = max(tri_x, [], 2).';
mesh.tri_y_min = min(tri_y, [], 2).';
mesh.tri_y_max = max(tri_y, [], 2).';
mesh.tri_z_min = min(tri_z, [], 2).';
mesh.tri_z_max = max(tri_z, [], 2).';
mesh.bounds = [xs(1), xs(end), ys(1), ys(end), min(Z(:)), max(Z(:))];
mesh.cell_h = H - 1;
mesh.cell_w = W - 1;
mesh.num_cells = mesh.cell_h * mesh.cell_w;
mesh.dx = xs(2) - xs(1);
mesh.dy = ys(2) - ys(1);
mesh.x_min = xs(1);
mesh.x_max = xs(end);
mesh.y_min = ys(1);
mesh.y_max = ys(end);
end

function blocked = ray_tracing(mesh, radar, targets, chunk_size, pair_batch_size)
%RAY_TRACING Exact terrain ray tracing using vectorized Moller-Trumbore.
n = size(targets, 2);
blocked = false(1, n);
if n == 0
    return;
end

chunk_size = max(1, round(chunk_size));
pair_batch_size = max(1, round(pair_batch_size));
eps_t = 1e-7;
for first = 1:chunk_size:n
    last = min(n, first + chunk_size - 1);
    cols = first:last;
    target_chunk = targets(:, cols);
    dir_chunk = target_chunk - radar;

    blocked_chunk = false(1, numel(cols));
    for local_idx = 1:numel(cols)
        tri_ids = ray_candidate_triangles_dda(mesh, radar(1), radar(2), ...
            target_chunk(1, local_idx), target_chunk(2, local_idx));
        if isempty(tri_ids)
            continue;
        end
        for pair_first = 1:pair_batch_size:numel(tri_ids)
            pair_last = min(numel(tri_ids), pair_first + pair_batch_size - 1);
            tri_batch = tri_ids(pair_first:pair_last);
            ray_origin = radar(:, ones(1, numel(tri_batch)));
            ray_dir = dir_chunk(:, local_idx * ones(1, numel(tri_batch)));
            [hit, t, ~, ~, ~] = moller_trumbore(ray_origin, ray_dir, ...
                mesh.v0(:, tri_batch), mesh.v1(:, tri_batch), mesh.v2(:, tri_batch));

            % Direction is the full radar-to-target vector, so t is the
            % segment fraction. Ignore hits at radar/target endpoints.
            if any(hit & t > eps_t & t < (1 - eps_t))
                blocked_chunk(local_idx) = true;
                break;
            end
        end
    end
    blocked(cols) = blocked_chunk;
end
end

function tri_ids = ray_candidate_triangles_dda(mesh, x0, y0, x1, y1)
%RAY_CANDIDATE_TRIANGLES_DDA DEM cells crossed by the ray's xy footprint.
x0 = min(max(x0, mesh.x_min), mesh.x_max);
x1 = min(max(x1, mesh.x_min), mesh.x_max);
y0 = min(max(y0, mesh.y_min), mesh.y_max);
y1 = min(max(y1, mesh.y_min), mesh.y_max);

ix = point_to_cell_index(x0, mesh.x_min, mesh.dx, mesh.cell_w);
iy = point_to_cell_index(y0, mesh.y_min, mesh.dy, mesh.cell_h);
ix_end = point_to_cell_index(x1, mesh.x_min, mesh.dx, mesh.cell_w);
iy_end = point_to_cell_index(y1, mesh.y_min, mesh.dy, mesh.cell_h);

dx_ray = x1 - x0;
dy_ray = y1 - y0;
if dx_ray > 0
    step_x = 1;
    next_x = mesh.x_min + ix * mesh.dx;
    t_max_x = (next_x - x0) / dx_ray;
    t_delta_x = mesh.dx / dx_ray;
elseif dx_ray < 0
    step_x = -1;
    next_x = mesh.x_min + (ix - 1) * mesh.dx;
    t_max_x = (next_x - x0) / dx_ray;
    t_delta_x = -mesh.dx / dx_ray;
else
    step_x = 0;
    t_max_x = inf;
    t_delta_x = inf;
end

if dy_ray > 0
    step_y = 1;
    next_y = mesh.y_min + iy * mesh.dy;
    t_max_y = (next_y - y0) / dy_ray;
    t_delta_y = mesh.dy / dy_ray;
elseif dy_ray < 0
    step_y = -1;
    next_y = mesh.y_min + (iy - 1) * mesh.dy;
    t_max_y = (next_y - y0) / dy_ray;
    t_delta_y = -mesh.dy / dy_ray;
else
    step_y = 0;
    t_max_y = inf;
    t_delta_y = inf;
end

max_steps = mesh.cell_h + mesh.cell_w + 4;
cells = zeros(1, max_steps, 'uint32');
count = 0;
for step_count = 1:max_steps
    count = count + 1;
    cells(count) = uint32(sub2ind([mesh.cell_h, mesh.cell_w], iy, ix));
    if ix == ix_end && iy == iy_end
        break;
    end
    if t_max_x < t_max_y
        ix = min(max(ix + step_x, 1), mesh.cell_w);
        t_max_x = t_max_x + t_delta_x;
    elseif t_max_y < t_max_x
        iy = min(max(iy + step_y, 1), mesh.cell_h);
        t_max_y = t_max_y + t_delta_y;
    else
        ix = min(max(ix + step_x, 1), mesh.cell_w);
        iy = min(max(iy + step_y, 1), mesh.cell_h);
        t_max_x = t_max_x + t_delta_x;
        t_max_y = t_max_y + t_delta_y;
    end
end

cells = unique(cells(1:count));
tri_ids = double([cells(:).', mesh.num_cells + cells(:).']);
end

function idx = point_to_cell_index(x, x_min, dx, max_idx)
idx = floor((x - x_min) / dx) + 1;
idx = min(max(idx, 1), max_idx);
end

function a = wrap_pi(a)
a = mod(a + pi, 2*pi) - pi;
a(a <= -pi) = a(a <= -pi) + 2*pi;
end

function validate_bundle(bundle, cfg)
actual_k = size(bundle.waypoints, 1);
waypoint_dim = 3 + double(cfg.export_time_channel);
assert(isa(bundle.elevation, 'single') && isequal(size(bundle.elevation), cfg.dem_size));
assert(isa(bundle.min_visible_alt, 'single') && isequal(size(bundle.min_visible_alt), cfg.dem_size));
assert(isa(bundle.thickness, 'single') && isequal(size(bundle.thickness), cfg.dem_size));
assert(isa(bundle.scene_input, 'single') && isequal(size(bundle.scene_input), [cfg.dem_size, 3]));
assert(isa(bundle.resolution_m, 'single') && isscalar(bundle.resolution_m));
assert(isa(bundle.origin_xy, 'single') && numel(bundle.origin_xy) == 2);
assert(isa(bundle.radar_positions, 'single') && size(bundle.radar_positions, 2) == 3);
assert(isa(bundle.cost_field, 'single') && isequal(size(bundle.cost_field), [cfg.dem_size, numel(cfg.z_layers_agl)]));
assert(isa(bundle.z_layers, 'single') && numel(bundle.z_layers) == numel(cfg.z_layers_agl));
assert(actual_k >= 1 && actual_k <= cfg.K, 'Invalid trajectory count in bundle.');
assert(isa(bundle.waypoints, 'single') && isequal(size(bundle.waypoints), [actual_k, cfg.N_samples, waypoint_dim]));
assert(isa(bundle.starts, 'single') && isequal(size(bundle.starts), [actual_k, 3]));
assert(isa(bundle.goals, 'single') && isequal(size(bundle.goals), [actual_k, 3]));
assert(isa(bundle.rcs_costs, 'single') && numel(bundle.rcs_costs) == actual_k);
assert(isa(bundle.planner_times, 'single') && numel(bundle.planner_times) == actual_k);
assert(isa(bundle.trajectory_dts, 'single') && numel(bundle.trajectory_dts) == actual_k);
assert(isa(bundle.durations, 'single') && numel(bundle.durations) == actual_k);
assert(isa(bundle.scenario_ids, 'int32') && numel(bundle.scenario_ids) == actual_k);
assert(isa(bundle.straight_rcs_costs, 'single') && numel(bundle.straight_rcs_costs) == actual_k);
assert(isa(bundle.straight_visible_ratios, 'single') && numel(bundle.straight_visible_ratios) == actual_k);
assert(isa(bundle.planned_visible_ratios, 'single') && numel(bundle.planned_visible_ratios) == actual_k);
assert(isa(bundle.radar_improvement_ratios, 'single') && numel(bundle.radar_improvement_ratios) == actual_k);
assert(isa(bundle.target_weights, 'single') && numel(bundle.target_weights) == actual_k);
assert(isa(bundle.visible_counts, 'int32') && numel(bundle.visible_counts) == actual_k);
assert(isa(bundle.visibility_budget_used, 'int32') && numel(bundle.visibility_budget_used) == actual_k);
assert(isa(bundle.feasibility_stage, 'int32') && numel(bundle.feasibility_stage) == actual_k);
assert(isa(bundle.region_heuristic, 'single') && isequal(size(bundle.region_heuristic), cfg.dem_size));
assert(isa(bundle.region_safe_mask, 'single') && isequal(size(bundle.region_safe_mask), cfg.dem_size));
assert(isa(bundle.region_threat_distance, 'single') && isequal(size(bundle.region_threat_distance), cfg.dem_size));
assert(min(bundle.region_heuristic(:)) >= 0 && max(bundle.region_heuristic(:)) <= 1, 'Region heuristic out of [0,1].');
assert(min(bundle.region_safe_mask(:)) >= 0 && max(bundle.region_safe_mask(:)) <= 1, 'Region safe_mask out of [0,1].');
assert(min(bundle.region_threat_distance(:)) >= 0, 'Region threat_distance contains negative values.');

names = fieldnames(bundle);
for i = 1:numel(names)
    data = bundle.(names{i});
    if ischar(data) || isstring(data)
        continue;
    end
    assert(all(isfinite(double(data(:)))), 'Non-finite values in bundle.%s.', names{i});
end
z_agl = bundle.waypoints(:, :, 3);
assert(all(z_agl(:) >= 0), 'Waypoint z_AGL contains underground points.');
assert(min(z_agl(:)) >= cfg.validation_min_agl, 'Waypoint z_AGL violates minimum validation clearance.');
assert(max(z_agl(:)) <= cfg.max_clearance + 1e-3, 'Waypoint z_AGL exceeds maximum configured clearance.');
if cfg.export_time_channel
    t = bundle.waypoints(:, :, 4);
    assert(all(abs(t(:, 1)) <= 1e-6), 'Trajectory time channel must start at zero.');
    dt_check = diff(t, 1, 2);
    assert(all(dt_check(:) > 0), 'Trajectory time channel must be strictly increasing.');
end
assert(min(bundle.cost_field(:)) >= 0 && max(bundle.cost_field(:)) <= 1, 'Cost field out of [0,1].');
validate_python_export_contract(bundle, cfg);
end

function validate_python_export_contract(bundle, cfg)
% This generator exists only to feed the Python training pipeline. Keep this
% contract aligned with data.loader.MapBundle and data.dataset.TrajectoryDataset.
H = cfg.dem_size(1);
W = cfg.dem_size(2);
res = double(bundle.resolution_m);
expected_xs = ((0:W-1) - (W - 1) / 2) * res;
expected_ys = ((0:H-1) - (H - 1) / 2) * res;
origin = double(bundle.origin_xy(:)');
assert(numel(origin) == 2, 'origin_xy must be [east_min, north_min].');
assert(abs(origin(1) - expected_xs(1)) < 1e-4 && abs(origin(2) - expected_ys(1)) < 1e-4, ...
    'origin_xy must equal the first centered DEM grid coordinate.');

xy_limit = [max(abs(expected_xs)), max(abs(expected_ys))] + 5 * res;
all_xy = [reshape(bundle.waypoints(:, :, 1), [], 1), reshape(bundle.waypoints(:, :, 2), [], 1); ...
          double(bundle.starts(:, 1:2)); double(bundle.goals(:, 1:2)); ...
          double(bundle.radar_positions(:, 1:2))];
assert(all(abs(all_xy(:, 1)) <= xy_limit(1)) && all(abs(all_xy(:, 2)) <= xy_limit(2)), ...
    'Export XY must be map-centered ENU and remain near map bounds.');

assert(isequal(size(bundle.scene_input), [H, W, 3]), ...
    'scene_input must be MATLAB [H,W,3]; Python loader converts it to [3,H,W].');
assert(isequal(size(bundle.cost_field), [H, W, numel(cfg.z_layers_agl)]), ...
    'cost_field must be MATLAB [H,W,Z]; Python loader restores [H,W,Z].');
assert(isequal(size(bundle.region_heuristic), [H, W]) && ...
       isequal(size(bundle.region_safe_mask), [H, W]) && ...
       isequal(size(bundle.region_threat_distance), [H, W]), ...
    'region labels must be MATLAB [H,W].');
waypoint_dim = 3 + double(cfg.export_time_channel);
assert(size(bundle.waypoints, 3) == waypoint_dim, 'waypoints must be [K,N,3/4].');
z_agl = bundle.waypoints(:, :, 3);
assert(all(z_agl(:) >= 0), 'waypoint z must be AGL, not ASL.');
end

function write_bundle_h5(out_file, bundle, cfg)
tmp_file = [out_file, '.tmp'];
if exist(tmp_file, 'file')
    delete(tmp_file);
end
if exist(out_file, 'file')
    delete(out_file);
end

write_h5_single(tmp_file, '/dem/elevation', bundle.elevation);
write_h5_single(tmp_file, '/dem/min_visible_alt', bundle.min_visible_alt);
write_h5_single(tmp_file, '/dem/thickness', bundle.thickness);
write_h5_single(tmp_file, '/dem/resolution_m', bundle.resolution_m);
write_h5_single(tmp_file, '/dem/origin_xy', bundle.origin_xy);
write_h5_single(tmp_file, '/scene/input', bundle.scene_input);
write_h5_single(tmp_file, '/radars/positions', bundle.radar_positions);
write_h5_single(tmp_file, '/cost_field/data', bundle.cost_field);
write_h5_single(tmp_file, '/cost_field/z_layers', bundle.z_layers);
write_h5_single(tmp_file, '/cost_field/resolution_m', bundle.cost_resolution_m);
write_h5_single(tmp_file, '/region/heuristic', bundle.region_heuristic);
write_h5_single(tmp_file, '/region/safe_mask', bundle.region_safe_mask);
write_h5_single(tmp_file, '/region/threat_distance', bundle.region_threat_distance);
write_h5_single(tmp_file, '/trajectories/waypoints', bundle.waypoints);
write_h5_single(tmp_file, '/trajectories/starts', bundle.starts);
write_h5_single(tmp_file, '/trajectories/goals', bundle.goals);
write_h5_single(tmp_file, '/trajectories/rcs_costs', bundle.rcs_costs);
write_h5_single(tmp_file, '/trajectories/planner_times', bundle.planner_times);
write_h5_single(tmp_file, '/trajectories/dt', bundle.trajectory_dts);
write_h5_single(tmp_file, '/trajectories/duration', bundle.durations);
write_h5_int32(tmp_file, '/trajectories/scenario_ids', bundle.scenario_ids);
write_h5_single(tmp_file, '/trajectories/straight_rcs_costs', bundle.straight_rcs_costs);
write_h5_single(tmp_file, '/trajectories/straight_visible_ratios', bundle.straight_visible_ratios);
write_h5_single(tmp_file, '/trajectories/planned_visible_ratios', bundle.planned_visible_ratios);
write_h5_single(tmp_file, '/trajectories/radar_improvement_ratios', bundle.radar_improvement_ratios);
write_h5_single(tmp_file, '/trajectories/target_weights', bundle.target_weights);
write_h5_int32(tmp_file, '/trajectories/visible_counts', bundle.visible_counts);
write_h5_int32(tmp_file, '/trajectories/visibility_budget_used', bundle.visibility_budget_used);
write_h5_int32(tmp_file, '/trajectories/feasibility_stage', bundle.feasibility_stage);
if isfield(bundle, 'bspline_control_points')
    write_h5_single(tmp_file, '/trajectories/bspline_control_points', bundle.bspline_control_points);
end
h5writeatt(tmp_file, '/', 'coordinate_convention', 'centered ENU; z in trajectories is AGL');
h5writeatt(tmp_file, '/', 'planner', char(bundle.planner_label));
h5writeatt(tmp_file, '/', 'planner_backend', char(bundle.planner_backend));
h5writeatt(tmp_file, '/', 'generated_trajectory', double(bundle.generated_trajectory));
h5writeatt(tmp_file, '/', 'path_source', char(bundle.path_source));
h5writeatt(tmp_file, '/', 'dt', double(bundle.dt));
h5writeatt(tmp_file, '/', 'dt_mode', char(bundle.dt_mode));
h5writeatt(tmp_file, '/', 'waypoint_channels', char(bundle.waypoint_channels));
h5writeatt(tmp_file, '/', 'v_max', double(bundle.v_max));
h5writeatt(tmp_file, '/', 'a_max', double(bundle.a_max));
h5writeatt(tmp_file, '/', 'start_goal_mode', char(bundle.start_goal_mode));
h5writeatt(tmp_file, '/', 'enforce_start_goal_reachability', double(bundle.enforce_start_goal_reachability));
h5writeatt(tmp_file, '/', 'scenario_sampling_mode', char(cfg.scenario_sampling_mode));
h5writeatt(tmp_file, '/', 'label_mode', char(cfg.label_mode));
h5writeatt(tmp_file, '/', 'region_label_type', char(cfg.region_label_type));
h5writeatt(tmp_file, '/', 'region_corridor_mode', char(cfg.region_corridor_mode));
h5writeatt(tmp_file, '/', 'region_sigma_min_m', double(cfg.region_sigma_min_m));
h5writeatt(tmp_file, '/', 'region_sigma_max_m', double(cfg.region_sigma_max_m));
h5writeatt(tmp_file, '/', 'region_sigma_threat_gain', double(cfg.region_sigma_threat_gain));
h5writeatt(tmp_file, '/', 'radar_constraint_mode', char(cfg.radar_constraint_mode));
h5writeatt(tmp_file, '/', 'scenario_ids', '1=easy,2=medium,3=hard');
h5writeatt(tmp_file, '/', 'feasibility_stage_ids', '0=strict,1=relaxed_budget,2=min_violation');
h5writeatt(tmp_file, '/', 'visible_cost_threshold', double(cfg.visible_cost_threshold));
h5writeatt(tmp_file, '/', 'max_planned_visible_ratio', double(cfg.max_planned_visible_ratio));
h5writeatt(tmp_file, '/', 'min_radar_improvement_ratio', double(cfg.min_radar_improvement_ratio));
h5writeatt(tmp_file, '/', 'radar_effective_range_frac', double(cfg.radar_effective_range_frac));
h5writeatt(tmp_file, '/', 'radar_range_softness_frac', double(cfg.radar_range_softness_frac));
h5writeatt(tmp_file, '/', 'local_solver', char(bundle.solver_backend));
h5writeatt(tmp_file, '/', 'casadi_version', char(bundle.casadi_version));
h5writeatt(tmp_file, '/', 'used_local_optimizer', double(bundle.used_local_optimizer));
h5writeatt(tmp_file, '/', 'num_trajectories', double(bundle.num_trajectories));
h5writeatt(tmp_file, '/', 'requested_K', double(bundle.requested_K));
h5writeatt(tmp_file, '/', 'partial', double(bundle.partial));
movefile(tmp_file, out_file);
end

function write_h5_single(filename, dataset, data)
data = single(data);
if isscalar(data)
    write_h5_scalar_single(filename, dataset, data);
    return;
elseif isvector(data)
    dims = size(data);
else
    dims = size(data);
end
h5create(filename, dataset, dims, 'Datatype', 'single');
h5write(filename, dataset, data);
end

function write_h5_scalar_single(filename, dataset, data)
[parent_path, dataset_name] = split_h5_path(dataset);
file_id = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
group_id = H5G.open(file_id, parent_path);
space_id = H5S.create('H5S_SCALAR');
file_type = H5T.copy('H5T_IEEE_F32LE');
mem_type = H5T.copy('H5T_NATIVE_FLOAT');
dset_id = H5D.create(group_id, dataset_name, file_type, space_id, 'H5P_DEFAULT');
cleanup = onCleanup(@() close_h5_handles(dset_id, mem_type, file_type, space_id, group_id, file_id));
H5D.write(dset_id, mem_type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', single(data));
end

function write_h5_int32(filename, dataset, data)
data = int32(data);
dims = size(data);
h5create(filename, dataset, dims, 'Datatype', 'int32');
h5write(filename, dataset, data);
end

function [parent_path, dataset_name] = split_h5_path(dataset)
parts = strsplit(dataset, '/');
parts = parts(~cellfun('isempty', parts));
dataset_name = parts{end};
if isscalar(parts)
    parent_path = '/';
else
    parent_path = ['/', strjoin(parts(1:end-1), '/')];
end
end

function close_h5_handles(dset_id, mem_type, file_type, space_id, group_id, file_id)
H5D.close(dset_id);
H5T.close(mem_type);
H5T.close(file_type);
H5S.close(space_id);
H5G.close(group_id);
H5F.close(file_id);
end

function status = maybe_start_pool(cfg)
status = struct('enabled', false, 'workers', 0, 'message', '');
if ~cfg.use_parallel
    status.message = 'disabled by cfg.use_parallel=false';
    return;
end
if ~can_use_parfor()
    status.message = 'unavailable: Parallel Computing Toolbox/parfor not detected';
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        if strcmpi(cfg.parallel_pool, 'threads')
            parpool('threads');
        else
            parpool;
        end
        pool = gcp('nocreate');
    end
    if isempty(pool)
        status.message = 'unavailable: pool did not start';
    else
        status.enabled = true;
        status.workers = pool.NumWorkers;
        status.message = sprintf('enabled using %s with %d workers', class(pool), pool.NumWorkers);
    end
catch ME
    warning('generate_radar_uav_dataset:ParallelUnavailable', ...
        'Parallel pool startup failed (%s). Continuing with serial fallback.', ME.message);
    status.message = sprintf('startup failed: %s', ME.message);
end
end

function tf = can_use_parfor()
tf = license('test', 'Distrib_Computing_Toolbox') && ...
    (exist('parfor', 'builtin') > 0 || exist('parfor', 'file') > 0);
end

function tf = can_use_gpu()
tf = false;
if ~license('test', 'Distrib_Computing_Toolbox')
    return;
end
if exist('gpuDeviceCount', 'file') ~= 2 && exist('gpuDeviceCount', 'builtin') ~= 5
    return;
end
try
    tf = gpuDeviceCount('available') > 0;
catch
    tf = false;
end
end

function status = maybe_select_gpu(cfg)
status = struct('enabled', false, 'index', NaN, 'name', '', 'message', '');
if ~isfield(cfg, 'use_gpu_cost_field') || ~cfg.use_gpu_cost_field
    status.message = 'disabled by cfg.use_gpu_cost_field=false';
    return;
end
if ~can_use_gpu()
    status.message = 'unavailable: Parallel Computing Toolbox/GPU device not detected';
    warning('generate_radar_uav_dataset:GpuUnavailable', ...
        'GPU cost field requested but unavailable. Continuing with CPU cost field.');
    return;
end
try
    count = gpuDeviceCount('available');
    index = max(1, min(count, round(double(cfg.gpu_device_index))));
    gpu = gpuDevice(index);
    status.enabled = true;
    status.index = index;
    status.name = gpu.Name;
    status.message = sprintf('enabled using GPU %d: %s', index, gpu.Name);
catch ME
    status.message = sprintf('startup failed: %s', ME.message);
    warning('generate_radar_uav_dataset:GpuUnavailable', ...
        'GPU cost field startup failed (%s). Continuing with CPU cost field.', ME.message);
end
end

function val = rand_range(range_vals)
range_vals = double(range_vals);
val = range_vals(1) + rand() * (range_vals(end) - range_vals(1));
end

function val = get_cfg(cfg, name, default_val)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = default_val;
end
end
