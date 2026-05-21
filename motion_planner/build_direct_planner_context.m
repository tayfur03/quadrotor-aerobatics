function ctx = build_direct_planner_context(x_start, x_goal, env, params)
%BUILD_DIRECT_PLANNER_CONTEXT Shared environment/query setup for direct planners.
%
% Inputs
%   x_start  [3x1] world position [N;E;Alt]
%   x_goal   [3x1] world position [N;E;Alt]
%   env      direct-planner environment struct
%   params   planner options
%
% Output
%   ctx struct with:
%     .terrain_input
%     .terrain_query
%     .radar_los_map
%     .risk_query
%     .vis_threshold
%     .risk_limit
%     .bounds
%     .point_clearance
%     .env_meta
%     .points_clear_of_terrain
%     .evaluate_world_path_cost
%     .world_to_ned

if nargin < 4
    params = struct();
end

validateattributes(x_start, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_start', 1);
validateattributes(x_goal, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_goal', 2);
x_start = x_start(:);
x_goal = x_goal(:);

if ~isstruct(env)
    error('build_direct_planner_context:InvalidEnv', 'env must be a struct.');
end

point_clearance = max(get_opt(params, 'fsm_point_clearance', get_opt(params, 'min_clearance', 5.0)), 0.0);
max_clearance = get_opt(params, 'max_clearance', get_opt(params, 'max_agl', inf));
if isempty(max_clearance)
    max_clearance = inf;
end
max_clearance = max(max_clearance, point_clearance);
terrain_source = [];
terrain_fields = {'terrain', 'terrain_map', 'fsm_terrain_map'};
for i = 1:numel(terrain_fields)
    if isfield(env, terrain_fields{i}) && ~isempty(env.(terrain_fields{i}))
        terrain_source = env.(terrain_fields{i});
        break;
    end
end
if isempty(terrain_source)
    error('build_direct_planner_context:MissingTerrain', ...
        'env must include terrain or terrain_map.');
end

terrain_input = terrain_source_to_struct(terrain_source, env);
terrain_query = make_terrain_query(terrain_input);
env_meta = struct();
env_meta.terrain_bounds = [terrain_input.N_vec(1), terrain_input.N_vec(end), ...
                           terrain_input.E_vec(1), terrain_input.E_vec(end)];

if isfield(env, 'radar_los_map') && ~isempty(env.radar_los_map)
    radar_los_map = logical(env.radar_los_map);
elseif isfield(env, 'risk_grid') && ~isempty(env.risk_grid)
    radar_los_map = [];
else
    error('build_direct_planner_context:MissingThreat', ...
        'env must include radar_los_map or risk_grid.');
end

if isfield(env, 'visibility_threshold') && ~isempty(env.visibility_threshold)
    vis_threshold = env.visibility_threshold;
else
    vis_threshold = 0.5;
end

if ~isfield(env, 'N_vec') || ~isfield(env, 'E_vec') || ~isfield(env, 'alt_vec')
    error('build_direct_planner_context:MissingAxes', ...
        'env must include N_vec, E_vec, and alt_vec.');
end

N_vec = env.N_vec(:)';
E_vec = env.E_vec(:)';
alt_vec = env.alt_vec(:)';

if isfield(env, 'risk_altitude_mode') && ~isempty(env.risk_altitude_mode)
    risk_altitude_mode = lower(char(env.risk_altitude_mode));
elseif isfield(env, 'altitude_mode') && ~isempty(env.altitude_mode)
    risk_altitude_mode = lower(char(env.altitude_mode));
else
    risk_altitude_mode = 'asl';
end

if isfield(env, 'risk_grid') && ~isempty(env.risk_grid)
    V = env.risk_grid;
    if isequal(size(V), [numel(E_vec), numel(N_vec), numel(alt_vec)])
        if isfield(env, 'risk_extrapolation') && ~isempty(env.risk_extrapolation)
            risk_extrapolation = char(env.risk_extrapolation);
        else
            risk_extrapolation = 'none';
        end
        F = griddedInterpolant({E_vec, N_vec, alt_vec}, V, 'linear', risk_extrapolation);
        if strcmpi(risk_altitude_mode, 'agl')
            risk_query = @(P) risk_query_agl(P, F, terrain_query);
        else
            risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);
        end
        if isempty(radar_los_map)
            radar_los_map = V > vis_threshold;
        end
    else
        error('build_direct_planner_context:RiskGridSizeMismatch', ...
            'env.risk_grid size must match [numel(E_vec) x numel(N_vec) x numel(alt_vec)].');
    end
elseif ~isempty(radar_los_map)
    if isequal(size(radar_los_map), [numel(E_vec), numel(N_vec), numel(alt_vec)])
        V = double(radar_los_map);
        F = griddedInterpolant({E_vec, N_vec, alt_vec}, V, 'nearest', 'none');
        risk_query = @(P) reshape(F(P(2, :), P(1, :), P(3, :)), 1, []);
    else
        error('build_direct_planner_context:LOSSizeMismatch', ...
            'env.radar_los_map size must match [numel(E_vec) x numel(N_vec) x numel(alt_vec)].');
    end
end

if isfield(env, 'bounds_world') && ~isempty(env.bounds_world)
    bounds = env.bounds_world(:)';
else
    bounds = [N_vec(1), N_vec(end), E_vec(1), E_vec(end), alt_vec(1), alt_vec(end)];
end
bounds = complete_bounds(bounds, x_start, x_goal);

strict_zero_risk_tol = max(0, get_opt(params, 'strict_zero_risk_tol', 1e-6));
stealth_weight = get_opt(params, 'stealth_weight', 2.0);
altitude_cost_weight = get_opt(params, 'gamma', 0.0);
preferred_agl = get_opt(params, 'preferred_agl', 50.0);
if isinf(stealth_weight)
    risk_limit = strict_zero_risk_tol;
else
    risk_limit = vis_threshold;
end

env_meta.grid_size = [numel(N_vec), numel(E_vec), numel(alt_vec)];
env_meta.bounds = bounds;

ctx = struct();
ctx.terrain_input = terrain_input;
ctx.terrain_query = terrain_query;
ctx.radar_los_map = radar_los_map;
ctx.risk_query = risk_query;
ctx.vis_threshold = vis_threshold;
ctx.risk_limit = risk_limit;
ctx.bounds = bounds;
ctx.point_clearance = point_clearance;
ctx.max_clearance = max_clearance;
ctx.env_meta = env_meta;
ctx.altitude_cost_weight = altitude_cost_weight;
ctx.preferred_agl = preferred_agl;
ctx.points_clear_of_terrain = @(P) points_clear_of_terrain(P, terrain_query, point_clearance, max_clearance);
ctx.evaluate_world_path_cost = @(path_world) evaluate_world_path_cost( ...
    path_world, risk_query, risk_limit, stealth_weight, altitude_cost_weight, preferred_agl, point_clearance, max_clearance, terrain_query);
ctx.world_to_ned = @world_to_ned;
end

function risk = risk_query_agl(P, F, terrain_query)
if isempty(P)
    risk = zeros(1, 0);
    return;
end
terrain_h = terrain_query(P(1, :).', P(2, :).');
terrain_h = reshape(terrain_h, 1, []);
agl = P(3, :) - terrain_h;
risk = reshape(F(P(2, :), P(1, :), agl), 1, []);
end

function bounds = complete_bounds(bounds, x_start, x_goal)
bounds = bounds(:)';
bounds(1) = min([bounds(1), x_start(1), x_goal(1)]);
bounds(2) = max([bounds(2), x_start(1), x_goal(1)]);
bounds(3) = min([bounds(3), x_start(2), x_goal(2)]);
bounds(4) = max([bounds(4), x_start(2), x_goal(2)]);
bounds(5) = min([bounds(5), x_start(3), x_goal(3)]);
bounds(6) = max([bounds(6), x_start(3), x_goal(3)]);
end

function terrain_input = terrain_source_to_struct(src, env)
terrain_input = struct();
if isa(src, 'terrain_map')
    terrain_input.Z = src.Z;
    terrain_input.N_vec = src.N_vec(:)';
    terrain_input.E_vec = src.E_vec(:)';
    terrain_input.dx = mean(diff(src.N_vec));
    terrain_input.dy = mean(diff(src.E_vec));
elseif isstruct(src)
    terrain_input = src;
    if ~isfield(terrain_input, 'N_vec') && isfield(terrain_input, 'x_vec')
        terrain_input.N_vec = terrain_input.x_vec(:)';
    end
    if ~isfield(terrain_input, 'E_vec') && isfield(terrain_input, 'y_vec')
        terrain_input.E_vec = terrain_input.y_vec(:)';
    end
    if ~isfield(terrain_input, 'dx') && isfield(terrain_input, 'N_vec')
        terrain_input.dx = mean(diff(terrain_input.N_vec));
    end
    if ~isfield(terrain_input, 'dy') && isfield(terrain_input, 'E_vec')
        terrain_input.dy = mean(diff(terrain_input.E_vec));
    end
else
    error('build_direct_planner_context:InvalidTerrainContext', ...
        'Unsupported terrain source type: %s', class(src));
end

if isfield(env, 'alt_vec') && ~isempty(env.alt_vec)
    terrain_input.alt_vec = env.alt_vec(:)';
end
terrain_input.x_vec = terrain_input.N_vec(:)';
terrain_input.y_vec = terrain_input.E_vec(:)';
end

function F = make_terrain_query(terrain_input)
N_vec = terrain_input.N_vec(:)';
E_vec = terrain_input.E_vec(:)';
Z = terrain_input.Z;
if isequal(size(Z), [numel(E_vec), numel(N_vec)])
    Z_NE = double(Z.');
elseif isequal(size(Z), [numel(N_vec), numel(E_vec)])
    Z_NE = double(Z);
else
    error('build_direct_planner_context:TerrainSizeMismatch', ...
        'terrain_input.Z size does not match terrain axes.');
end
F = griddedInterpolant({N_vec, E_vec}, Z_NE, 'linear', 'nearest');
end

function tf = points_clear_of_terrain(P, terrain_query, clearance, max_clearance)
if isempty(P)
    tf = false(1, 0);
    return;
end
h = terrain_query(P(1, :).', P(2, :).');
h = reshape(h, 1, []);
agl = P(3, :) - h;
tf = isfinite(h) & (agl >= clearance) & (agl <= max_clearance);
end

function cost = evaluate_world_path_cost(path_world, risk_query, risk_limit, stealth_weight, altitude_cost_weight, preferred_agl, point_clearance, max_clearance, terrain_query)
if isempty(path_world) || size(path_world, 2) < 2
    cost = Inf;
    return;
end

cost = 0;
for i = 1:(size(path_world, 2) - 1)
    p0 = path_world(:, i);
    p1 = path_world(:, i + 1);
    seg_mid = 0.5 * (p0 + p1);
    seg_len = norm(p1 - p0);
    if ~all(points_clear_of_terrain([p0, p1, seg_mid], terrain_query, point_clearance, max_clearance))
        cost = Inf;
        return;
    end
    risk_mid = risk_query(seg_mid);
    if ~isfinite(risk_mid) || risk_mid > risk_limit
        cost = Inf;
        return;
    end
    cost_alt = compute_altitude_cost(seg_mid, terrain_query, altitude_cost_weight, preferred_agl, seg_len);
    if isinf(stealth_weight)
        cost = cost + seg_len + cost_alt;
    else
        cost = cost + seg_len + stealth_weight * seg_len * max(risk_mid, 0) + cost_alt;
    end
end
end

function cost_alt = compute_altitude_cost(P, terrain_query, altitude_cost_weight, preferred_agl, seg_len)
if altitude_cost_weight <= 0
    cost_alt = 0;
    return;
end
h = terrain_query(P(1, :).', P(2, :).');
agl = P(3, :)' - h;
agl_scale = max(preferred_agl, 1);
err_norm = (agl - preferred_agl) ./ agl_scale;
cost_alt = altitude_cost_weight * seg_len * mean(err_norm .^ 2);
end

function x_ned = world_to_ned(x_world)
x_ned = [x_world(1); x_world(2); -x_world(3)];
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
