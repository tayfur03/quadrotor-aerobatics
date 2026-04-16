function [T_fwd, T_bwd, F_field, meta] = compute_fmm_heuristic( ...
    terrain_map, radar_los_map, x_start, x_goal, params)
%COMPUTE_FMM_HEURISTIC Build forward/backward FMM heuristics on a 3D grid.
%
% Inputs
%   terrain_map   struct or terrain_map-like object carrying terrain heights
%   radar_los_map logical 3D grid, true where visible to at least one radar
%   x_start       [3x1] NED position [N;E;D], D negative-down
%   x_goal        [3x1] NED position [N;E;D], D negative-down
%   params        optional fields:
%     .w_mask         visible-cell penalty weight       (default 1.0)
%     .w_alt          altitude penalty weight           (default 0.1)
%     .h_opt          preferred AGL altitude [m]        (default 50)
%     .h_scale        altitude penalty scale [m]        (default 200)
%     .dz             vertical spacing [m]              (default 50)
%     .eps0           cost floor                        (default 1e-3)
%     .F_obs          obstacle speed floor              (default 1e-4)
%     .mem_limit_gb   allocation guard [GB]             (default 2.0)
%     .terrain_margin obstacle margin above terrain [m] (default 5.0)
%
% Outputs
%   T_fwd        [nN x nE x nZ] arrival time from x_start
%   T_bwd        [nN x nE x nZ] arrival time from x_goal
%   F_field      [nN x nE x nZ] speed field
%   meta         metadata and warm-start output

if nargin < 5
    params = struct();
end

validateattributes(x_start, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_start', 3);
validateattributes(x_goal, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'x_goal', 4);
x_start = x_start(:);
x_goal = x_goal(:);

assert(x_start(3) < 0, ...
    'compute_fmm_heuristic:InvalidNEDStart', ...
    'x_start must use NED with negative-down D (expected x_start(3) < 0).');
assert(x_goal(3) < 0, ...
    'compute_fmm_heuristic:InvalidNEDGoal', ...
    'x_goal must use NED with negative-down D (expected x_goal(3) < 0).');

w_mask = get_opt(params, 'w_mask', 1.0);
w_alt = get_opt(params, 'w_alt', 0.1);
h_opt = get_opt(params, 'h_opt', 50.0);
h_scale = max(get_opt(params, 'h_scale', 200.0), 1e-9);
dz = max(get_opt(params, 'dz', 50.0), 1e-6);
eps0 = max(get_opt(params, 'eps0', 1e-3), 1e-9);
F_obs = max(get_opt(params, 'F_obs', 1e-4), 1e-9);
mem_limit_gb = max(get_opt(params, 'mem_limit_gb', 2.0), 0.05);
terrain_margin = max(get_opt(params, 'terrain_margin', 5.0), 0.0);

[grid_x, grid_y, grid_z, dx, dy, terrain_NE, terrain_query] = ...
    normalize_terrain_inputs(terrain_map, x_start, x_goal, dz);

nN = numel(grid_x);
nE = numel(grid_y);
nZ = numel(grid_z);

meta = struct();
meta.grid_x = grid_x;
meta.grid_y = grid_y;
meta.grid_z = grid_z;
meta.dx = dx;
meta.dy = dy;
meta.dz = dz;
meta.fmm_time_s = 0;
meta.mem_used_gb = 0;
meta.fallback = false;
meta.warmstart_path = [];
meta.warmstart_cost = Inf;

radar_los_NE = normalize_los_map(radar_los_map, nN, nE, nZ);
[grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, crop_meta] = ...
    crop_fmm_domain(grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, x_start, x_goal, params);
terrain_query = griddedInterpolant({grid_x, grid_y}, terrain_NE, 'linear', 'nearest');

nN = numel(grid_x);
nE = numel(grid_y);
nZ = numel(grid_z);
meta.grid_x = grid_x;
meta.grid_y = grid_y;
meta.grid_z = grid_z;
meta.crop_indices = crop_meta;

mem_gb = 3 * double(nN) * double(nE) * double(nZ) * 4 / 1e9;
if mem_gb > mem_limit_gb
    warning('compute_fmm_heuristic:MemoryLimit', ...
        ['FMM grid would allocate %.2f GB (> %.2f GB limit). ', ...
         'Returning empty heuristic and falling back to Euclidean.'], ...
        mem_gb, mem_limit_gb);
    T_fwd = [];
    T_bwd = [];
    F_field = [];
    meta.mem_used_gb = mem_gb;
    meta.fallback = true;
    return;
end

[N_grid, E_grid] = ndgrid(grid_x, grid_y);
terrain_slice = terrain_query(N_grid, E_grid);
terrain_slice(~isfinite(terrain_slice)) = max(terrain_NE(:));

F_field = zeros(nN, nE, nZ, 'single');
cost_field = zeros(nN, nE, nZ, 'single');

for k = 1:nZ
    altitude_msl = grid_z(k);
    agl = altitude_msl - terrain_slice;
    alt_pen = max(0.0, agl - h_opt) / h_scale;

    if isempty(radar_los_NE)
        vis_pen = zeros(nN, nE);
    else
        masked = ~radar_los_NE(:, :, k);
        vis_pen = 1.0 - double(masked);
    end

    c_slice = w_mask .* vis_pen + w_alt .* alt_pen + eps0;
    obstacle_mask = altitude_msl <= (terrain_slice + terrain_margin);
    c_slice(obstacle_mask) = 1.0 / F_obs;

    cost_field(:, :, k) = single(c_slice);
    F_field(:, :, k) = single(1.0 ./ c_slice);
end

% Map NED [N;E;D] to FMM grid indices:
%   i = round((x_N - grid_x(1))/dx) + 1
%   j = round((x_E - grid_y(1))/dy) + 1
%   k = round((-x_D - grid_z(1))/dz) + 1
start_idx = ned_to_grid_index(x_start, grid_x, grid_y, grid_z, dx, dy, dz);
goal_idx = ned_to_grid_index(x_goal, grid_x, grid_y, grid_z, dx, dy, dz);

t_fmm = tic;
T_fwd = run_fmm_pass(F_field, start_idx, dx, dy, dz);
T_bwd = run_fmm_pass(F_field, goal_idx, dx, dy, dz);
meta.fmm_time_s = toc(t_fmm);

w = whos('T_fwd', 'T_bwd', 'F_field');
meta.mem_used_gb = sum([w.bytes]) / 1e9;

[warm_path_idx, ok_path] = descend_warmstart(T_fwd, goal_idx, start_idx, F_field, terrain_slice, grid_z, terrain_margin);
if ok_path
    meta.warmstart_path = idx_path_to_ned(warm_path_idx, grid_x, grid_y, grid_z);
    meta.warmstart_cost = integrate_path_cost(meta.warmstart_path, cost_field, meta);
else
    meta.warmstart_path = [];
    meta.warmstart_cost = Inf;
end
end

function [grid_x, grid_y, grid_z, dx, dy, terrain_NE, terrain_query] = ...
        normalize_terrain_inputs(terrain_map, x_start, x_goal, dz)
% Returns terrain on [N x E] grids and a griddedInterpolant F(N,E).

if isa(terrain_map, 'terrain_map')
    grid_x = terrain_map.N_vec(:)';
    grid_y = terrain_map.E_vec(:)';
    Z_src = terrain_map.Z;
elseif isstruct(terrain_map)
    if isfield(terrain_map, 'N_vec')
        grid_x = terrain_map.N_vec(:)';
    elseif isfield(terrain_map, 'x_vec')
        grid_x = terrain_map.x_vec(:)';
    else
        error('compute_fmm_heuristic:MissingNorthAxis', ...
            'terrain_map must include N_vec or x_vec.');
    end

    if isfield(terrain_map, 'E_vec')
        grid_y = terrain_map.E_vec(:)';
    elseif isfield(terrain_map, 'y_vec')
        grid_y = terrain_map.y_vec(:)';
    else
        error('compute_fmm_heuristic:MissingEastAxis', ...
            'terrain_map must include E_vec or y_vec.');
    end

    if ~isfield(terrain_map, 'Z') || isempty(terrain_map.Z)
        error('compute_fmm_heuristic:MissingTerrainZ', ...
            'terrain_map must include a non-empty Z field.');
    end
    Z_src = terrain_map.Z;
else
    error('compute_fmm_heuristic:InvalidTerrainInput', ...
        'terrain_map must be a struct or terrain_map object.');
end

if any(diff(grid_x) <= 0) || any(diff(grid_y) <= 0)
    error('compute_fmm_heuristic:NonMonotonicAxes', ...
        'Terrain axes must be strictly increasing.');
end

if isa(terrain_map, 'terrain_map')
    dx = mean(diff(grid_x));
    dy = mean(diff(grid_y));
else
    dx = get_struct_opt(terrain_map, {'dx'}, mean(diff(grid_x)));
    dy = get_struct_opt(terrain_map, {'dy'}, mean(diff(grid_y)));
end
dx = max(abs(dx), 1e-6);
dy = max(abs(dy), 1e-6);

if isequal(size(Z_src), [numel(grid_y), numel(grid_x)])
    terrain_NE = double(Z_src.');
elseif isequal(size(Z_src), [numel(grid_x), numel(grid_y)])
    terrain_NE = double(Z_src);
else
    error('compute_fmm_heuristic:TerrainSizeMismatch', ...
        'terrain_map.Z size does not match the horizontal axes.');
end

terrain_query = griddedInterpolant({grid_x, grid_y}, terrain_NE, 'linear', 'nearest');

if isa(terrain_map, 'terrain_map')
    grid_z = derive_alt_grid([], terrain_NE, x_start, x_goal, dz);
else
    grid_z = get_struct_opt(terrain_map, {'alt_vec', 'grid_z', 'z_vec'}, []);
    if isempty(grid_z)
        grid_z = derive_alt_grid([], terrain_NE, x_start, x_goal, dz);
    else
        grid_z = grid_z(:)';
    end
end

if any(diff(grid_z) <= 0)
    error('compute_fmm_heuristic:NonMonotonicAltitudeAxis', ...
        'Altitude axis must be strictly increasing.');
end
end

function grid_z = derive_alt_grid(grid_z_hint, terrain_NE, x_start, x_goal, dz)
if nargin >= 1 && ~isempty(grid_z_hint)
    grid_z = grid_z_hint(:)';
    return;
end

terrain_min = min(terrain_NE(:));
terrain_max = max(terrain_NE(:));
start_alt = -x_start(3);
goal_alt = -x_goal(3);
alt_min = floor(min([terrain_min, start_alt, goal_alt]) / dz) * dz;
alt_max = ceil(max([terrain_max + 300, start_alt, goal_alt]) / dz) * dz;
if alt_max <= alt_min
    alt_max = alt_min + dz;
end
grid_z = alt_min:dz:alt_max;
if numel(grid_z) < 2
    grid_z = [alt_min, alt_min + dz];
end
end

function radar_los_NE = normalize_los_map(radar_los_map, nN, nE, nZ)
if isempty(radar_los_map)
    radar_los_NE = [];
    return;
end

if ~islogical(radar_los_map)
    radar_los_map = radar_los_map ~= 0;
end

sz = size(radar_los_map);
if isequal(sz, [nN, nE, nZ])
    radar_los_NE = radar_los_map;
elseif isequal(sz, [nE, nN, nZ])
    radar_los_NE = permute(radar_los_map, [2, 1, 3]);
else
    warning('compute_fmm_heuristic:LOSSizeMismatch', ...
        ['radar_los_map size [%d %d %d] does not match either [nN nE nZ] ', ...
         'or [nE nN nZ]. Ignoring radar masking.'], sz(1), sz(2), sz(3));
    radar_los_NE = [];
end
end

function idx = ned_to_grid_index(x_ned, grid_x, grid_y, grid_z, dx, dy, dz)
altitude_msl = -x_ned(3);
i = round((x_ned(1) - grid_x(1)) / dx) + 1;
j = round((x_ned(2) - grid_y(1)) / dy) + 1;
k = round((altitude_msl - grid_z(1)) / dz) + 1;
idx = [clamp_index(i, numel(grid_x)); ...
       clamp_index(j, numel(grid_y)); ...
       clamp_index(k, numel(grid_z))];
end

function T = run_fmm_pass(F_field, seed_idx, dx, dy, dz)
sz = size(F_field);
T = inf(sz, 'single');
accepted = false(sz);

T(seed_idx(1), seed_idx(2), seed_idx(3)) = 0;
max_heap = max(1024, ceil(prod(sz) * 0.02));
heap_key = inf(max_heap, 1);
heap_i = zeros(max_heap, 1);
heap_j = zeros(max_heap, 1);
heap_k = zeros(max_heap, 1);
heap_size = 0;
[heap_key, heap_i, heap_j, heap_k, heap_size] = ...
    heap_push(heap_key, heap_i, heap_j, heap_k, heap_size, 0.0, seed_idx(1), seed_idx(2), seed_idx(3));

while heap_size > 0
    [~, i, j, k, heap_key, heap_i, heap_j, heap_k, heap_size] = ...
        heap_pop(heap_key, heap_i, heap_j, heap_k, heap_size);

    if accepted(i, j, k)
        continue;
    end

    if ~isfinite(T(i, j, k))
        continue;
    end

    accepted(i, j, k) = true;
    neigh = [i-1, j,   k;
             i+1, j,   k;
             i,   j-1, k;
             i,   j+1, k;
             i,   j,   k-1;
             i,   j,   k+1];

    for n = 1:size(neigh, 1)
        ii = neigh(n, 1);
        jj = neigh(n, 2);
        kk = neigh(n, 3);
        if ii < 1 || ii > sz(1) || jj < 1 || jj > sz(2) || kk < 1 || kk > sz(3)
            continue;
        end

        if accepted(ii, jj, kk)
            continue;
        end

        T_new = solve_local_update(T, accepted, ii, jj, kk, double(F_field(ii, jj, kk)), dx, dy, dz);
        if T_new + 1e-9 < double(T(ii, jj, kk))
            T(ii, jj, kk) = single(T_new);
            [heap_key, heap_i, heap_j, heap_k, heap_size] = ...
                heap_push(heap_key, heap_i, heap_j, heap_k, heap_size, T_new, ii, jj, kk);
        end
    end
end
end

function [grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, crop_meta] = ...
        crop_fmm_domain(grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, x_start, x_goal, params)
crop_meta = struct('i_idx', 1:numel(grid_x), 'j_idx', 1:numel(grid_y), 'k_idx', 1:numel(grid_z));

if isfield(params, 'fmm_bounds_world') && ~isempty(params.fmm_bounds_world)
    b = params.fmm_bounds_world(:)';
else
    return;
end

pad_xy = get_opt(params, 'fmm_bounds_pad_xy', 3 * max(mean(diff(grid_x)), mean(diff(grid_y))));
pad_z = get_opt(params, 'fmm_bounds_pad_z', 2 * mean(diff(grid_z)));

N_min = min([b(1), x_start(1), x_goal(1)]) - pad_xy;
N_max = max([b(2), x_start(1), x_goal(1)]) + pad_xy;
E_min = min([b(3), x_start(2), x_goal(2)]) - pad_xy;
E_max = max([b(4), x_start(2), x_goal(2)]) + pad_xy;
A_start = -x_start(3);
A_goal = -x_goal(3);
Z_min = min([b(5), A_start, A_goal]) - pad_z;
Z_max = max([b(6), A_start, A_goal]) + pad_z;

i_idx = find(grid_x >= N_min & grid_x <= N_max);
j_idx = find(grid_y >= E_min & grid_y <= E_max);
k_idx = find(grid_z >= Z_min & grid_z <= Z_max);

if isempty(i_idx), i_idx = 1:numel(grid_x); end
if isempty(j_idx), j_idx = 1:numel(grid_y); end
if isempty(k_idx), k_idx = 1:numel(grid_z); end

grid_x = grid_x(i_idx);
grid_y = grid_y(j_idx);
grid_z = grid_z(k_idx);
terrain_NE = terrain_NE(i_idx, j_idx);
if ~isempty(radar_los_NE)
    radar_los_NE = radar_los_NE(i_idx, j_idx, k_idx);
end

crop_meta.i_idx = i_idx;
crop_meta.j_idx = j_idx;
crop_meta.k_idx = k_idx;
end

function [heap_key, heap_i, heap_j, heap_k, heap_size] = ...
        heap_push(heap_key, heap_i, heap_j, heap_k, heap_size, key, i, j, k)
heap_size = heap_size + 1;
if heap_size > numel(heap_key)
    new_cap = max(heap_size, 2 * numel(heap_key));
    heap_key(new_cap, 1) = inf;
    heap_i(new_cap, 1) = 0;
    heap_j(new_cap, 1) = 0;
    heap_k(new_cap, 1) = 0;
end

idx = heap_size;
heap_key(idx) = key;
heap_i(idx) = i;
heap_j(idx) = j;
heap_k(idx) = k;

while idx > 1
    parent = floor(idx / 2);
    if heap_key(parent) <= heap_key(idx)
        break;
    end
    [heap_key(parent), heap_key(idx)] = deal(heap_key(idx), heap_key(parent));
    [heap_i(parent), heap_i(idx)] = deal(heap_i(idx), heap_i(parent));
    [heap_j(parent), heap_j(idx)] = deal(heap_j(idx), heap_j(parent));
    [heap_k(parent), heap_k(idx)] = deal(heap_k(idx), heap_k(parent));
    idx = parent;
end
end

function [key, i, j, k, heap_key, heap_i, heap_j, heap_k, heap_size] = ...
        heap_pop(heap_key, heap_i, heap_j, heap_k, heap_size)
key = heap_key(1);
i = heap_i(1);
j = heap_j(1);
k = heap_k(1);

heap_key(1) = heap_key(heap_size);
heap_i(1) = heap_i(heap_size);
heap_j(1) = heap_j(heap_size);
heap_k(1) = heap_k(heap_size);
heap_size = heap_size - 1;

idx = 1;
while true
    left = 2 * idx;
    right = left + 1;
    if left > heap_size
        break;
    end

    smallest = left;
    if right <= heap_size && heap_key(right) < heap_key(left)
        smallest = right;
    end
    if heap_key(idx) <= heap_key(smallest)
        break;
    end

    [heap_key(idx), heap_key(smallest)] = deal(heap_key(smallest), heap_key(idx));
    [heap_i(idx), heap_i(smallest)] = deal(heap_i(smallest), heap_i(idx));
    [heap_j(idx), heap_j(smallest)] = deal(heap_j(smallest), heap_j(idx));
    [heap_k(idx), heap_k(smallest)] = deal(heap_k(smallest), heap_k(idx));
    idx = smallest;
end
end

function T_new = solve_local_update(T, accepted, i, j, k, speed, dx, dy, dz)
if speed <= 0 || ~isfinite(speed)
    T_new = inf;
    return;
end

vals = [min_neighbor(T, accepted, i-1, j, k, i+1, j, k), ...
        min_neighbor(T, accepted, i, j-1, k, i, j+1, k), ...
        min_neighbor(T, accepted, i, j, k-1, i, j, k+1)];
h = [dx, dy, dz];
valid = isfinite(vals);

if ~any(valid)
    T_new = inf;
    return;
end

slowness = 1.0 / speed;
axis_idx = find(valid);

while ~isempty(axis_idx)
    vals_use = vals(axis_idx);
    h_use = h(axis_idx);
    inv_h2 = 1.0 ./ (h_use .^ 2);

    if numel(axis_idx) == 1
        T_new = vals_use(1) + slowness * h_use(1);
        return;
    end

    A_coef = sum(inv_h2);
    B_coef = sum(vals_use .* inv_h2);
    C_coef = sum((vals_use .^ 2) .* inv_h2) - slowness^2;
    disc = B_coef^2 - A_coef * C_coef;

    if disc >= 0
        T_candidate = (B_coef + sqrt(disc)) / A_coef;
        if T_candidate > max(vals_use) + 1e-12
            T_new = T_candidate;
            return;
        end
    end

    contrib = (vals_use .^ 2) .* inv_h2;
    [~, drop_local] = max(contrib);
    axis_idx(drop_local) = [];
end

T_new = inf;
end

function val = min_neighbor(T, accepted, i1, j1, k1, i2, j2, k2)
sz = size(T);
vals = inf(1, 2);

if i1 >= 1 && i1 <= sz(1) && j1 >= 1 && j1 <= sz(2) && k1 >= 1 && k1 <= sz(3)
    if accepted(i1, j1, k1)
        vals(1) = T(i1, j1, k1);
    end
end
if i2 >= 1 && i2 <= sz(1) && j2 >= 1 && j2 <= sz(2) && k2 >= 1 && k2 <= sz(3)
    if accepted(i2, j2, k2)
        vals(2) = T(i2, j2, k2);
    end
end

val = min(vals);
end

function [path_idx, ok] = descend_warmstart(T_fwd, goal_idx, start_idx, F_field, terrain_slice, grid_z, terrain_margin)
sz = size(T_fwd);
max_steps = 4 * prod(sz);
path_idx = zeros(3, max_steps);
path_idx(:, 1) = goal_idx(:);
count = 1;
curr = goal_idx(:);
ok = false;

for step = 2:max_steps
    if all(abs(curr - start_idx(:)) <= 0)
        ok = true;
        break;
    end

    [next_idx, found] = best_descending_neighbor(T_fwd, curr, F_field, terrain_slice, grid_z, terrain_margin);
    if ~found
        break;
    end
    if T_fwd(next_idx(1), next_idx(2), next_idx(3)) >= T_fwd(curr(1), curr(2), curr(3)) - 1e-9
        break;
    end

    count = count + 1;
    path_idx(:, count) = next_idx;
    curr = next_idx;

    if all(curr == start_idx(:))
        ok = true;
        break;
    end
end

path_idx = path_idx(:, 1:count);
if ok
    path_idx = fliplr(path_idx);
end
end

function [next_idx, found] = best_descending_neighbor(T_fwd, curr, F_field, terrain_slice, grid_z, terrain_margin)
sz = size(T_fwd);
best_val = inf;
next_idx = curr(:);
found = false;

for di = -1:1
    for dj = -1:1
        for dk = -1:1
            if di == 0 && dj == 0 && dk == 0
                continue;
            end
            ii = curr(1) + di;
            jj = curr(2) + dj;
            kk = curr(3) + dk;
            if ii < 1 || ii > sz(1) || jj < 1 || jj > sz(2) || kk < 1 || kk > sz(3)
                continue;
            end
            if ~isfinite(T_fwd(ii, jj, kk)) || F_field(ii, jj, kk) <= 0
                continue;
            end
            if grid_z(kk) <= terrain_slice(ii, jj) + terrain_margin
                continue;
            end
            val = double(T_fwd(ii, jj, kk));
            if val < best_val
                best_val = val;
                next_idx = [ii; jj; kk];
                found = true;
            end
        end
    end
end
end

function path_ned = idx_path_to_ned(path_idx, grid_x, grid_y, grid_z)
n_pts = size(path_idx, 2);
path_ned = zeros(3, n_pts);
for m = 1:n_pts
    ii = path_idx(1, m);
    jj = path_idx(2, m);
    kk = path_idx(3, m);
    path_ned(:, m) = [grid_x(ii); grid_y(jj); -grid_z(kk)];
end
end

function total_cost = integrate_path_cost(path_ned, cost_field, meta)
if size(path_ned, 2) < 2
    total_cost = Inf;
    return;
end

total_cost = 0;
for k = 1:(size(path_ned, 2) - 1)
    p0 = path_ned(:, k);
    p1 = path_ned(:, k + 1);
    seg_len = norm(p1 - p0);
    c0 = sample_cost(cost_field, p0, meta);
    c1 = sample_cost(cost_field, p1, meta);
    if ~isfinite(c0) || ~isfinite(c1)
        total_cost = Inf;
        return;
    end
    total_cost = total_cost + 0.5 * (c0 + c1) * seg_len;
end
end

function c = sample_cost(cost_field, x_ned, meta)
alt = -x_ned(3);
fx = (x_ned(1) - meta.grid_x(1)) / meta.dx + 1;
fy = (x_ned(2) - meta.grid_y(1)) / meta.dy + 1;
fz = (alt - meta.grid_z(1)) / meta.dz + 1;
c = trilinear_sample(cost_field, fx, fy, fz);
end

function val = trilinear_sample(V, fx, fy, fz)
sz = size(V);
if fx < 1 || fy < 1 || fz < 1 || ...
   fx > sz(1) || fy > sz(2) || fz > sz(3)
    val = Inf;
    return;
end

i0 = max(1, min(sz(1) - 1, floor(fx)));
j0 = max(1, min(sz(2) - 1, floor(fy)));
k0 = max(1, min(sz(3) - 1, floor(fz)));
i1 = min(sz(1), i0 + 1);
j1 = min(sz(2), j0 + 1);
k1 = min(sz(3), k0 + 1);

tx = max(0, min(1, fx - i0));
ty = max(0, min(1, fy - j0));
tz = max(0, min(1, fz - k0));

c000 = double(V(i0, j0, k0));
c100 = double(V(i1, j0, k0));
c010 = double(V(i0, j1, k0));
c110 = double(V(i1, j1, k0));
c001 = double(V(i0, j0, k1));
c101 = double(V(i1, j0, k1));
c011 = double(V(i0, j1, k1));
c111 = double(V(i1, j1, k1));

c00 = (1 - tx) * c000 + tx * c100;
c10 = (1 - tx) * c010 + tx * c110;
c01 = (1 - tx) * c001 + tx * c101;
c11 = (1 - tx) * c011 + tx * c111;
c0 = (1 - ty) * c00 + ty * c10;
c1 = (1 - ty) * c01 + ty * c11;
val = (1 - tz) * c0 + tz * c1;
end

function idx = clamp_index(idx, n)
idx = max(1, min(n, idx));
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end

function val = get_struct_opt(S, names, default)
val = default;
for i = 1:numel(names)
    if isfield(S, names{i}) && ~isempty(S.(names{i}))
        val = S.(names{i});
        return;
    end
end
end
