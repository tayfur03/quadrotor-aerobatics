function [T_fwd, T_bwd, F_field, meta] = compute_fsm_heuristic( ...
    terrain_map, radar_los_map, x_start, x_goal, params)
%COMPUTE_FSM_HEURISTIC Build forward/backward anisotropic FSM heuristics.
%
% Inputs
%   terrain_map   struct or terrain_map-like object with terrain heights
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
%     .max_sweeps     maximum FSM sweep iterations      (default 48)
%     .tol            convergence tolerance             (default 1e-3)
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
    'compute_fsm_heuristic:InvalidNEDStart', ...
    'x_start must use NED with negative-down D (expected x_start(3) < 0).');
assert(x_goal(3) < 0, ...
    'compute_fsm_heuristic:InvalidNEDGoal', ...
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
max_sweeps = max(4, round(get_opt(params, 'max_sweeps', 48)));
tol = max(get_opt(params, 'tol', 1e-3), 1e-9);

[grid_x, grid_y, grid_z, dx, dy, terrain_NE] = normalize_terrain_inputs(terrain_map, x_start, x_goal, dz);
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
meta.fsm_time_s = 0;
meta.mem_used_gb = 0;
meta.fallback = false;
meta.warmstart_path = [];
meta.warmstart_cost = Inf;
meta.max_sweeps = max_sweeps;
meta.tol = tol;
meta.parallel_used = false;

radar_los_NE = normalize_los_map(radar_los_map, nN, nE, nZ);
[grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, crop_meta] = ...
    crop_fsm_domain(grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, x_start, x_goal, params);

nN = numel(grid_x);
nE = numel(grid_y);
nZ = numel(grid_z);
meta.grid_x = grid_x;
meta.grid_y = grid_y;
meta.grid_z = grid_z;
meta.crop_indices = crop_meta;

mem_gb = 3 * double(nN) * double(nE) * double(nZ) * 4 / 1e9;
if mem_gb > mem_limit_gb
    warning('compute_fsm_heuristic:MemoryLimit', ...
        ['FSM grid would allocate %.2f GB (> %.2f GB limit). ', ...
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
terrain_slice = interp2(grid_y, grid_x, terrain_NE, E_grid, N_grid, 'linear', max(terrain_NE(:)));

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

start_idx = ned_to_grid_index(x_start, grid_x, grid_y, grid_z, dx, dy, dz);
goal_idx = ned_to_grid_index(x_goal, grid_x, grid_y, grid_z, dx, dy, dz);

t_fsm = tic;
[T_fwd, sweeps_fwd, delta_fwd, T_bwd, sweeps_bwd, delta_bwd, parallel_used] = ...
    run_fsm_dual_pass(F_field, start_idx, goal_idx, dx, dy, dz, max_sweeps, tol);
meta.fsm_time_s = toc(t_fsm);
meta.sweeps_fwd = sweeps_fwd;
meta.sweeps_bwd = sweeps_bwd;
meta.final_delta_fwd = delta_fwd;
meta.final_delta_bwd = delta_bwd;
meta.parallel_used = parallel_used;

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

function [grid_x, grid_y, grid_z, dx, dy, terrain_NE] = normalize_terrain_inputs(terrain_map, x_start, x_goal, dz)
if isa(terrain_map, 'terrain_map')
    grid_x = terrain_map.N_vec(:)';
    grid_y = terrain_map.E_vec(:)';
    Z_src = terrain_map.Z;
    dx = mean(diff(grid_x));
    dy = mean(diff(grid_y));
elseif isstruct(terrain_map)
    if isfield(terrain_map, 'N_vec')
        grid_x = terrain_map.N_vec(:)';
    elseif isfield(terrain_map, 'x_vec')
        grid_x = terrain_map.x_vec(:)';
    else
        error('compute_fsm_heuristic:MissingNorthAxis', 'terrain_map must include N_vec or x_vec.');
    end

    if isfield(terrain_map, 'E_vec')
        grid_y = terrain_map.E_vec(:)';
    elseif isfield(terrain_map, 'y_vec')
        grid_y = terrain_map.y_vec(:)';
    else
        error('compute_fsm_heuristic:MissingEastAxis', 'terrain_map must include E_vec or y_vec.');
    end

    if ~isfield(terrain_map, 'Z') || isempty(terrain_map.Z)
        error('compute_fsm_heuristic:MissingTerrainZ', 'terrain_map must include a non-empty Z field.');
    end
    Z_src = terrain_map.Z;
    dx = get_struct_opt(terrain_map, {'dx'}, mean(diff(grid_x)));
    dy = get_struct_opt(terrain_map, {'dy'}, mean(diff(grid_y)));
else
    error('compute_fsm_heuristic:InvalidTerrainInput', ...
        'terrain_map must be a struct or terrain_map object.');
end

if isequal(size(Z_src), [numel(grid_y), numel(grid_x)])
    terrain_NE = double(Z_src.');
elseif isequal(size(Z_src), [numel(grid_x), numel(grid_y)])
    terrain_NE = double(Z_src);
else
    error('compute_fsm_heuristic:TerrainSizeMismatch', ...
        'terrain_map.Z size does not match the horizontal axes.');
end

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
    warning('compute_fsm_heuristic:LOSSizeMismatch', ...
        ['radar_los_map size [%d %d %d] does not match either [nN nE nZ] ', ...
         'or [nE nN nZ]. Ignoring radar masking.'], sz(1), sz(2), sz(3));
    radar_los_NE = [];
end
end

function [grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, crop_meta] = ...
        crop_fsm_domain(grid_x, grid_y, grid_z, terrain_NE, radar_los_NE, x_start, x_goal, params)
crop_meta = struct('i_idx', 1:numel(grid_x), 'j_idx', 1:numel(grid_y), 'k_idx', 1:numel(grid_z));

if isfield(params, 'fsm_bounds_world') && ~isempty(params.fsm_bounds_world)
    b = params.fsm_bounds_world(:)';
elseif isfield(params, 'bounds_world') && ~isempty(params.bounds_world)
    b = params.bounds_world(:)';
else
    return;
end

pad_xy = get_opt(params, 'fsm_bounds_pad_xy', 3 * max(mean(diff(grid_x)), mean(diff(grid_y))));
pad_z = get_opt(params, 'fsm_bounds_pad_z', 2 * mean(diff(grid_z)));

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

function idx = ned_to_grid_index(x_ned, grid_x, grid_y, grid_z, dx, dy, dz)
altitude_msl = -x_ned(3);
i = round((x_ned(1) - grid_x(1)) / dx) + 1;
j = round((x_ned(2) - grid_y(1)) / dy) + 1;
k = round((altitude_msl - grid_z(1)) / dz) + 1;
idx = [clamp_index(i, numel(grid_x)); ...
       clamp_index(j, numel(grid_y)); ...
       clamp_index(k, numel(grid_z))];
end

function [T_fwd, sweeps_fwd, delta_fwd, T_bwd, sweeps_bwd, delta_bwd, parallel_used] = ...
        run_fsm_dual_pass(F_field, start_idx, goal_idx, dx, dy, dz, max_sweeps, tol)
parallel_used = false;

if can_use_parallel_fsm()
    try
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('threads');
        end
        future_fwd = parfeval(pool, @run_fsm_pass, 3, F_field, start_idx, dx, dy, dz, max_sweeps, tol);
        future_bwd = parfeval(pool, @run_fsm_pass, 3, F_field, goal_idx, dx, dy, dz, max_sweeps, tol);
        [T_fwd, sweeps_fwd, delta_fwd] = fetchOutputs(future_fwd);
        [T_bwd, sweeps_bwd, delta_bwd] = fetchOutputs(future_bwd);
        parallel_used = true;
        return;
    catch ME
        warning('compute_fsm_heuristic:ParallelFallback', ...
            'Parallel FSM solve unavailable (%s). Reverting to serial passes.', ME.message);
    end
end

[T_fwd, sweeps_fwd, delta_fwd] = run_fsm_pass(F_field, start_idx, dx, dy, dz, max_sweeps, tol);
[T_bwd, sweeps_bwd, delta_bwd] = run_fsm_pass(F_field, goal_idx, dx, dy, dz, max_sweeps, tol);
end

function tf = can_use_parallel_fsm()
tf = false;
if exist('parfeval', 'file') ~= 2 || exist('gcp', 'file') ~= 2
    return;
end
try
    tf = license('test', 'Distrib_Computing_Toolbox') || ~isempty(gcp('nocreate'));
catch
    tf = false;
end
end

function [T, sweeps_done, last_delta] = run_fsm_pass(F_field, seed_idx, dx, dy, dz, max_sweeps, tol)
sz = size(F_field);
T = inf(sz, 'single');
T(seed_idx(1), seed_idx(2), seed_idx(3)) = 0;

dirs = [...
     1,  1,  1;
    -1,  1,  1;
     1, -1,  1;
    -1, -1,  1;
     1,  1, -1;
    -1,  1, -1;
     1, -1, -1;
    -1, -1, -1];

last_delta = inf;
sweeps_done = 0;
for s = 1:max_sweeps
    T_prev = T;
    for d = 1:size(dirs, 1)
        T = directional_update(T, F_field, seed_idx, dirs(d, :), dx, dy, dz);
    end

    dif = abs(T(:) - T_prev(:));
    last_delta = max(dif(isfinite(dif)));
    if isempty(last_delta) || ~isfinite(last_delta)
        last_delta = inf;
    end
    sweeps_done = s;
    if last_delta < tol
        break;
    end
end
end

function T_out = directional_update(T_in, F_field, seed_idx, dir_sign, dx, dy, dz)
T_work = T_in;

n1 = size(T_work, 1);
n2 = size(T_work, 2);
n3 = size(T_work, 3);
i_range = sweep_range(n1, dir_sign(1));
j_range = sweep_range(n2, dir_sign(2));
k_range = sweep_range(n3, dir_sign(3));
inv_dx2 = 1.0 / (dx * dx);
inv_dy2 = 1.0 / (dy * dy);
inv_dz2 = 1.0 / (dz * dz);

for k = k_range
    for j = j_range
        for i = i_range
            if i == seed_idx(1) && j == seed_idx(2) && k == seed_idx(3)
                T_work(i, j, k) = 0;
                continue;
            end

            speed = double(F_field(i, j, k));
            if speed <= 0 || ~isfinite(speed)
                continue;
            end

            a = inf;
            if i > 1
                a = double(T_work(i - 1, j, k));
            end
            if i < n1
                a = min(a, double(T_work(i + 1, j, k)));
            end

            b = inf;
            if j > 1
                b = double(T_work(i, j - 1, k));
            end
            if j < n2
                b = min(b, double(T_work(i, j + 1, k)));
            end

            c = inf;
            if k > 1
                c = double(T_work(i, j, k - 1));
            end
            if k < n3
                c = min(c, double(T_work(i, j, k + 1)));
            end

            valid_a = isfinite(a);
            valid_b = isfinite(b);
            valid_c = isfinite(c);
            if ~(valid_a || valid_b || valid_c)
                continue;
            end

            slowness2 = 1.0 / (speed * speed);
            t_candidate = inf;

            if valid_a && valid_b && valid_c
                A3 = inv_dx2 + inv_dy2 + inv_dz2;
                B3 = a * inv_dx2 + b * inv_dy2 + c * inv_dz2;
                C3 = a * a * inv_dx2 + b * b * inv_dy2 + c * c * inv_dz2 - slowness2;
                disc3 = B3 * B3 - A3 * C3;
                if disc3 >= 0
                    cand3 = (B3 + sqrt(disc3)) / A3;
                    if cand3 > max([a, b, c]) + 1e-12
                        t_candidate = cand3;
                    end
                end
            end

            if ~isfinite(t_candidate)
                if valid_a && valid_b
                    cand = solve_two_axis(a, b, inv_dx2, inv_dy2, slowness2);
                    if cand > max(a, b) + 1e-12
                        t_candidate = min(t_candidate, cand);
                    end
                end
                if valid_a && valid_c
                    cand = solve_two_axis(a, c, inv_dx2, inv_dz2, slowness2);
                    if cand > max(a, c) + 1e-12
                        t_candidate = min(t_candidate, cand);
                    end
                end
                if valid_b && valid_c
                    cand = solve_two_axis(b, c, inv_dy2, inv_dz2, slowness2);
                    if cand > max(b, c) + 1e-12
                        t_candidate = min(t_candidate, cand);
                    end
                end
            end

            if ~isfinite(t_candidate)
                if valid_a
                    t_candidate = min(t_candidate, a + dx / speed);
                end
                if valid_b
                    t_candidate = min(t_candidate, b + dy / speed);
                end
                if valid_c
                    t_candidate = min(t_candidate, c + dz / speed);
                end
            end

            if t_candidate < double(T_work(i, j, k))
                T_work(i, j, k) = single(t_candidate);
            end
        end
    end
end

T_out = T_work;
end

function idx_range = sweep_range(n, sign_dir)
if sign_dir >= 0
    idx_range = 1:n;
else
    idx_range = n:-1:1;
end
end

function cand = solve_two_axis(v1, v2, inv_h12, inv_h22, slowness2)
A = inv_h12 + inv_h22;
B = v1 * inv_h12 + v2 * inv_h22;
C = v1 * v1 * inv_h12 + v2 * v2 * inv_h22 - slowness2;
disc = B * B - A * C;
if disc < 0
    cand = inf;
else
    cand = (B + sqrt(disc)) / A;
end
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
    if all(curr == start_idx(:))
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
if fx < 1 || fy < 1 || fz < 1 || fx > sz(1) || fy > sz(2) || fz > sz(3)
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
