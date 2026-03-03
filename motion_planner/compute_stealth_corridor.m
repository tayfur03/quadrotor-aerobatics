function corridor = compute_stealth_corridor(V, start_coord, goal_coord, visibility_threshold, alpha, params)
%COMPUTE_STEALTH_CORRIDOR Build a deterministic safe sampling tube in 3D.
%
% Inputs
%   V                    [n1 x n2 x n3] visibility/risk voxel grid
%   start_coord          [3x1] start coordinate
%   goal_coord           [3x1] goal coordinate
%   visibility_threshold scalar threshold for "safe" voxels
%   alpha                scalar cost sensitivity
%   params (optional) fields:
%     .N_vec, .E_vec, .alt_vec  : world axes (if provided, V must be [nE x nN x nAlt])
%     .max_sweeps               : FSM sweep count (default: 56)
%     .tol                      : convergence tolerance (default: 1e-3)
%     .backtrack_step_vox       : backtracking step in voxels (default: 0.8)
%     .max_backtrack_points     : max backtracking points (default: 6000)
%     .r_min                    : min corridor radius (m if axes provided, else vox)
%     .r_max                    : max corridor radius (m if axes provided, else vox)
%     .lambda                   : distance-transform scale factor (default: 1.8)
%     .blocked_penalty          : large additive RHS on unsafe voxels (default: 1e3)
%
% Output
%   corridor struct:
%     .skeleton        [3xM] [N;E;Alt] if axes provided, otherwise voxel coords
%     .skeleton_idx    [3xM] voxel coordinates [i;j;k]
%     .radii           [1xM] radius values (m or voxels)
%     .radii_voxel     [1xM] radius values in voxels
%     .T               arrival-time map
%     .V               input visibility map
%     .safe_mask       safe-voxel mask
%     .bounds_world    [Nmin Nmax Emin Emax Altmin Altmax] (or voxel bounds)
%     .grid            axis metadata
%
% Safety/certification rationale:
% A deterministic geometric tube is auditable and repeatable. Unlike a
% learned GAN initializer, this produces explicit geometry (skeleton+radii)
% that can be independently verified in certification-oriented workflows.

if nargin < 6
    params = struct();
end

validateattributes(V, {'numeric', 'logical'}, {'nonempty', 'real', '3d'}, mfilename, 'V', 1);
validateattributes(start_coord, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'start_coord', 2);
validateattributes(goal_coord, {'numeric'}, {'real', 'vector', 'numel', 3}, mfilename, 'goal_coord', 3);
validateattributes(visibility_threshold, {'numeric'}, {'real', 'scalar'}, mfilename, 'visibility_threshold', 4);
validateattributes(alpha, {'numeric'}, {'real', 'scalar', 'finite', 'nonnegative'}, mfilename, 'alpha', 5);

V = double(V);
sz = size(V);
n1 = sz(1);
n2 = sz(2);
n3 = sz(3);

has_axes = isfield(params, 'N_vec') && isfield(params, 'E_vec') && isfield(params, 'alt_vec') && ...
           ~isempty(params.N_vec) && ~isempty(params.E_vec) && ~isempty(params.alt_vec);

if has_axes
    N_vec = params.N_vec(:)';
    E_vec = params.E_vec(:)';
    alt_vec = params.alt_vec(:)';

    if n1 ~= numel(E_vec) || n2 ~= numel(N_vec) || n3 ~= numel(alt_vec)
        error(['compute_stealth_corridor:SizeMismatch ', ...
               'With world axes, V must be [numel(E_vec) x numel(N_vec) x numel(alt_vec)].']);
    end

    start_idx = world_to_voxel(start_coord(:), N_vec, E_vec, alt_vec);
    goal_idx = world_to_voxel(goal_coord(:), N_vec, E_vec, alt_vec);

    dE = median(diff(E_vec));
    dN = median(diff(N_vec));
    dA = median(diff(alt_vec));
    voxel_scale = mean([abs(dE), abs(dN), abs(dA)]);
else
    N_vec = [];
    E_vec = [];
    alt_vec = [];
    start_idx = start_coord(:);
    goal_idx = goal_coord(:);
    voxel_scale = 1.0;
end

start_idx = clamp_to_grid(start_idx, sz);
goal_idx = clamp_to_grid(goal_idx, sz);
start_cell = round(start_idx);
goal_cell = round(goal_idx);

safe_mask = V <= visibility_threshold;
safe_mask(start_cell(1), start_cell(2), start_cell(3)) = true;
safe_mask(goal_cell(1), goal_cell(2), goal_cell(3)) = true;

max_sweeps = get_opt(params, 'max_sweeps', 56);
tol = get_opt(params, 'tol', 1e-3);
blocked_penalty = get_opt(params, 'blocked_penalty', 1e3);
h = max(get_opt(params, 'voxel_step', 1.0), 1e-6);

default_r_min = 1.5 * voxel_scale;
default_r_max = 12.0 * voxel_scale;
r_min = get_opt(params, 'r_min', default_r_min);
r_max = get_opt(params, 'r_max', default_r_max);
lambda = get_opt(params, 'lambda', 1.8);

backtrack_step_vox = get_opt(params, 'backtrack_step_vox', 0.8);
max_backtrack_points = round(get_opt(params, 'max_backtrack_points', 6000));
max_backtrack_points = max(32, max_backtrack_points);

rhs = exp(-alpha .* V);
rhs(~safe_mask) = rhs(~safe_mask) + blocked_penalty;
rhs = max(rhs, 1e-9);

T = inf(sz);
T(start_cell(1), start_cell(2), start_cell(3)) = 0.0;

sweep_dirs = [...
     1,  1,  1; ...
    -1,  1,  1; ...
     1, -1,  1; ...
    -1, -1,  1; ...
     1,  1, -1; ...
    -1,  1, -1; ...
     1, -1, -1; ...
    -1, -1, -1];

converged = false;
last_delta = inf;
sweeps_done = 0;

for s = 1:max_sweeps
    T_prev = T;
    for d = 1:size(sweep_dirs, 1)
        T = directional_update(T, rhs, safe_mask, sweep_dirs(d, :), h);
        T(start_cell(1), start_cell(2), start_cell(3)) = 0.0;
    end

    dif = abs(T(safe_mask) - T_prev(safe_mask));
    if isempty(dif)
        last_delta = 0.0;
    else
        last_delta = max(dif);
    end
    sweeps_done = s;
    if last_delta < tol
        converged = true;
        break;
    end
end

[dTi, dTj, dTk] = gradient(T);

path_idx = nan(3, max_backtrack_points);
path_idx(:, 1) = goal_idx;
path_count = 1;
p = goal_idx;
reached_start = false;

for it = 2:max_backtrack_points
    if norm(p - start_idx) <= 1.25
        path_idx(:, it) = start_idx;
        path_count = it;
        reached_start = true;
        break;
    end

    g = sample_gradient(dTi, dTj, dTk, p);
    g_norm = norm(g);
    curr_T = sample_field(T, p);

    if ~isfinite(g_norm) || g_norm < 1e-10 || ~isfinite(curr_T)
        p_next = descend_neighbor(T, safe_mask, p);
    else
        p_step = p - backtrack_step_vox * (g / g_norm);
        p_step = clamp_to_grid(p_step, sz);
        T_step = sample_field(T, p_step);
        if ~isfinite(T_step) || (T_step > curr_T - 1e-8)
            p_next = descend_neighbor(T, safe_mask, p);
        else
            p_next = p_step;
        end
    end

    if norm(p_next - p) < 1e-6
        p_next = descend_neighbor(T, safe_mask, p);
        if norm(p_next - p) < 1e-6
            path_count = it - 1;
            break;
        end
    end

    p = p_next;
    path_idx(:, it) = p;
    path_count = it;
end

if ~reached_start && path_count < max_backtrack_points
    path_count = path_count + 1;
    path_idx(:, path_count) = start_idx;
end

path_idx = path_idx(:, 1:path_count);
path_idx = fliplr(path_idx);

if isempty(which('bwdist'))
    error('compute_stealth_corridor:MissingToolbox', ...
          'bwdist is required (Image Processing Toolbox).');
end

D = bwdist(~safe_mask);
D_path = interpn(D, path_idx(1, :), path_idx(2, :), path_idx(3, :), 'linear', 0);

r_min_vox = max(r_min / voxel_scale, 0.5);
r_max_vox = max(r_max / voxel_scale, r_min_vox);
radii_voxel = min(max(lambda .* D_path, r_min_vox), r_max_vox);
radii = radii_voxel .* voxel_scale;

if has_axes
    skeleton_world = voxel_to_world(path_idx, N_vec, E_vec, alt_vec);
else
    skeleton_world = path_idx;
end

bounds_world = corridor_bounds(skeleton_world, radii);

grid_info = struct();
grid_info.has_world_axes = has_axes;
grid_info.N_vec = N_vec;
grid_info.E_vec = E_vec;
grid_info.alt_vec = alt_vec;
grid_info.size = sz;
grid_info.voxel_scale = voxel_scale;
grid_info.start_idx = start_idx;
grid_info.goal_idx = goal_idx;

meta = struct();
meta.max_sweeps = max_sweeps;
meta.sweeps_done = sweeps_done;
meta.converged = converged;
meta.final_delta = last_delta;

corridor = struct();
corridor.skeleton = skeleton_world;
corridor.skeleton_idx = path_idx;
corridor.radii = radii;
corridor.radii_voxel = radii_voxel;
corridor.T = T;
corridor.V = V;
corridor.safe_mask = safe_mask;
corridor.visibility_threshold = visibility_threshold;
corridor.alpha = alpha;
corridor.bounds_world = bounds_world;
corridor.grid = grid_info;
corridor.meta = meta;
end

function T_out = directional_update(T_in, rhs, safe_mask, dir_sign, h)
T_work = T_in;
rhs_work = rhs;
safe_work = safe_mask;

if dir_sign(1) < 0
    T_work = flip(T_work, 1);
    rhs_work = flip(rhs_work, 1);
    safe_work = flip(safe_work, 1);
end
if dir_sign(2) < 0
    T_work = flip(T_work, 2);
    rhs_work = flip(rhs_work, 2);
    safe_work = flip(safe_work, 2);
end
if dir_sign(3) < 0
    T_work = flip(T_work, 3);
    rhs_work = flip(rhs_work, 3);
    safe_work = flip(safe_work, 3);
end

left = cat(1, T_work(1, :, :), T_work(1:end-1, :, :));
right = cat(1, T_work(2:end, :, :), T_work(end, :, :));
a = min(left, right);

down = cat(2, T_work(:, 1, :), T_work(:, 1:end-1, :));
up = cat(2, T_work(:, 2:end, :), T_work(:, end, :));
b = min(down, up);

back = cat(3, T_work(:, :, 1), T_work(:, :, 1:end-1));
front = cat(3, T_work(:, :, 2:end), T_work(:, :, end));
c = min(back, front);

T_cand = solve_eikonal_vec(a, b, c, rhs_work, h);
T_work(safe_work) = min(T_work(safe_work), T_cand(safe_work));

if dir_sign(1) < 0
    T_work = flip(T_work, 1);
end
if dir_sign(2) < 0
    T_work = flip(T_work, 2);
end
if dir_sign(3) < 0
    T_work = flip(T_work, 3);
end

T_out = T_work;
end

function u = solve_eikonal_vec(a, b, c, f, h)
abc = cat(4, a, b, c);
abc = sort(abc, 4, 'ascend');
a1 = abc(:, :, :, 1);
b1 = abc(:, :, :, 2);
c1 = abc(:, :, :, 3);

fh = f .* h;
u = a1 + fh;

mask2 = u > b1;
disc2 = max(0.0, 2.0 .* (fh .^ 2) - (a1 - b1) .^ 2);
u2 = 0.5 .* (a1 + b1 + sqrt(disc2));
u(mask2) = u2(mask2);

mask3 = u > c1;
sabc = a1 + b1 + c1;
disc3 = max(0.0, sabc .^ 2 - 3.0 .* (a1 .^ 2 + b1 .^ 2 + c1 .^ 2 - fh .^ 2));
u3 = (sabc + sqrt(disc3)) ./ 3.0;
u(mask3) = u3(mask3);

u = max(u, a1);
u(~isfinite(u)) = inf;
end

function idx = world_to_voxel(p_world, N_vec, E_vec, alt_vec)
idx = zeros(3, 1);
idx(1) = interp1(E_vec, 1:numel(E_vec), p_world(2), 'linear', 'extrap');
idx(2) = interp1(N_vec, 1:numel(N_vec), p_world(1), 'linear', 'extrap');
idx(3) = interp1(alt_vec, 1:numel(alt_vec), p_world(3), 'linear', 'extrap');
end

function p_world = voxel_to_world(p_idx, N_vec, E_vec, alt_vec)
p_world = zeros(3, size(p_idx, 2));
p_world(1, :) = interp1(1:numel(N_vec), N_vec, p_idx(2, :), 'linear', 'extrap');
p_world(2, :) = interp1(1:numel(E_vec), E_vec, p_idx(1, :), 'linear', 'extrap');
p_world(3, :) = interp1(1:numel(alt_vec), alt_vec, p_idx(3, :), 'linear', 'extrap');
end

function p = clamp_to_grid(p, sz)
p = p(:);
p(1) = min(max(p(1), 1), sz(1));
p(2) = min(max(p(2), 1), sz(2));
p(3) = min(max(p(3), 1), sz(3));
end

function g = sample_gradient(dTi, dTj, dTk, p)
g = zeros(3, 1);
g(1) = interpn(dTi, p(1), p(2), p(3), 'linear', 0);
g(2) = interpn(dTj, p(1), p(2), p(3), 'linear', 0);
g(3) = interpn(dTk, p(1), p(2), p(3), 'linear', 0);
end

function v = sample_field(F, p)
v = interpn(F, p(1), p(2), p(3), 'linear', inf);
end

function p_next = descend_neighbor(T, safe_mask, p_curr)
sz = size(T);
p0 = round(p_curr(:));

i0 = p0(1);
j0 = p0(2);
k0 = p0(3);

i_min = max(1, i0 - 1);
i_max = min(sz(1), i0 + 1);
j_min = max(1, j0 - 1);
j_max = min(sz(2), j0 + 1);
k_min = max(1, k0 - 1);
k_max = min(sz(3), k0 + 1);

[I, J, K] = ndgrid(i_min:i_max, j_min:j_max, k_min:k_max);
I = I(:);
J = J(:);
K = K(:);

is_self = (I == i0) & (J == j0) & (K == k0);
I = I(~is_self);
J = J(~is_self);
K = K(~is_self);

lin = sub2ind(sz, I, J, K);
vals = T(lin);
vals(~safe_mask(lin)) = inf;

[best_val, best_idx] = min(vals);
curr_val = T(i0, j0, k0);

if isfinite(best_val) && best_val < curr_val
    p_next = [I(best_idx); J(best_idx); K(best_idx)];
else
    p_next = p_curr(:);
end
end

function b = corridor_bounds(centers, radii)
N_min = min(centers(1, :) - radii);
N_max = max(centers(1, :) + radii);
E_min = min(centers(2, :) - radii);
E_max = max(centers(2, :) + radii);
A_min = min(centers(3, :) - radii);
A_max = max(centers(3, :) + radii);
b = [N_min, N_max, E_min, E_max, A_min, A_max];
end

function val = get_opt(params, name, default)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
