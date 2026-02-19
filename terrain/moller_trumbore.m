function [hit, t, u, v, intersection_point] = moller_trumbore(ray_origin, ray_dir, v0, v1, v2, epsilon)
%MOLLER_TRUMBORE Vectorized Moller-Trumbore ray-triangle intersection.
%
% Batch mode:
% - ray_origin, ray_dir: [3 x N]
% - v0, v1, v2: [3 x 1] or [3 x N]
%   * [3 x 1] triangle inputs are broadcast to all rays.
%   * [3 x N] triangle inputs are paired column-wise with rays.
%
% Outputs:
% - hit: [1 x N] logical
% - t, u, v: [1 x N]
% - intersection_point: [3 x N], NaN where hit is false
%
% Notes:
% - No scalar per-ray loops or scalar if-branches are used.
% - Validation is handled with logical masks.

if nargin < 6
    epsilon = 1e-8;
end

if size(ray_origin, 1) ~= 3 || size(ray_dir, 1) ~= 3 || ...
   size(v0, 1) ~= 3 || size(v1, 1) ~= 3 || size(v2, 1) ~= 3
    error('moller_trumbore:ShapeError', ...
        'Inputs must have 3 rows: ray_origin, ray_dir, v0, v1, v2 are [3 x N].');
end

n_cols = max([size(ray_origin, 2), size(ray_dir, 2), size(v0, 2), size(v1, 2), size(v2, 2)]);

ray_origin = expand_to_n(ray_origin, n_cols, 'ray_origin');
ray_dir = expand_to_n(ray_dir, n_cols, 'ray_dir');
v0 = expand_to_n(v0, n_cols, 'v0');
v1 = expand_to_n(v1, n_cols, 'v1');
v2 = expand_to_n(v2, n_cols, 'v2');

e1 = v1 - v0;
e2 = v2 - v0;

% p = cross(ray_dir, e2, 1)
p = cross(ray_dir, e2, 1);
det = sum(e1 .* p, 1);

det_mask = abs(det) > epsilon;

inv_det = zeros(1, n_cols, 'like', det);
inv_det(det_mask) = 1.0 ./ det(det_mask);

tvec = ray_origin - v0;
u = sum(tvec .* p, 1) .* inv_det;

q = cross(tvec, e1, 1);
v = sum(ray_dir .* q, 1) .* inv_det;
t = sum(e2 .* q, 1) .* inv_det;

valid_mask = det_mask & (u >= 0.0) & (u <= 1.0) & ...
             (v >= 0.0) & ((u + v) <= 1.0) & (t > epsilon);

hit = valid_mask;

t(~valid_mask) = inf;
u(~valid_mask) = 0;
v(~valid_mask) = 0;

intersection_point = nan(3, n_cols);
if any(valid_mask)
    intersection_point(:, valid_mask) = ray_origin(:, valid_mask) + ...
        ray_dir(:, valid_mask) .* t(valid_mask);
end

end

function x = expand_to_n(x, n_cols, name)
% Expand [3x1] -> [3xN], or validate [3xN].

if size(x, 2) == 1 && n_cols > 1
    x = x(:, ones(1, n_cols));
elseif size(x, 2) ~= n_cols
    error('moller_trumbore:BatchSizeMismatch', ...
        '%s has %d columns; expected 1 or %d.', name, size(x, 2), n_cols);
end

end
