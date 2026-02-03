function [hit, t, u, v, intersection_point] = moller_trumbore(ray_origin, ray_dir, v0, v1, v2, epsilon)
%MOLLER_TRUMBORE Ray-triangle intersection using Möller-Trumbore algorithm
%
% Computes the intersection of a ray with a triangle using the efficient
% Möller-Trumbore algorithm. This is the industry-standard method for
% ray-triangle intersection in computer graphics and simulation.
%
% Reference:
%   Möller, T., & Trumbore, B. (1997). "Fast, minimum storage ray-triangle
%   intersection." Journal of Graphics Tools, 2(1), 21-28.
%
% Inputs:
%   ray_origin - [3x1] Ray starting point [x; y; z]
%   ray_dir    - [3x1] Ray direction vector (should be normalized for t to be distance)
%   v0, v1, v2 - [3x1] Triangle vertex positions in counter-clockwise order
%   epsilon    - (optional) Tolerance for parallel detection (default: 1e-8)
%
% Outputs:
%   hit        - true if ray intersects triangle
%   t          - Distance along ray to intersection (ray_origin + t*ray_dir)
%   u, v       - Barycentric coordinates of intersection point
%                (point = (1-u-v)*v0 + u*v1 + v*v2)
%   intersection_point - [3x1] The actual intersection point coordinates
%
% Barycentric coordinates (u, v, w) where w = 1 - u - v:
%   - w = weight for v0
%   - u = weight for v1
%   - v = weight for v2
%   - Valid intersection: u >= 0, v >= 0, u + v <= 1
%
% Example:
%   % Define a triangle
%   v0 = [0; 0; 0];
%   v1 = [1; 0; 0];
%   v2 = [0; 1; 0];
%
%   % Ray from above, pointing down
%   origin = [0.2; 0.2; 1];
%   direction = [0; 0; -1];
%
%   [hit, t, u, v, point] = moller_trumbore(origin, direction, v0, v1, v2);
%   % hit = true, t = 1, point = [0.2; 0.2; 0]
%
% See also: terrain_mesh, radar_ray_caster

% Default epsilon
if nargin < 6
    epsilon = 1e-8;
end

% Ensure column vectors
ray_origin = ray_origin(:);
ray_dir = ray_dir(:);
v0 = v0(:);
v1 = v1(:);
v2 = v2(:);

% Initialize outputs
hit = false;
t = inf;
u = 0;
v = 0;
intersection_point = nan(3, 1);

% Edge vectors
e1 = v1 - v0;
e2 = v2 - v0;

% Calculate determinant (cross product of direction and e2)
p = cross(ray_dir, e2);
det = dot(e1, p);

% If determinant is near zero, ray lies in plane of triangle or is parallel
if abs(det) < epsilon
    return;
end

inv_det = 1.0 / det;

% Calculate distance from v0 to ray origin
tvec = ray_origin - v0;

% Calculate u parameter and test bounds
u = dot(tvec, p) * inv_det;
if u < 0.0 || u > 1.0
    return;
end

% Calculate v parameter and test bounds
q = cross(tvec, e1);
v = dot(ray_dir, q) * inv_det;
if v < 0.0 || (u + v) > 1.0
    return;
end

% Calculate t (distance along ray)
t = dot(e2, q) * inv_det;

% Check if intersection is in positive ray direction
if t > epsilon
    hit = true;
    intersection_point = ray_origin + t * ray_dir;
else
    t = inf;
end

end
