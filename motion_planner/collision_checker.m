function [is_collision, min_dist, nearest_obs_idx] = collision_checker(position, obstacles, safety_margin)
% COLLISION_CHECKER Check if a position collides with any obstacle
%
% Inputs:
%   position      - [3x1] Current position in NED frame [N; E; D]
%   obstacles     - Cell array of obstacle structs, each with:
%                   - type: 'sphere' or 'box'
%                   - For sphere: center [3x1], radius (scalar)
%                   - For box: center [3x1], half_size [3x1]
%   safety_margin - (optional) Additional safety distance, default 0.5m
%
% Outputs:
%   is_collision    - true if position is inside any obstacle (+ margin)
%   min_dist        - Minimum distance to nearest obstacle surface
%   nearest_obs_idx - Index of nearest obstacle (0 if no obstacles)
%
% Example:
%   obstacles = {struct('type', 'sphere', 'center', [5;2;-2], 'radius', 1.0)};
%   [collision, dist, idx] = collision_checker([5;2;-2], obstacles, 0.5);

    if nargin < 3
        safety_margin = 0.5;  % Default 0.5m safety margin
    end

    % Initialize outputs
    is_collision = false;
    min_dist = inf;
    nearest_obs_idx = 0;

    if isempty(obstacles)
        return;
    end

    n_obs = length(obstacles);

    for i = 1:n_obs
        obs = obstacles{i};

        switch obs.type
            case 'sphere'
                dist = distance_to_sphere(position, obs.center, obs.radius);
            case 'box'
                dist = distance_to_box(position, obs.center, obs.half_size);
            otherwise
                warning('Unknown obstacle type: %s', obs.type);
                continue;
        end

        if dist < min_dist
            min_dist = dist;
            nearest_obs_idx = i;
        end

        if dist < safety_margin
            is_collision = true;
        end
    end
end

function dist = distance_to_sphere(point, center, radius)
% Distance from point to sphere surface (negative if inside)
    dist = norm(point - center) - radius;
end

function dist = distance_to_box(point, center, half_size)
% Distance from point to axis-aligned box surface (negative if inside)
% half_size is [hx; hy; hz] - half the box dimensions

    % Vector from center to point
    d = abs(point - center);

    % Clamp to box surface
    clamped = max(d - half_size, 0);

    % Distance to surface
    dist_outside = norm(clamped);

    if dist_outside > 0
        % Point is outside the box
        dist = dist_outside;
    else
        % Point is inside the box - return negative distance
        % (distance to nearest face, negated)
        dist_to_faces = half_size - d;
        dist = -min(dist_to_faces);
    end
end
