function smoothed_path = smooth_path_geometric(path, terrain, iterations)
%SMOOTH_PATH_GEOMETRIC Smooths a path using iterative corner cutting
%
% Uses Chaikin's algorithm or similar to cut corners, resulting in a
% denser but smoother path that is easier for the trajectory optimizer
% to track.
%
% Inputs:
%   path       - [3 x N] waypoints
%   terrain    - terrain_map object (for safety checks)
%   iterations - Number of smoothing iterations (default: 2)
%
% Outputs:
%   smoothed_path - [3 x M] smoothed waypoints
%
% Author: Quadrotor Terrain Following Project

if nargin < 3
    iterations = 2;
end

current_path = path;
min_clearance = 5.0; % Hard-coded safety margin in smoothing

for k = 1:iterations
    n_pts = size(current_path, 2);
    if n_pts < 3
        break;
    end

    new_path = zeros(3, 2*n_pts - 2);

    % Keep start point
    new_path(:, 1) = current_path(:, 1);

    idx = 2;
    for i = 1:n_pts-1
        p1 = current_path(:, i);
        p2 = current_path(:, i+1);

        % Chaikin's: Q = 0.75*P1 + 0.25*P2
        %            R = 0.25*P1 + 0.75*P2

        q = 0.75 * p1 + 0.25 * p2;
        r = 0.25 * p1 + 0.75 * p2;

        % Safety Check: Ensure new points are above terrain
        h_q = terrain.get_height(q(1), q(2));
        h_r = terrain.get_height(r(1), r(2));

        if -q(3) < h_q + min_clearance
            q(3) = -(h_q + min_clearance);
        end

        if -r(3) < h_r + min_clearance
            r(3) = -(h_r + min_clearance);
        end

        if i == 1
            % First segment, keep p1 (already added)
            new_path(:, idx) = r;
            idx = idx + 1;
        elseif i == n_pts - 1
            % Last segment
            new_path(:, idx) = q;
            idx = idx + 1;
        else
            new_path(:, idx) = q;
            new_path(:, idx+1) = r;
            idx = idx + 2;
        end
    end

    % Keep end point
    new_path(:, idx) = current_path(:, end);

    % Resize (trim unused)
    current_path = new_path(:, 1:idx);
end

smoothed_path = current_path;

end
