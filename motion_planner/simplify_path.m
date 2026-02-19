function simplified = simplify_path(path, tolerance)
%SIMPLIFY_PATH Reduce waypoints using Douglas-Peucker algorithm
%
% Removes intermediate waypoints that are within 'tolerance' distance
% from the line connecting their neighbors. This reduces zigzag while
% preserving the overall path shape.
%
% Inputs:
%   path      - [3 x N] waypoints
%   tolerance - Maximum perpendicular distance to keep a point
%
% Output:
%   simplified - [3 x M] simplified path (M <= N)

    if size(path, 2) <= 2
        simplified = path;
        return;
    end
    
    % Douglas-Peucker recursive simplification
    keep = douglas_peucker(path, 1, size(path, 2), tolerance);
    
    % Always keep first and last
    keep(1) = true;
    keep(end) = true;
    
    simplified = path(:, keep);
end

function keep = douglas_peucker(path, start_idx, end_idx, tolerance)
    n = size(path, 2);
    keep = false(1, n);
    
    if end_idx - start_idx < 2
        keep(start_idx) = true;
        keep(end_idx) = true;
        return;
    end
    
    % Find point with maximum distance from line
    p1 = path(:, start_idx);
    p2 = path(:, end_idx);
    line_vec = p2 - p1;
    line_len = norm(line_vec);
    
    if line_len < 1e-6
        keep(start_idx) = true;
        keep(end_idx) = true;
        return;
    end
    
    line_unit = line_vec / line_len;
    
    max_dist = 0;
    max_idx = start_idx;
    
    for i = (start_idx + 1):(end_idx - 1)
        pt = path(:, i);
        vec_to_pt = pt - p1;
        proj_len = dot(vec_to_pt, line_unit);
        proj_pt = p1 + proj_len * line_unit;
        dist = norm(pt - proj_pt);
        
        if dist > max_dist
            max_dist = dist;
            max_idx = i;
        end
    end
    
    if max_dist > tolerance
        % Recursively simplify both halves
        keep1 = douglas_peucker(path, start_idx, max_idx, tolerance);
        keep2 = douglas_peucker(path, max_idx, end_idx, tolerance);
        keep = keep1 | keep2;
    else
        % Keep only endpoints of this segment
        keep(start_idx) = true;
        keep(end_idx) = true;
    end
end
