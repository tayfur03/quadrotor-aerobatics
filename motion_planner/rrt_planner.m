function [path, success, tree] = rrt_planner(start, goal, obstacles, bounds, params)
% RRT_PLANNER Rapidly-exploring Random Tree path planner
%
% Implements the RRT algorithm for 3D path planning with obstacle avoidance.
%
% Inputs:
%   start     - [3x1] Start position in NED frame
%   goal      - [3x1] Goal position in NED frame
%   obstacles - Cell array of obstacle structs (from obstacle_manager)
%   bounds    - Struct with fields: min [3x1], max [3x1] - workspace bounds
%   params    - (optional) Struct with RRT parameters:
%               - max_iter: Maximum iterations (default: 5000)
%               - step_size: Extension step size in meters (default: 0.5)
%               - goal_bias: Probability of sampling goal (default: 0.1)
%               - goal_tol: Distance to consider goal reached (default: 0.5)
%               - safety_margin: Collision check margin (default: 0.3)
%
% Outputs:
%   path    - [3xM] Matrix of waypoints from start to goal
%   success - true if path found, false otherwise
%   tree    - Struct with tree data for visualization
%
% Example:
%   start = [0; 0; -2];
%   goal = [10; 5; -3];
%   obstacles = {struct('type', 'sphere', 'center', [5;2;-2.5], 'radius', 1.5)};
%   bounds = struct('min', [-5; -5; -10], 'max', [15; 10; 0]);
%   [path, success] = rrt_planner(start, goal, obstacles, bounds);

    % Default parameters
    if nargin < 5
        params = struct();
    end

    max_iter = get_param(params, 'max_iter', 5000);
    step_size = get_param(params, 'step_size', 0.5);
    goal_bias = get_param(params, 'goal_bias', 0.1);
    goal_tol = get_param(params, 'goal_tol', 0.5);
    safety_margin = get_param(params, 'safety_margin', 0.3);

    % Ensure column vectors
    start = start(:);
    goal = goal(:);

    % Initialize tree with start node
    % Tree structure: nodes (positions), parents (indices)
    tree.nodes = start;
    tree.parents = 0;  % Root has no parent

    success = false;
    goal_node_idx = -1;

    for iter = 1:max_iter
        % Sample random point (with goal bias)
        if rand() < goal_bias
            q_rand = goal;
        else
            q_rand = sample_random_point(bounds);
        end

        % Find nearest node in tree
        [nearest_idx, nearest_node] = find_nearest(tree.nodes, q_rand);

        % Extend tree toward random point
        q_new = extend(nearest_node, q_rand, step_size);

        % Check if new node is collision-free
        if is_path_collision_free(nearest_node, q_new, obstacles, safety_margin)
            % Add new node to tree
            tree.nodes = [tree.nodes, q_new];
            tree.parents = [tree.parents, nearest_idx];

            new_node_idx = size(tree.nodes, 2);

            % Check if goal reached
            if norm(q_new - goal) < goal_tol
                success = true;
                goal_node_idx = new_node_idx;
                break;
            end
        end
    end

    % Extract path by backtracking from goal
    if success
        path = extract_path(tree, goal_node_idx);
        % Add exact goal position
        path = [path, goal];
        fprintf('RRT: Path found in %d iterations with %d waypoints\n', iter, size(path, 2));
    else
        path = start;  % Return just start if no path found
        fprintf('RRT: No path found after %d iterations\n', max_iter);
    end
end

function val = get_param(params, name, default)
% Get parameter value or default
    if isfield(params, name)
        val = params.(name);
    else
        val = default;
    end
end

function q = sample_random_point(bounds)
% Sample uniformly random point within bounds
    range = bounds.max - bounds.min;
    q = bounds.min + range .* rand(3, 1);
end

function [nearest_idx, nearest_node] = find_nearest(nodes, q)
% Find nearest node in tree to query point
    n_nodes = size(nodes, 2);
    min_dist = inf;
    nearest_idx = 1;

    for i = 1:n_nodes
        d = norm(nodes(:, i) - q);
        if d < min_dist
            min_dist = d;
            nearest_idx = i;
        end
    end

    nearest_node = nodes(:, nearest_idx);
end

function q_new = extend(q_near, q_rand, step_size)
% Extend from q_near toward q_rand by step_size
    direction = q_rand - q_near;
    dist = norm(direction);

    if dist < step_size
        q_new = q_rand;
    else
        q_new = q_near + (step_size / dist) * direction;
    end
end

function collision_free = is_path_collision_free(q1, q2, obstacles, safety_margin)
% Check if straight line path from q1 to q2 is collision-free
%
% Uses line sampling to check for collisions along the path

    if isempty(obstacles)
        collision_free = true;
        return;
    end

    % Sample points along the line
    dist = norm(q2 - q1);
    n_samples = max(2, ceil(dist / 0.1));  % Sample every 0.1m

    for i = 0:n_samples
        alpha = i / n_samples;
        q_sample = q1 + alpha * (q2 - q1);

        [is_collision, ~, ~] = collision_checker(q_sample, obstacles, safety_margin);

        if is_collision
            collision_free = false;
            return;
        end
    end

    collision_free = true;
end

function path = extract_path(tree, goal_idx)
% Extract path from tree by backtracking from goal to root
    path = [];
    current_idx = goal_idx;

    while current_idx > 0
        path = [tree.nodes(:, current_idx), path]; %#ok<AGROW>
        current_idx = tree.parents(current_idx);
    end
end
