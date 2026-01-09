function [path, success, tree] = rrt_planner(start, goal, obstacles, bounds, params)
% RRT_PLANNER Rapidly-exploring Random Tree path planner (RRT and RRT*)
%
% Implements both RRT and RRT* algorithms for 3D path planning with
% obstacle avoidance. RRT* provides asymptotically optimal paths by
% rewiring the tree to minimize path cost.
%
% Inputs:
%   start     - [3x1] Start position in NED frame
%   goal      - [3x1] Goal position in NED frame
%   obstacles - Cell array of obstacle structs (from obstacle_manager)
%   bounds    - Struct with fields: min [3x1], max [3x1] - workspace bounds
%   params    - (optional) Struct with RRT parameters:
%               - algorithm: 'rrt' or 'rrt_star' (default: 'rrt')
%               - max_iter: Maximum iterations (default: 5000)
%               - step_size: Extension step size in meters (default: 0.5)
%               - goal_bias: Probability of sampling goal (default: 0.1)
%               - goal_tol: Distance to consider goal reached (default: 0.5)
%               - safety_margin: Collision check margin (default: 0.3)
%               - rewire_radius: Radius for RRT* rewiring (default: 1.5)
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
%   params.algorithm = 'rrt_star';
%   [path, success] = rrt_planner(start, goal, obstacles, bounds, params);

    % Default parameters
    if nargin < 5
        params = struct();
    end

    algorithm = get_param(params, 'algorithm', 'rrt');
    max_iter = get_param(params, 'max_iter', 5000);
    step_size = get_param(params, 'step_size', 0.5);
    goal_bias = get_param(params, 'goal_bias', 0.1);
    goal_tol = get_param(params, 'goal_tol', 0.5);
    safety_margin = get_param(params, 'safety_margin', 0.3);
    rewire_radius = get_param(params, 'rewire_radius', 1.5);

    use_rrt_star = strcmpi(algorithm, 'rrt_star') || strcmpi(algorithm, 'rrt*');

    % Ensure column vectors
    start = start(:);
    goal = goal(:);

    % Initialize tree with start node
    % Tree structure: nodes (positions), parents (indices), costs (path cost from start)
    tree.nodes = start;
    tree.parents = 0;  % Root has no parent
    tree.costs = 0;    % Cost to reach start is 0

    success = false;
    goal_node_idx = -1;
    best_goal_cost = inf;

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

            if use_rrt_star
                % RRT*: Find best parent among nearby nodes
                [tree, new_node_idx] = rrt_star_extend(tree, q_new, nearest_idx, ...
                    obstacles, safety_margin, rewire_radius);
            else
                % Standard RRT: Just add the node
                tree.nodes = [tree.nodes, q_new];
                tree.parents = [tree.parents, nearest_idx];
                new_cost = tree.costs(nearest_idx) + norm(q_new - nearest_node);
                tree.costs = [tree.costs, new_cost];
                new_node_idx = size(tree.nodes, 2);
            end

            % Check if goal reached
            dist_to_goal = norm(q_new - goal);
            if dist_to_goal < goal_tol
                goal_cost = tree.costs(new_node_idx) + dist_to_goal;

                if ~use_rrt_star
                    % Standard RRT: Stop at first solution
                    success = true;
                    goal_node_idx = new_node_idx;
                    break;
                else
                    % RRT*: Keep improving until max_iter
                    if goal_cost < best_goal_cost
                        success = true;
                        goal_node_idx = new_node_idx;
                        best_goal_cost = goal_cost;
                    end
                end
            end
        end
    end

    % Extract path by backtracking from goal
    if success
        path = extract_path(tree, goal_node_idx);
        % Add exact goal position
        path = [path, goal];

        if use_rrt_star
            fprintf('RRT*: Path found with cost %.2f in %d iterations (%d waypoints)\n', ...
                best_goal_cost, iter, size(path, 2));
        else
            fprintf('RRT: Path found in %d iterations with %d waypoints\n', iter, size(path, 2));
        end
    else
        path = start;  % Return just start if no path found
        fprintf('%s: No path found after %d iterations\n', upper(algorithm), max_iter);
    end
end

function [tree, new_node_idx] = rrt_star_extend(tree, q_new, nearest_idx, obstacles, safety_margin, rewire_radius)
% RRT* extension with optimal parent selection and tree rewiring

    n_nodes = size(tree.nodes, 2);

    % Find all nodes within rewire_radius
    near_indices = [];
    for i = 1:n_nodes
        if norm(tree.nodes(:, i) - q_new) < rewire_radius
            near_indices = [near_indices, i]; %#ok<AGROW>
        end
    end

    % Find best parent (minimum cost to reach q_new)
    best_parent = nearest_idx;
    best_cost = tree.costs(nearest_idx) + norm(q_new - tree.nodes(:, nearest_idx));

    for i = 1:length(near_indices)
        idx = near_indices(i);
        potential_cost = tree.costs(idx) + norm(q_new - tree.nodes(:, idx));

        if potential_cost < best_cost
            % Check if path to this node is collision-free
            if is_path_collision_free(tree.nodes(:, idx), q_new, obstacles, safety_margin)
                best_parent = idx;
                best_cost = potential_cost;
            end
        end
    end

    % Add new node with best parent
    tree.nodes = [tree.nodes, q_new];
    tree.parents = [tree.parents, best_parent];
    tree.costs = [tree.costs, best_cost];
    new_node_idx = size(tree.nodes, 2);

    % Rewire: Check if any nearby nodes benefit from going through q_new
    for i = 1:length(near_indices)
        idx = near_indices(i);
        if idx == best_parent
            continue;  % Skip the parent
        end

        new_cost_via_qnew = best_cost + norm(tree.nodes(:, idx) - q_new);

        if new_cost_via_qnew < tree.costs(idx)
            % Check if path from q_new to this node is collision-free
            if is_path_collision_free(q_new, tree.nodes(:, idx), obstacles, safety_margin)
                % Rewire: change parent to new node
                tree.parents(idx) = new_node_idx;
                % Update cost for this node and propagate to children
                tree = propagate_cost_update(tree, idx, new_cost_via_qnew);
            end
        end
    end
end

function tree = propagate_cost_update(tree, node_idx, new_cost)
% Propagate cost update to node and all its descendants

    cost_diff = tree.costs(node_idx) - new_cost;
    tree.costs(node_idx) = new_cost;

    % Find all children and update their costs
    n_nodes = size(tree.nodes, 2);
    for i = 1:n_nodes
        if tree.parents(i) == node_idx
            tree = propagate_cost_update(tree, i, tree.costs(i) - cost_diff);
        end
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
