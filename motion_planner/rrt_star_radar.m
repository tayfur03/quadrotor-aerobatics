function [path, info] = rrt_star_radar(start, goal, terrain, threat, params)
%RRT_STAR_RADAR Radar-aware RRT* path planner (Upgraded)
%
% Features:
%   - Dynamic RCS cost (stealthier when appearing smaller to radar)
%   - Log-Barrier radar cost function
%   - Informed RRT* sampling (ellipsoidal)
%   - Shadow Zone sampling bias
%
% Inputs:
%   start   - [3x1] start position [N; E; D] in NED frame
%   goal    - [3x1] goal position [N; E; D] in NED frame
%   terrain - terrain_map object
%   threat  - threat_map object
%   params  - Struct with planning parameters
%
% Outputs:
%   path    - [3 x M] planned path waypoints in NED
%   info    - Struct with planning statistics

tic;

% Default parameters
if nargin < 5
    params = struct();
end

max_iter = get_param(params, 'max_iter', 3000); % Reduced from 5000 due to better sampling
base_step_size = get_param(params, 'base_step_size', get_param(params, 'step_size', 50));
max_step_size = get_param(params, 'max_step_size', base_step_size * 3);
goal_bias = get_param(params, 'goal_bias', 0.1);
rewire_radius = get_param(params, 'rewire_radius', base_step_size * 3);  % Must be > step_size
alpha = get_param(params, 'alpha', 1.2);        % Distance weight
beta = get_param(params, 'beta', 100);         % Radar cost weight
gamma = get_param(params, 'gamma', 0.1);        % Climb-rate cost weight
min_clearance = get_param(params, 'min_clearance', 20);
max_climb = get_param(params, 'max_climb', 30);
max_descent = get_param(params, 'max_descent', 30);

% Shadow sampling bias (probability to reject visible points)
shadow_bias = 0.7;

% Get bounds
if isfield(params, 'bounds')
    bounds = params.bounds;
else
    tb = terrain.bounds;
    bounds = [tb(1), tb(2), tb(3), tb(4), -500, -min_clearance];
end

start = start(:);
goal = goal(:);

% Initialize tree
nodes = start;
parents = 0;
costs = 0;      % Cost from start to node

goal_reached = false;
goal_node_idx = -1;  % Persistent goal node index (reused, never duplicated)
goal_tolerance = base_step_size;

c_min = norm(goal - start);
x_center = (start + goal) / 2;
C = rotation_matrix_3d(start, goal);

% Main RRT* loop
for iter = 1:max_iter

    % Informed Sampling Logic
    c_best = inf;
    if goal_reached
        c_best = costs(goal_node_idx); % Current best cost
    end

    % Sample random point
    if rand < goal_bias
        q_rand = goal;
    else
        % Informed sampling with shadow bias
        q_rand = sample_informed(bounds, terrain, threat, min_clearance, ...
            start, goal, c_best, c_min, x_center, C, shadow_bias);
    end

    % Find nearest node
    [q_near_idx, ~] = find_nearest(nodes, q_rand);
    q_near = nodes(:, q_near_idx);

    % Steer (adaptive step size based on terrain slope)
    q_new = steer_adaptive(q_near, q_rand, base_step_size, max_step_size, ...
        max_climb, max_descent, terrain);

    % Check validity
    if ~is_valid_point(q_new, terrain, min_clearance, bounds)
        continue;
    end

    % Check edge validity (collision)
    if ~is_edge_valid(q_near, q_new, terrain, min_clearance)
        continue;
    end

    % Compute Cost (Dynamic RCS + Log-Barrier + Climb-Rate)
    edge_cost = compute_edge_cost(q_near, q_new, terrain, threat, ...
        alpha, beta, gamma);

    if isinf(edge_cost)
        continue; % Skip if P_det is too high (near 1.0)
    end

    new_cost = costs(q_near_idx) + edge_cost;

    % Rewire step
    near_indices = find_near(nodes, q_new, rewire_radius);

    best_parent = q_near_idx;
    best_cost = new_cost;

    % Connect to best parent
    for i = 1:length(near_indices)
        idx = near_indices(i);
        q_candidate = nodes(:, idx);

        if is_edge_valid(q_candidate, q_new, terrain, min_clearance)
            cand_edge_cost = compute_edge_cost(q_candidate, q_new, ...
                terrain, threat, alpha, beta, gamma);
            cand_cost = costs(idx) + cand_edge_cost;

            if cand_cost < best_cost
                best_parent = idx;
                best_cost = cand_cost;
            end
        end
    end

    % Add node
    new_idx = size(nodes, 2) + 1;
    nodes(:, new_idx) = q_new;
    parents(new_idx) = best_parent;
    costs(new_idx) = best_cost;

    % Rewire neighbors
    for i = 1:length(near_indices)
        idx = near_indices(i);
        if idx == best_parent
            continue;
        end
        q_near_node = nodes(:, idx);

        if is_edge_valid(q_new, q_near_node, terrain, min_clearance)
            rew_edge_cost = compute_edge_cost(q_new, q_near_node, ...
                terrain, threat, alpha, beta, gamma);
            rew_cost = costs(new_idx) + rew_edge_cost;

            if rew_cost < costs(idx)
                cost_delta = costs(idx) - rew_cost;
                parents(idx) = new_idx;
                costs(idx) = rew_cost;
                % Inline cost propagation (can't use helper — MATLAB pass-by-value)
                prop_queue = find(parents == idx);
                while ~isempty(prop_queue)
                    cidx = prop_queue(1);
                    prop_queue(1) = [];
                    costs(cidx) = costs(cidx) - cost_delta;
                    prop_queue = [prop_queue, find(parents == cidx)];
                end
            end
        end
    end

    % Check goal
    dist_to_goal = norm(q_new - goal);
    if dist_to_goal < goal_tolerance
        % Add exact connection to goal if valid
        if is_edge_valid(q_new, goal, terrain, min_clearance)
            last_edge = compute_edge_cost(q_new, goal, terrain, threat, alpha, beta, gamma);
            total_goal_cost = costs(new_idx) + last_edge;

            if ~goal_reached
                % First time reaching goal: add goal node to tree
                goal_node_idx = size(nodes, 2) + 1;
                nodes(:, goal_node_idx) = goal;
                parents(goal_node_idx) = new_idx;
                costs(goal_node_idx) = total_goal_cost;
                goal_reached = true;
                fprintf('Goal reached at iter %d! Cost: %.2f\n', iter, total_goal_cost);
            elseif total_goal_cost < costs(goal_node_idx)
                % Better path found: UPDATE existing goal node (don't add new one)
                parents(goal_node_idx) = new_idx;
                costs(goal_node_idx) = total_goal_cost;
                fprintf('Goal improved at iter %d! Cost: %.2f\n', iter, total_goal_cost);
            end
        end
    end

    if mod(iter, 500) == 0
        fprintf('RRT* iter %d, tree: %d, best cost: %.2f\n', iter, size(nodes, 2), ternary(isinf(c_best), NaN, c_best));
    end
end

% Extract path
if goal_reached
    path = extract_path(nodes, parents, goal_node_idx);
else
    [closest_idx, ~] = find_nearest(nodes, goal);
    path = extract_path(nodes, parents, closest_idx);
    warning('Goal not reached.');
end

planning_time = toc;

% Stats
path_length = compute_path_length(path);
path_risk = compute_path_risk(path, threat);

info.success = goal_reached;
info.iterations = iter;
info.tree_size = size(nodes, 2);
if goal_reached && goal_node_idx > 0
    info.path_cost = costs(goal_node_idx);
else
    info.path_cost = costs(end);
end
info.path_length = path_length;
info.path_risk = path_risk;
info.planning_time = planning_time;

fprintf('RRT* Finished. Time: %.2fs. Length: %.1fm. Risk: %.2f.\n', planning_time, path_length, path_risk);

end

%% Helper Functions

function q = sample_informed(bounds, terrain, threat, min_clearance, start, goal, c_best, c_min, x_center, C, shadow_bias)
% Informed sampling with shadow bias

max_attempts = 100;

for k = 1:max_attempts
    if isinf(c_best)
        % Uniform sampling in bounds
        N = bounds(1) + rand * (bounds(2) - bounds(1));
        E = bounds(3) + rand * (bounds(4) - bounds(3));
        D = bounds(5) + rand * (bounds(6) - bounds(5));
        cand = [N; E; D];
    else
        % Ellipsoidal sampling
        r = [c_best/2; sqrt(c_best^2 - c_min^2)/2; sqrt(c_best^2 - c_min^2)/2];
        L = diag(r);
        x_ball = sample_unit_ball();
        cand = C * L * x_ball + x_center;
    end

    % Check bounds
    if cand(1) < bounds(1) || cand(1) > bounds(2) || ...
            cand(2) < bounds(3) || cand(2) > bounds(4) || ...
            cand(3) < bounds(5) || cand(3) > bounds(6)
        continue;
    end

    % Check shadow bias
    % If point is visible to radar, reject it with probability shadow_bias
    % unless we really struggle to find points
    if ~isempty(threat) && threat.computed && rand < shadow_bias
        alt = -cand(3);
        terrain_h = terrain.get_height(cand(1), cand(2));
        if alt > terrain_h + min_clearance
            visible_to_any = false;
            for i = 1:length(threat.radars)
                radar = threat.radars{i};
                if threat.los.has_los(radar.position, [cand(1); cand(2); alt])
                    visible_to_any = true;
                    break;
                end
            end

            if visible_to_any
                continue; % Reject visible point to favor shadow
            end
        end
    end

    % Check terrain
    terrain_h = terrain.get_height(cand(1), cand(2));
    alt = -cand(3);
    if alt >= terrain_h + min_clearance
        q = cand;
        return;
    end
end

% Fallback if sampling fails
N = bounds(1) + rand * (bounds(2) - bounds(1));
E = bounds(3) + rand * (bounds(4) - bounds(3));
q = [N; E; -(terrain.get_height(N, E) + min_clearance + 20)];
end

function x = sample_unit_ball()
% Uniform sample in 3D unit ball
while true
    x = 2*rand(3,1) - 1;
    if norm(x) <= 1
        return;
    end
end
end

function R = rotation_matrix_3d(a, b)
% Rotation matrix aligning X-axis with vector (b-a)
u = (b - a) / norm(b - a);
id1 = [1; 0; 0];

v = cross(id1, u);
s = norm(v);
c = dot(id1, u);

if s < 1e-6
    R = eye(3);
else
    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + vx + vx^2 * ((1 - c) / s^2);
end
end

function cost = compute_edge_cost(q1, q2, terrain, threat, alpha, beta, gamma)
% Edge cost = distance + radar (log-barrier + dynamic RCS) + climb-rate
%
% The climb-rate cost penalizes the climb angle squared, encouraging
% smooth altitude profiles regardless of terrain. This replaces the
% former preferred_alt penalty which was terrain-agnostic.

dist = norm(q2 - q1);
if dist < 1e-6
    cost = 0;
    return;
end
cost_dist = alpha * dist;

% --- Radar cost (log-barrier + dynamic RCS) ---
cost_radar = 0;

if ~isempty(threat) && threat.computed
    n_samples = max(2, ceil(dist / 20));

    % Flight direction vector (unit)
    v_path = (q2 - q1) / dist;

    sum_log_prob = 0;

    for i = 1:n_samples
        t = (i - 1) / (n_samples - 1);
        q = q1 + t * (q2 - q1);
        alt = -q(3);

        % Base detection probability (for static RCS)
        P_base = 0;
        try
            P_base = threat.get_risk(q(1), q(2), alt);
        catch
            P_base = 0;
        end

        % Dynamic RCS Factor - Find nearest radar and compute aspect angle
        best_dist = inf;
        nearest_radar_idx = 0;
        for r = 1:length(threat.radars)
            r_pos_ned = [threat.radars{r}.position(1); threat.radars{r}.position(2); -threat.radars{r}.position(3)];
            d = norm(r_pos_ned - q);
            if d < best_dist
                best_dist = d;
                nearest_radar_idx = r;
            end
        end

        rcs_mod = 1.0;
        if nearest_radar_idx > 0
            r_pos_ned = [threat.radars{nearest_radar_idx}.position(1); ...
                threat.radars{nearest_radar_idx}.position(2); ...
                -threat.radars{nearest_radar_idx}.position(3)];
            v_los = (r_pos_ned - q) / norm(r_pos_ned - q);
            sin_theta = norm(cross(v_path, v_los));
            rcs_mod = 0.1 + 0.9 * (sin_theta^2);
        end

        P_effective = P_base * rcs_mod;

        % Log-Barrier Cost: C = -log(1 - P)
        P_effective = min(P_effective, 0.99); % Clamp to avoid singularity
        radar_penalty = -log(1 - P_effective);
        sum_log_prob = sum_log_prob + radar_penalty;
    end

    cost_radar = beta * sum_log_prob * (dist / n_samples);
end

% --- Climb-rate cost (AGL-based: penalizes changes relative to terrain) ---
cost_climb = 0;
if gamma > 0
    % Use AGL (above-ground-level) so terrain-following gets no penalty
    alt1 = -q1(3);
    alt2 = -q2(3);
    terrain_h1 = terrain.get_height(q1(1), q1(2));
    terrain_h2 = terrain.get_height(q2(1), q2(2));
    agl1 = alt1 - terrain_h1;
    agl2 = alt2 - terrain_h2;

    delta_agl = abs(agl2 - agl1);  % AGL change, not MSL
    horiz_dist = norm(q2(1:2) - q1(1:2));

    if horiz_dist > 1e-3
        climb_angle = atan2(delta_agl, horiz_dist);  % [rad]
    else
        climb_angle = pi/2;  % Pure vertical move
    end

    % Quadratic penalty: gentle AGL changes ~free, steep AGL changes expensive
    cost_climb = gamma * dist * (climb_angle)^2;
end

cost = cost_dist + cost_radar + cost_climb;
end

% Standard Helper Functions (find_nearest, steer, etc)
function [idx, dist] = find_nearest(nodes, q)
diffs = nodes - q;
distances = vecnorm(diffs, 2, 1);
[dist, idx] = min(distances);
end

function indices = find_near(nodes, q, radius)
diffs = nodes - q;
distances = vecnorm(diffs, 2, 1);
indices = find(distances <= radius);
end

function q_new = steer_adaptive(q_from, q_to, base_step, max_step, ...
                                 max_climb, max_descent, terrain)
% STEER_ADAPTIVE - Step size adapts to local terrain slope
%
% On steep terrain, uses shorter steps so the tree can follow the surface.
% On flat terrain, uses longer steps for efficiency.
% Logistic blending centered at ~20 deg slope.

direction = q_to - q_from;
dist = norm(direction);
if dist < 1e-6, q_new = q_to; return; end
direction = direction / dist;

% --- Adaptive step size based on local terrain slope ---
h0 = terrain.get_height(q_from(1), q_from(2));
probe_dist = min(10, base_step / 5);  % Small probe for gradient
h_north = terrain.get_height(q_from(1) + probe_dist, q_from(2));
h_east  = terrain.get_height(q_from(1), q_from(2) + probe_dist);

slope_n = (h_north - h0) / probe_dist;
slope_e = (h_east  - h0) / probe_dist;
slope_mag = sqrt(slope_n^2 + slope_e^2);  % Terrain gradient magnitude
slope_angle = atan(slope_mag);             % [rad]

% Logistic blend: flat → max_step, steep → base_step
blend = 1 / (1 + exp(-10 * (slope_angle - deg2rad(20))));
step_size = max_step * (1 - blend) + base_step * blend;

% Apply step
step = min(dist, step_size) * direction;

% --- Climb angle clamping ---
horizontal_dist = norm(step(1:2));
vertical_dist = step(3);

if horizontal_dist > 1e-6
    max_climb_dist = horizontal_dist * tand(max_climb);
    max_descent_dist = horizontal_dist * tand(max_descent);
    if -vertical_dist > max_climb_dist
        step(3) = -max_climb_dist;
    elseif vertical_dist > max_descent_dist
        step(3) = max_descent_dist;
    end
end

q_new = q_from + step;
end

function valid = is_valid_point(q, terrain, min_clearance, bounds)
if q(1)<bounds(1) || q(1)>bounds(2) || q(2)<bounds(3) || q(2)>bounds(4) || q(3)<bounds(5) || q(3)>bounds(6)
    valid = false; return;
end
alt = -q(3);
if alt < terrain.get_height(q(1), q(2)) + min_clearance
    valid = false; return;
end
valid = true;
end

function valid = is_edge_valid(q1, q2, terrain, min_clearance)
dist = norm(q2 - q1);
n_samples = max(3, ceil(dist / 10));
for i = 1:n_samples
    t = (i-1)/(n_samples-1);
    q = q1 + t*(q2-q1);
    if -q(3) < terrain.get_height(q(1), q(2)) + min_clearance
        valid = false; return;
    end
end
valid = true;
end

function path = extract_path(nodes, parents, goal_idx)
path = [];
idx = goal_idx;
while idx > 0
    path = [nodes(:, idx), path];
    idx = parents(idx);
end
end

function len = compute_path_length(path)
len = 0;
for i=2:size(path,2), len = len + norm(path(:,i)-path(:,i-1)); end
end

function risk = compute_path_risk(path, threat)
risk = 0;
if isempty(threat) || ~threat.computed, return; end
for i=1:size(path,2)
    risk = risk + threat.get_risk(path(1,i), path(2,i), -path(3,i));
end
end

function val = get_param(params, name, default)
if isfield(params, name), val = params.(name); else, val = default; end
end

function result = ternary(condition, true_val, false_val)
if condition, result = true_val; else, result = false_val; end
end