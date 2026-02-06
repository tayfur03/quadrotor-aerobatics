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
step_size = get_param(params, 'step_size', 50);
goal_bias = get_param(params, 'goal_bias', 0.1);
rewire_radius = get_param(params, 'rewire_radius', 100);
alpha = get_param(params, 'alpha', 1.2);        % Distance weight
beta = get_param(params, 'beta', 100);         % Radar cost weight
gamma = get_param(params, 'gamma', 0.1);        % Altitude cost weight
%lambda = get_param(params,'lambda', 100);       % Rate of Altitude Change weight
min_clearance = get_param(params, 'min_clearance', 20);
preferred_alt = get_param(params, 'preferred_alt', 50);
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
goal_node_idx = -1;
goal_tolerance = step_size;

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

    % Steer
    q_new = steer(q_near, q_rand, step_size, max_climb, max_descent);

    % Check validity
    if ~is_valid_point(q_new, terrain, min_clearance, bounds)
        continue;
    end

    % Check edge validity (collision)
    if ~is_edge_valid(q_near, q_new, terrain, min_clearance)
        continue;
    end

    % Compute Cost (Dynamic RCS + Log-Barrier)
    edge_cost = compute_edge_cost(q_near, q_new, terrain, threat, ...
        alpha, beta, gamma, preferred_alt);

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
                terrain, threat, alpha, beta, gamma, preferred_alt);
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
                terrain, threat, alpha, beta, gamma, preferred_alt);
            rew_cost = costs(new_idx) + rew_edge_cost;

            if rew_cost < costs(idx)
                parents(idx) = new_idx;
                costs(idx) = rew_cost;
                % Propagate cost updates? (Skipped for performance in this demo)
            end
        end
    end

    % Check goal
    dist_to_goal = norm(q_new - goal);
    if dist_to_goal < goal_tolerance
        % Add exact connection to goal if valid
        if is_edge_valid(q_new, goal, terrain, min_clearance)
            last_edge = compute_edge_cost(q_new, goal, terrain, threat, alpha, beta, gamma, preferred_alt);
            total_goal_cost = costs(new_idx) + last_edge;

            % Add goal node explicitly (or just mark this node as close enough)
            % Let's add the goal node to the tree to be precise
            goal_idx_in_tree = size(nodes, 2) + 1;
            nodes(:, goal_idx_in_tree) = goal;
            parents(goal_idx_in_tree) = new_idx;
            costs(goal_idx_in_tree) = total_goal_cost;
            disp(total_goal_cost)

            if ~goal_reached || total_goal_cost < costs(goal_node_idx)
                goal_reached = true;
                goal_node_idx = goal_idx_in_tree;
                
                fprintf('Goal reached at iter %d! Cost: %.2f\n', iter, total_goal_cost);
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
disp(goal_node_idx)
info.path_cost = ternary(goal_reached, costs(goal_node_idx), costs(end));
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

function cost = compute_edge_cost(q1, q2, terrain, threat, alpha, beta, gamma, preferred_alt)
% Updated Cost Function

dist = norm(q2 - q1);
if dist < 1e-6
    cost = 0;
    return;
end
cost_dist = alpha * dist;

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

        % Dynamic RCS Factor
        RCS_factor = 1.0;

        % Iterate radars to find worst-case alignment?
        % Or average alignment?
        % Let's consider the dominant radar (max P_det)

        max_P_scaled = 0;

        for r = 1:length(threat.radars)
            radar = threat.radars{r};
            % Vector from UAV to Radar
            % radar.position is [N; E; Alt] (Alt is +up)
            % q is [N; E; D] -> Position is [N; E; -D]
            pos_uav = [q(1); q(2); -q(3)];
            v_to_radar = radar.position - pos_uav;
            d_radar = norm(v_to_radar);

            if d_radar > 1e-3
                v_to_radar = v_to_radar / d_radar;

                % Angle theta between velocity and LOS
                cos_theta = dot(v_path, v_to_radar); % But v_path is in NED?
                % v_path is in NED [N; E; D]
                % radar.position is usually [N; E; Alt].
                % We need consistent frames.
                % RRT nodes D is +Down.
                % radar.position(3) is Altitude (+Up).
                % So let's convert radar pos to NED?
                % No, let's convert everything to NED.

                p_radar_ned = [radar.position(1); radar.position(2); -radar.position(3)];
                v_los_ned = p_radar_ned - q; % Vector from UAV to Radar in NED
                v_los_ned = v_los_ned / norm(v_los_ned);

                % Angle calculation
                % Facing radar: v_path aligns with v_los_ned (theta = 0)
                % or opposite (theta = 180).
                % Stealthiest when facing -> RCS min at 0, 180.
                % RCS max at 90.

                sin_theta = norm(cross(v_path, v_los_ned));

                % Factor: 0.1 (min) to 1.0 (max)
                % RCS_mod = 0.1 + 0.9 * sin_theta^2;
                % Let's use a simpler multiplier for P_det
                % P_det_new ~ P_det_old * RCS_mod

                % We don't have P_det individual per radar here easily without querying radar object.
                % threat.get_risk returns max(P_i).
                % Let's approximate: apply the factor to the aggregate risk based on nearest radar?
                % Better: Re-calculate P_det is too expensive.
                % Assumption: The returned risk is dominated by one radar.
                % We apply the aspect angle of the nearest radar?
                % Or assume risk is from the radar we are looking at.

                % Let's average the aspect factor weighted by distance?
                % Or just take the worst case?

                current_factor = 0.2 + 0.8 * (sin_theta^2); % 20% to 100%

                % How to know WHICH radar provides the risk?
                % We don't.
                % Let's compute a weighted factor based on distances (closer radar matters more).
                % weight = 1 / d_radar^2
                % Refined approach: Just calculate factor for ALL radars and average it?
                % No, if I'm far from radar A but close to B, B's angle matters.

                % Let's use the aspect angle relative to the "Threat Center of Mass" or similar?
                % No, let's just re-query detection prob efficiently?
                % No, too slow.

                % Heuristic: Find nearest radar. Use its aspect angle.
                % This is fast and likely correct (risk is dominated by nearest).
            end
        end

        % Find nearest radar
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

        % Log-Barrier Cost
        % C = -log(1 - P)
        % If P -> 1, Cost -> Inf.
        % But we need to handle P >= 1 (clamped).
        P_effective = min(P_effective, 0.99); % Clamp to avoiding singularity

        radar_penalty = -log(1 - P_effective);

        sum_log_prob = sum_log_prob + radar_penalty;
    end

    cost_radar = beta * sum_log_prob * (dist / n_samples);
end

% Altitude cost
cost_alt = 0;
if gamma > 0
    % ... (Same as before)
    n_samples = max(2, ceil(dist / 20));
    sum_alt = 0;
    for i = 1:n_samples
        t = (i - 1) / (n_samples - 1);
        q = q1 + t * (q2 - q1);
        alt = -q(3);
        terrain_h = terrain.get_height(q(1), q(2));
        agl = alt - terrain_h;
        sum_alt = sum_alt + abs(agl - preferred_alt)/preferred_alt;
    end
    cost_alt = gamma * sum_alt * dist / n_samples;
end

cost = cost_dist + cost_radar + cost_alt;
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

function q_new = steer(q_from, q_to, step_size, max_climb, max_descent)
direction = q_to - q_from;
dist = norm(direction);
if dist < 1e-6, q_new = q_to; return; end

direction = direction / dist;
if dist > step_size, direction = direction * step_size; else, direction = direction * dist; end

horizontal_dist = norm(direction(1:2));
vertical_dist = direction(3);

if horizontal_dist > 1e-6
    max_climb_dist = horizontal_dist * tand(max_climb);
    max_descent_dist = horizontal_dist * tand(max_descent);
    if -vertical_dist > max_climb_dist
        direction(3) = -max_climb_dist;
    elseif vertical_dist > max_descent_dist
        direction(3) = max_descent_dist;
    end
end
q_new = q_from + direction;
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