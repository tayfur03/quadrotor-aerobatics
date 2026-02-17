function [path, info] = rrt_star_radar(start, goal, terrain, threat, params)
%RRT_STAR_RADAR Radar-aware RRT* path planner (Upgraded)
%
% Features:
%   - Dynamic RCS cost (stealthier when appearing smaller to radar)
%   - Log-Barrier radar cost function
%   - Informed RRT* sampling (ellipsoidal)
%   - Shadow Zone sampling bias
%   - Optional hard radar visibility constraint (binary safe/unsafe)
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
beta = get_param(params, 'beta', 100);          % Radar cost weight
gamma = get_param(params, 'gamma', 0.3);        % Preferred-AGL penalty weight
min_clearance = get_param(params, 'min_clearance', 20);
preferred_agl = get_param(params, 'preferred_agl', max(min_clearance + 20, 60));
max_climb = get_param(params, 'max_climb', 30);
max_descent = get_param(params, 'max_descent', 30);

% Shadow sampling bias (probability to reject visible points)
shadow_bias = get_param(params, 'shadow_bias', 0.7);
animate = get_param(params, 'animate', false);
plot_interval = max(1, round(get_param(params, 'plot_interval', 200)));

% Binary visibility mode (visible points/edges become invalid)
radar_hard_constraint = get_param(params, 'radar_hard_constraint', false);
radar_visibility_threshold = get_param(params, 'radar_visibility_threshold', 0.5);

% Optional live plotting handles
anim_axes = [];
if animate
    if isfield(params, 'anim_axes') && ~isempty(params.anim_axes) && isgraphics(params.anim_axes, 'axes')
        anim_axes = params.anim_axes;
    elseif isfield(params, 'anim_fig') && ~isempty(params.anim_fig) && isgraphics(params.anim_fig, 'figure')
        anim_axes = get(params.anim_fig, 'CurrentAxes');
    end

    if isempty(anim_axes) || ~isgraphics(anim_axes, 'axes')
        anim_fig = figure('Name', 'RRT* Live', 'Position', [120, 120, 1000, 760]);
        anim_axes = axes('Parent', anim_fig);
        hold(anim_axes, 'on'); grid(anim_axes, 'on');
        xlabel(anim_axes, 'North [m]');
        ylabel(anim_axes, 'East [m]');
        zlabel(anim_axes, 'Altitude [m]');
        view(anim_axes, 35, 30);
    end
    hold(anim_axes, 'on');
    h_tree = plot3(anim_axes, nan, nan, nan, '-', 'Color', [0.15 0.55 1.0], 'LineWidth', 0.8, ...
        'DisplayName', 'RRT* Tree');
    h_best = plot3(anim_axes, nan, nan, nan, 'm-', 'LineWidth', 2.4, 'DisplayName', 'Best Path');
else
    h_tree = [];
    h_best = [];
end

% Get bounds
if isfield(params, 'bounds')
    bounds = params.bounds;
else
    tb = terrain.bounds;
    % Altitude bounds in NED-D should follow terrain elevation scale.
    if isprop(terrain, 'Z')
        terrain_min = min(terrain.Z(:));
        terrain_max = max(terrain.Z(:));
    else
        % Conservative fallback when Z is unavailable
        terrain_min = 0;
        terrain_max = 500;
    end
    max_flight_alt = get_param(params, 'max_flight_alt', 500);
    min_alt = terrain_min + min_clearance;
    max_alt = terrain_max + max_flight_alt;
    bounds = [tb(1), tb(2), tb(3), tb(4), -max_alt, -min_alt];
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
    goal_updated_this_iter = false;

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
            start, goal, c_best, c_min, x_center, C, shadow_bias, ...
            radar_hard_constraint, radar_visibility_threshold);
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

    % Compute Cost (Distance + Radar + Preferred AGL)
    edge_cost = compute_edge_cost(q_near, q_new, terrain, threat, ...
        alpha, beta, gamma, preferred_agl, radar_hard_constraint, radar_visibility_threshold);

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
                terrain, threat, alpha, beta, gamma, preferred_agl, radar_hard_constraint, radar_visibility_threshold);
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
                terrain, threat, alpha, beta, gamma, preferred_agl, radar_hard_constraint, radar_visibility_threshold);
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
            last_edge = compute_edge_cost(q_new, goal, terrain, threat, alpha, beta, gamma, preferred_agl, ...
                radar_hard_constraint, radar_visibility_threshold);
            total_goal_cost = costs(new_idx) + last_edge;

            if ~goal_reached
                % First time reaching goal: add goal node to tree
                goal_node_idx = size(nodes, 2) + 1;
                nodes(:, goal_node_idx) = goal;
                parents(goal_node_idx) = new_idx;
                costs(goal_node_idx) = total_goal_cost;
                goal_reached = true;
                fprintf('Goal reached at iter %d! Cost: %.2f\n', iter, total_goal_cost);
                goal_updated_this_iter = true;
            elseif total_goal_cost < costs(goal_node_idx)
                % Better path found: UPDATE existing goal node (don't add new one)
                parents(goal_node_idx) = new_idx;
                costs(goal_node_idx) = total_goal_cost;
                fprintf('Goal improved at iter %d! Cost: %.2f\n', iter, total_goal_cost);
                goal_updated_this_iter = true;
            end
        end
    end

    if mod(iter, 500) == 0
        fprintf('RRT* iter %d, tree: %d, best cost: %.2f\n', iter, size(nodes, 2), ternary(isinf(c_best), NaN, c_best));
    end

    if animate && isgraphics(anim_axes, 'axes') && (mod(iter, plot_interval) == 0 || goal_updated_this_iter)
        [x_tree, y_tree, z_tree] = build_tree_lines(nodes, parents);
        set(h_tree, 'XData', x_tree, 'YData', y_tree, 'ZData', z_tree);

        if goal_reached && goal_node_idx > 0
            path_live = extract_path(nodes, parents, goal_node_idx);
            set(h_best, 'XData', path_live(1, :), 'YData', path_live(2, :), 'ZData', -path_live(3, :));
        end

        title(anim_axes, sprintf('RRT* Live | iter=%d | tree=%d | best cost=%.2f', ...
            iter, size(nodes, 2), ternary(isinf(c_best), NaN, c_best)));
        drawnow limitrate nocallbacks;
    end
end

if animate && isgraphics(anim_axes, 'axes')
    [x_tree, y_tree, z_tree] = build_tree_lines(nodes, parents);
    set(h_tree, 'XData', x_tree, 'YData', y_tree, 'ZData', z_tree);
    drawnow;
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

function q = sample_informed(bounds, terrain, threat, min_clearance, start, goal, c_best, c_min, x_center, C, shadow_bias, radar_hard_constraint, radar_visibility_threshold)
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

    % Check shadow bias / hard visibility constraint
    if ~isempty(threat) && threat.computed
        alt = -cand(3);
        terrain_h = terrain.get_height(cand(1), cand(2));
        if alt > terrain_h + min_clearance
            point_risk = 0;
            try
                point_risk = threat.get_risk(cand(1), cand(2), alt);
            catch
                point_risk = 0;
            end
            visible_to_any = point_risk >= radar_visibility_threshold;

            if radar_hard_constraint && visible_to_any
                continue; % Never sample visible points in hard mode
            end

            if ~radar_hard_constraint && rand < shadow_bias && visible_to_any
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

function [x, y, z] = build_tree_lines(nodes, parents)
% Build NaN-separated line arrays for fast tree plotting.

n_nodes = size(nodes, 2);
if n_nodes < 2
    x = nan; y = nan; z = nan;
    return;
end

child_idx = 2:n_nodes;
par_idx = parents(child_idx);
valid = par_idx > 0;
child_idx = child_idx(valid);
par_idx = par_idx(valid);

if isempty(child_idx)
    x = nan; y = nan; z = nan;
    return;
end

n_edges = numel(child_idx);
x = nan(1, 3 * n_edges);
y = nan(1, 3 * n_edges);
z = nan(1, 3 * n_edges);

x(1:3:end) = nodes(1, par_idx);
x(2:3:end) = nodes(1, child_idx);
y(1:3:end) = nodes(2, par_idx);
y(2:3:end) = nodes(2, child_idx);
z(1:3:end) = -nodes(3, par_idx);
z(2:3:end) = -nodes(3, child_idx);
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

function cost = compute_edge_cost(q1, q2, terrain, threat, alpha, beta, gamma, preferred_agl, radar_hard_constraint, radar_visibility_threshold)
% Edge cost = distance + radar (log-barrier + dynamic RCS) + preferred AGL

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

        % In hard mode, any visible sample invalidates the edge
        if radar_hard_constraint && P_base >= radar_visibility_threshold
            cost = inf;
            return;
        end

        if radar_hard_constraint
            continue;
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

% --- Preferred altitude (AGL) cost ---
cost_alt = 0;
if gamma > 0
    n_alt_samples = max(2, ceil(dist / 20));
    err_sq_sum = 0;
    agl_scale = max(preferred_agl, 1);
    for i = 1:n_alt_samples
        t = (i - 1) / (n_alt_samples - 1);
        q = q1 + t * (q2 - q1);
        alt = -q(3);
        terrain_h = terrain.get_height(q(1), q(2));
        agl = alt - terrain_h;
        err_norm = (agl - preferred_agl) / agl_scale;
        err_sq_sum = err_sq_sum + err_norm^2;
    end
    mean_err_sq = err_sq_sum / n_alt_samples;
    cost_alt = gamma * dist * mean_err_sq;
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
