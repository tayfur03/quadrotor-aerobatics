function [waypoints, timePoints, velBC, accBC, maneuver_info] = aerobatic_maneuvers(maneuver_type, params)
%AEROBATIC_MANEUVERS Generate waypoints and boundary conditions for aerobatic maneuvers.
%
% Generates trajectory data compatible with MATLAB's minsnappolytraj function.
% Based on aerobatic maneuvers from Tal et al. (2023) "Aerobatic Trajectory
% Generation for a VTOL Fixed-Wing Aircraft Using Differential Flatness".
%
% Inputs:
%   maneuver_type - String specifying the maneuver:
%                   'straight_line'   - Straight flight between two points (for yaw testing)
%                   'vertical_loop'   - Vertical loop in N-D plane (Fig. 3 of paper)
%                   'barrel_roll'     - Helical barrel roll
%                   'immelmann'       - Half loop + half roll (direction reversal)
%                   'split_s'         - Half roll + half loop (opposite of Immelmann)
%                   'climbing_turn'   - Banked turn with altitude change
%                   'figure_eight'    - Horizontal figure-8 pattern
%                   'hippodrome'      - Hippodrome (Stadium) shaped circuit
%
%   params - Struct with maneuver parameters:
%            .radius       - Loop/turn radius [m] (default: 3)
%            .V            - Nominal speed [m/s] (default: 5)
%            .start_pos    - Starting position [3x1] (default: [0;0;-3])
%            .end_pos      - End position for straight_line [3x1]
%            .direction    - 'cw' or 'ccw' for loop direction (default: 'cw')
%            .num_rolls    - Number of rolls for barrel roll (default: 1)
%            .oscillate    - true/false for back-and-forth motion (default: false)
%            .num_laps     - Number of back-and-forth laps (default: 1)
%
% Outputs:
%   waypoints     - [3 x N] waypoint positions in NED frame
%   timePoints    - [1 x N] time at each waypoint
%   velBC         - [3 x N] velocity boundary conditions (NaN = free)
%   accBC         - [3 x N] acceleration boundary conditions (NaN = free)
%   maneuver_info - Struct with additional info:
%                   .name          - Maneuver name
%                   .total_time    - Total maneuver duration
%                   .yaw_mode      - Recommended yaw planner mode
%                   .yaw_params    - Parameters for yaw planner
%                   .description   - Text description
%
% Example:
%   params.radius = 1.0;  % 1m radius as in paper
%   params.V = 3.5;       % ~3.5 m/s as in paper
%   [wp, tp, vBC, aBC, info] = aerobatic_maneuvers('vertical_loop', params);
%   [xref, vref, aref, jref, sref] = minsnappolytraj(wp, tp, 1000, ...
%       'VelocityBoundaryCondition', vBC, ...
%       'AccelerationBoundaryCondition', aBC);

% Default parameters
if nargin < 2
    params = struct();
end

radius = get_param(params, 'radius', 1.0);
V = get_param(params, 'V', 3.5);
start_pos = get_param(params, 'start_pos', [0; 0; -3]);
direction = get_param(params, 'direction', 'cw');
num_rolls = get_param(params, 'num_rolls', 1);
oscillate = get_param(params, 'oscillate', false);
num_laps = get_param(params, 'num_laps', 1);

start_pos = start_pos(:);

switch lower(maneuver_type)
    case 'straight_line'
        end_pos = get_param(params, 'end_pos', start_pos + [10; 0; 0]);
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_straight_line(start_pos, end_pos, V, oscillate, num_laps);

    case 'vertical_loop'
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_vertical_loop(start_pos, radius, V, direction);

    case 'barrel_roll'
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_barrel_roll(start_pos, radius, V, num_rolls);

    case 'immelmann'
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_immelmann(start_pos, radius, V);

    case 'split_s'
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_split_s(start_pos, radius, V);

    case 'climbing_turn'
        climb_height = get_param(params, 'climb_height', 1.0);
        turn_angle = get_param(params, 'turn_angle', 270);
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_climbing_turn(start_pos, radius, V, climb_height, turn_angle);

    case 'figure_eight'
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_figure_eight(start_pos, radius, V);

    case 'hippodrome'
        length_N = get_param(params, 'length', 12.0);
        width_E  = get_param(params, 'width', 6.0);
        [waypoints, timePoints, velBC, accBC, maneuver_info] = ...
            generate_hippodrome(start_pos, length_N, width_E, V);

    otherwise
        error('Unknown maneuver type: %s', maneuver_type);
end
end

%% ============== STRAIGHT LINE (for yaw mode testing) ==============

function [wp, tp, velBC, accBC, info] = generate_straight_line(start, finish, V, oscillate, num_laps)
%GENERATE_STRAIGHT_LINE Straight flight between two points
%
% Useful for testing different yaw modes:
%   - constant: fixed yaw during flight
%   - tangent: yaw aligned with velocity direction
%   - constant_rate: rotating yaw during flight
%   - poi: tracking a point of interest
%
% With oscillate=true, performs back-and-forth motion.

start = start(:);
finish = finish(:);

dist = norm(finish - start);
dir_vec = (finish - start) / dist;

if oscillate && num_laps > 1
    % Back-and-forth motion
    % Each lap: start -> finish -> start
    wp_list = start;
    for lap = 1:num_laps
        % Add midpoint for smoother turns
        wp_list = [wp_list, finish, start];
    end
    wp = wp_list;
    n_wp = size(wp, 2);

    % Time allocation
    segment_time = dist / V;
    total_time = segment_time * 2 * num_laps;
    tp = linspace(0, total_time, n_wp);

    % Velocity BCs: alternating directions
    velBC = nan(3, n_wp);
    velBC(:, 1) = [0; 0; 0];  % Start from hover
    for i = 2:n_wp-1
        if mod(i, 2) == 0
            velBC(:, i) = dir_vec * V;   % Forward
        else
            velBC(:, i) = -dir_vec * V;  % Backward
        end
    end
    velBC(:, end) = [0; 0; 0];  % End at hover

    info.description = sprintf('Back-and-forth %.1fm x %d laps at %.1fm/s', dist, num_laps, V);
else
    % Simple point-to-point
    % 5 waypoints for smooth trajectory: start, 1/4, 1/2, 3/4, end
    fractions = [0, 0.25, 0.5, 0.75, 1.0];
    wp = start + fractions .* (finish - start);

    % Time allocation based on distance
    total_time = dist / V + 2.0;  % Add time for accel/decel
    tp = linspace(0, total_time, length(fractions)); % Distance-proportional for straight line

    % Velocity BCs - Free interior waypoints for smoother snap optimization
    velBC = nan(3, 5);
    velBC(:, 1) = [0; 0; 0];           % Start from hover
    velBC(:, end) = [0; 0; 0];         % End at hover

    info.description = sprintf('Straight line %.1fm at %.1fm/s', dist, V);
end

% Acceleration BCs
n_wp = size(wp, 2);
accBC = nan(3, n_wp);
accBC(:, 1) = [0; 0; 0];    % Start from rest
accBC(:, end) = [0; 0; 0];  % End at rest

info.name = 'Straight Line';
info.total_time = tp(end);
info.yaw_mode = 'constant';  % Default, but user should override
info.yaw_params = struct();
end

%% ============== VERTICAL LOOP (Fig. 3 of paper) ==============

function [wp, tp, velBC, accBC, info] = generate_vertical_loop(start, R, V, direction)
%GENERATE_VERTICAL_LOOP Vertical loop in N-D plane with approach/exit segments
%
% Structure:
%   1. APPROACH: Start from hover, accelerate to V(1) toward loop entry
%   2. LOOP: 4 quadrants (30, 90, 180, 270, 360 deg)
%   3. EXIT: Decelerate from V(end) back to hover
%
% V can be a scalar (constant speed) or a vector of length 6 specifying speed at:
%   [0 deg, 30 deg, 90 deg, 180 deg, 270 deg, 360 deg]

% Define angles
if strcmpi(direction, 'ccw')
    theta_deg = [0, -90, -180, -270, -360];
else
    % Standard Pull-Up Loop: 0 (Entry) -> 360 (Exit)
    theta_deg = [0, 90, 180, 270, 360];
end
theta = deg2rad(theta_deg);
n_points_loop = length(theta); % Should be 6

% Handle Velocity Input
if isscalar(V)
    V_prof = repmat(V, 1, n_points_loop);
else
    if length(V) ~= n_points_loop
        error('For vertical loop, V must be scalar or vector of length %d (corresponding to angles: %s)', ...
            n_points_loop, mat2str(theta_deg));
    end
    V_prof = V;
end

V_entry = V_prof(1);
V_exit_start = V_prof(end);

% Approach/exit distance based on entry/exit speeds
a_accel = 5.0;  % Acceleration during approach [m/s^2]
approach_dist = V_entry^2 / (2 * a_accel);
approach_dist = max(approach_dist, 1.5);  % Minimum 1.5m approach

approach_time = V_entry / a_accel;
approach_time = max(approach_time, 1.0);  % Minimum 1s approach

% Loop bottom position (where the loop starts, 0 deg)
loop_bottom = start + [approach_dist; 0; 0];

% Circle center is above the loop bottom (in NED, up is -D)
center = loop_bottom + [0; 0; -R];

% === APPROACH WAYPOINTS ===
wp_approach = [start, loop_bottom];
vel_approach = [[0; 0; 0], [V_entry; 0; 0]];
acc_approach = [[0; 0; 0],  nan(3,1)];
t_approach = [0, approach_time];

% === LOOP WAYPOINTS ===
% Waypoint positions on circle (N-D plane)
% Note: wp_loop includes 0 deg point which duplicates loop_bottom
N_loop = center(1) + R * sin(theta);
E_loop = center(2) * ones(size(theta));
D_loop = center(3) + R * cos(theta);

wp_loop = [N_loop; E_loop; D_loop];

% Velocity Constraints (Tangent to circle)
v_N = V_prof .* cos(theta);
v_D = -V_prof .* sin(theta);
v_E = zeros(size(theta));

vel_loop = [v_N; v_E; v_D];

% Acceleration Constraints: Free (NaN)
acc_loop = nan(3, n_points_loop);

% === EXIT WAYPOINTS ===
loop_exit = wp_loop(:, end);
% Calculate exit distance based on V_exit_start
exit_dist = V_exit_start^2 / (2 * a_accel);
exit_dist = max(exit_dist, 1.5);
wp_exit = loop_exit + [exit_dist; 0; 0];
vel_exit = [0; 0; 0];
acc_exit = [0; 0; 0];

% === COMBINE ALL WAYPOINTS ===
% wp_approach ends at loop_bottom (0 deg)
% wp_loop starts at 0 deg, so wp_loop(:,1) is duplicate of wp_approach(:,end)
% We use wp_approach then wp_loop(:,2:end)
wp = [wp_approach, wp_loop(:,2:end), wp_exit];
velBC = [vel_approach, vel_loop(:,2:end), vel_exit];
accBC = [acc_approach, acc_loop(:,2:end), acc_exit];

% === TIME ALLOCATION ===
% Calculate time for loop segments based on varying speed
t_loop_segments = zeros(1, n_points_loop);
t_loop_segments(1) = 0; % Relative to loop start

for i = 2:n_points_loop
    % Arc length between points
    d_theta = abs(theta(i) - theta(i-1));
    seg_dist = d_theta * R;

    % Average speed in segment
    v_avg = (V_prof(i) + V_prof(i-1)) / 2;
    if v_avg < 0.1, v_avg = 0.1; end % Avoid divide by zero

    dt_seg = seg_dist / v_avg;
    t_loop_segments(i) = t_loop_segments(i-1) + dt_seg;
end

% Exit time
exit_time = V_exit_start / a_accel;
exit_time = max(exit_time, 1.0);

% Assemble time points
% t_approach ends at approach_time
t_loop_abs = approach_time + t_loop_segments;
t_exit_abs = t_loop_abs(end) + exit_time;

% tp corresponds to [start, loop 0 (shared), loop 30, ..., loop 360, exit]
tp = [t_approach(1), t_loop_abs, t_exit_abs];

total_time = tp(end);

a_cent_max = max(V_prof.^2) / R;
info.name = 'Vertical Loop (Variable V)';
info.total_time = total_time;
info.yaw_mode = 'coordinated';
info.yaw_params = struct('psi_constant', 0);
info.description = sprintf('Vertical loop R=%.1fm, Max V=%.1fm/s', R, max(V_prof));
info.center = center;
info.radius = R;
end

%% ============== BARREL ROLL ==============

function [wp, tp, velBC, accBC, info] = generate_barrel_roll(start, R, V, num_rolls)
%GENERATE_BARREL_ROLL Helical barrel roll trajectory with approach/exit
%
% The barrel roll is a helical path: forward motion + rotation around velocity axis
% The drone rolls while maintaining forward progress.
%
% Structure:
%   1. APPROACH: Hover -> accelerate to cruise speed V
%   2. ROLL: Helical barrel roll maneuver
%   3. EXIT: Decelerate back to hover
%
% Key insight from paper: For rolling motion, yaw changes at constant rate.

% Approach/exit parameters
a_accel = 2.0;  % Acceleration [m/s^2]
approach_dist = V^2 / (2 * a_accel);
approach_dist = max(approach_dist, 1.5);
approach_time = V / a_accel;
approach_time = max(approach_time, 1.0);

% Roll parameters
helix_pitch = 2 * pi * R;  % Forward distance per roll
roll_time = (2 * pi * R) / V;
maneuver_time = roll_time * num_rolls;
omega = 2 * pi * num_rolls / maneuver_time;  % Angular rate

% === APPROACH WAYPOINTS ===
wp_approach = [start, start + [approach_dist/2; 0; 0]];
vel_approach = [[0; 0; 0], [V/2; 0; 0]];
acc_approach = [[0; 0; 0], nan(3,1)];

% === ROLL WAYPOINTS ===
roll_start = start + [approach_dist; 0; 0];
points_per_roll = 4;
n_roll_points = num_rolls * points_per_roll + 1;
theta = linspace(0, 2*pi*num_rolls, n_roll_points);

% Helix parametrization
N_roll = roll_start(1) + (V * maneuver_time / (2*pi*num_rolls)) * theta;
E_roll = roll_start(2) + R * sin(theta);
D_roll = roll_start(3) - R * (1 - cos(theta));

wp_roll = [N_roll; E_roll; D_roll];

% Velocity: tangent to helix
vN_roll = V * ones(1, n_roll_points);
vE_roll = R * omega * cos(theta);
vD_roll = R * omega * sin(theta);
vel_roll = [vN_roll; vE_roll; vD_roll];

acc_roll = nan(3, n_roll_points);

% === EXIT WAYPOINTS ===
roll_end = wp_roll(:, end);
wp_exit = [roll_end + [approach_dist/2; 0; 0], roll_end + [approach_dist; 0; 0]];
vel_exit = [[V/2; 0; 0], [0; 0; 0]];
acc_exit = [nan(3,1), [0; 0; 0]];

% === COMBINE ===
wp = [wp_approach, wp_roll, wp_exit];
velBC = [vel_approach, vel_roll, vel_exit];
accBC = [acc_approach, acc_roll, acc_exit];

% === TIME ALLOCATION ===
t_approach = [0, approach_time/2];
t_roll = approach_time + linspace(0, maneuver_time, n_roll_points);
t_exit = t_roll(end) + [approach_time/2, approach_time];

tp = [t_approach, t_roll, t_exit];
total_time = tp(end);

info.name = sprintf('Barrel Roll x%d', num_rolls);
info.total_time = total_time;
info.yaw_mode = 'tangent';
info.yaw_params = struct('last_psi', 0);
info.description = sprintf('%d barrel roll(s), R=%.1fm at %.1fm/s, with approach/exit', ...
    num_rolls, R, V);
end

%% ============== IMMELMANN TURN ==============

function [wp, tp, velBC, accBC, info] = generate_immelmann(start, R, V)
%GENERATE_IMMELMANN Immelmann turn with approach/exit
%
% Structure:
%   1. APPROACH: Hover -> accelerate North
%   2. HALF LOOP: Pull up 180 degrees
%   3. EXIT: Decelerate flying South (direction reversed, altitude +2R)

% Approach/exit parameters
a_accel = 2.0;
approach_dist = V^2 / (2 * a_accel);
approach_dist = max(approach_dist, 1.5);
approach_time = V / a_accel;
approach_time = max(approach_time, 1.0);

% === APPROACH ===
% Removed intermediate waypoint for smoother trajectory
wp_approach = [start, start + [approach_dist; 0; 0]];
vel_approach = [[0; 0; 0], [V; 0; 0]];
acc_approach = [[0; 0; 0], nan(3,1)];

% === HALF LOOP ===
loop_entry = start + [approach_dist; 0; 0];
center = loop_entry + [0; 0; -R];

theta_deg = [0, 90, 180];
theta = deg2rad(theta_deg);

N_loop = center(1) + R * sin(theta);
E_loop = center(2) * ones(size(theta));
D_loop = center(3) + R * cos(theta);
wp_loop = [N_loop; E_loop; D_loop];

vN_loop = V * cos(theta);
vE_loop = zeros(size(theta));
vD_loop = -V * sin(theta);
vel_loop = [vN_loop; vE_loop; vD_loop];

a_cent = V^2 / R;
% Optimization: Leave acceleration free (NaN) to allow smooth roll entry/exit
acc_loop = nan(3, 3);
% acc_loop(:, 1) = [0; 0; -a_cent];
% acc_loop(:, 2) = [-a_cent; 0; 0];
% acc_loop(:, 3) = [0; 0; a_cent];

% === EXIT (flying South, -N direction) ===
loop_exit = wp_loop(:, end);
% Removed intermediate waypoint
wp_exit = [loop_exit + [-approach_dist; 0; 0]];
vel_exit = [[0; 0; 0]];
acc_exit = [[0; 0; 0]];

% === COMBINE ===
% Note: wp_loop includes entry/exit points.
% We need to merge carefully.
% wp_approach: [Start, Entry]
% wp_loop: [Entry, Mid, Exit]
% wp_exit: [Stop]
%
% Result: [Start, Entry, Mid, Exit, Stop]

wp = [wp_approach(:,1), wp_loop, wp_exit];
velBC = [vel_approach(:,1), vel_loop, vel_exit];
accBC = [acc_approach(:,1), acc_loop, acc_exit];

% === TIME ===
half_loop_time = (pi * R) / V;
t_approach = 0; % Start time
t_entry = approach_time;
t_loop_rel = [0, half_loop_time/2, half_loop_time];
t_loop = t_entry + t_loop_rel;
t_exit = t_loop(end) + approach_time;

tp = [t_approach, t_loop, t_exit];

info.name = 'Immelmann Turn';
info.total_time = tp(end);
info.yaw_mode = 'coordinated'; % Changed to coordinated for smooth yaw handling
info.yaw_params = struct('psi_constant', 0);
info.description = sprintf('Immelmann R=%.1fm, V=%.1fm/s, gains %.1fm altitude', R, V, 2*R);
end

%% ============== SPLIT-S ==============

function [wp, tp, velBC, accBC, info] = generate_split_s(start, R, V)
%GENERATE_SPLIT_S Split-S with approach/exit
%
% Structure:
%   1. APPROACH: Hover -> accelerate North
%   2. HALF LOOP DOWN: Pull through downward 180 degrees
%   3. EXIT: Decelerate flying South (direction reversed, altitude -2R)

% Approach/exit parameters
a_accel = 2.0;
approach_dist = V^2 / (2 * a_accel);
approach_dist = max(approach_dist, 1.5);
approach_time = V / a_accel;
approach_time = max(approach_time, 1.0);

% === APPROACH ===
% Removed intermediate waypoint
wp_approach = [start, start + [approach_dist; 0; 0]];
vel_approach = [[0; 0; 0], [V; 0; 0]];
% Constraint: Force Y and Z acceleration to 0 at entry to prevent 'anticipation' curve
% X acceleration is left free (NaN)
acc_approach = [[0; 0; 0], [nan; 0; 0]];

% === HALF LOOP DOWN ===
loop_entry = start + [approach_dist; 0; 0];
center = loop_entry + [0; 0; R];  % Center below

% Corrected: 0 -> 90 -> 180 (Forward loop down)
theta_deg = [0, 90, 180];
theta = deg2rad(theta_deg);

N_loop = center(1) + R * sin(theta);
E_loop = center(2) * ones(size(theta));
D_loop = center(3) - R * cos(theta);
wp_loop = [N_loop; E_loop; D_loop];

vN_loop = V * cos(theta);
vE_loop = zeros(size(theta));
vD_loop = V * sin(theta);
vel_loop = [vN_loop; vE_loop; vD_loop];

a_cent = V^2 / R;
% Optimization: Leave acceleration free (NaN) to allow smooth roll entry
acc_loop = nan(3, 3);
% acc_loop(:, 1) = [0; 0; a_cent];
% acc_loop(:, 2) = [-a_cent; 0; 0];
% acc_loop(:, 3) = [0; 0; -a_cent];

% === EXIT (flying South) ===
loop_exit = wp_loop(:, end);
% Removed intermediate waypoint
wp_exit = [loop_exit + [-approach_dist; 0; 0]];
vel_exit = [[0; 0; 0]];
acc_exit = [[0; 0; 0]];

% === COMBINE ===
wp = [wp_approach(:,1), wp_loop, wp_exit];
velBC = [vel_approach(:,1), vel_loop, vel_exit];
accBC = [acc_approach(:,1), acc_loop, acc_exit];

% === TIME ===
half_loop_time = (pi * R) / V;
t_approach = 0;
t_entry = approach_time;
t_loop_rel = [0, half_loop_time/2, half_loop_time];
t_loop = t_entry + t_loop_rel;
t_exit = t_loop(end) + approach_time;

tp = [t_approach, t_loop, t_exit];

info.name = 'Split-S';
info.total_time = tp(end);
info.yaw_mode = 'coordinated'; % Changed to coordinated
info.yaw_params = struct();
info.description = sprintf('Split-S R=%.1fm at %.1fm/s, loses %.1fm altitude', R, V, 2*R);
end

%% ============== CLIMBING TURN ==============

function [wp, tp, velBC, accBC, info] = generate_climbing_turn(start, R, V, dh, turn_angle)
%GENERATE_CLIMBING_TURN Climbing turn with approach/exit
%
% Structure:
%   1. APPROACH: Hover -> accelerate North
%   2. TURN: Banked turn with altitude change
%   3. EXIT: Decelerate to hover

% Approach/exit parameters
a_accel = 2.0;
approach_dist = V^2 / (2 * a_accel);
approach_dist = max(approach_dist, 1.5);
approach_time = V / a_accel;
approach_time = max(approach_time, 1.0);

% === APPROACH ===
wp_approach = [start, start + [approach_dist/2; 0; 0]];
vel_approach = [[0; 0; 0], [V/2; 0; 0]];
acc_approach = [[0; 0; 0], nan(3,1)];

% === TURN ===
turn_entry = start + [approach_dist; 0; 0];
turn_rad = deg2rad(turn_angle);

% Number of waypoints based on turn angle
if turn_angle <= 90
    n_turn = 2;
elseif turn_angle <= 180
    n_turn = 3;
elseif turn_angle <= 270
    n_turn = 4;
else
    n_turn = 5;
end

theta = linspace(0, turn_rad, n_turn);

% Circle center to the East of entry
center_N = turn_entry(1);
center_E = turn_entry(2) + R;

N_turn = center_N + R * sin(theta);
E_turn = center_E - R * cos(theta);
D_turn = turn_entry(3) + linspace(0, -dh, n_turn);
wp_turn = [N_turn; E_turn; D_turn];

% Velocity tangent to turn
arc_length = turn_rad * R;
horizontal_speed = V * cos(atan2(dh, arc_length));
vertical_speed = -dh / (arc_length / V);

vN_turn = horizontal_speed * cos(theta);
vE_turn = horizontal_speed * sin(theta);
vD_turn = vertical_speed * ones(size(theta));
vel_turn = [vN_turn; vE_turn; vD_turn];

% Centripetal acceleration
a_cent = horizontal_speed^2 / R;
acc_turn = nan(3, n_turn);
for i = 1:n_turn
    acc_dir_N = center_N - N_turn(i);
    acc_dir_E = center_E - E_turn(i);
    acc_norm = sqrt(acc_dir_N^2 + acc_dir_E^2);
    if acc_norm > 1e-6
        acc_turn(1, i) = a_cent * acc_dir_N / acc_norm;
        acc_turn(2, i) = a_cent * acc_dir_E / acc_norm;
    end
    acc_turn(3, i) = 0;
end

% === EXIT ===
turn_exit = wp_turn(:, end);
exit_dir = vel_turn(:, end) / norm(vel_turn(:, end));
wp_exit = [turn_exit + exit_dir * approach_dist/2, turn_exit + exit_dir * approach_dist];
vel_exit = [exit_dir * V/2, [0; 0; 0]];
acc_exit = [nan(3,1), [0; 0; 0]];

% === COMBINE ===
wp = [wp_approach, wp_turn, wp_exit];
velBC = [vel_approach, vel_turn, vel_exit];
accBC = [acc_approach, acc_turn, acc_exit];

% === TIME ===
total_arc = sqrt(arc_length^2 + dh^2);
turn_time = total_arc / V;
t_approach = [0, approach_time/2];
t_turn = approach_time + linspace(0, turn_time, n_turn);
t_exit = t_turn(end) + [approach_time/2, approach_time];
tp = [t_approach, t_turn, t_exit];

info.name = sprintf('Climbing Turn %.0f deg', turn_angle);
info.total_time = tp(end);
info.yaw_mode = 'tangent';
info.yaw_params = struct();
info.description = sprintf('%.0f deg turn, R=%.1fm, climb=%.1fm, with approach/exit', turn_angle, R, dh);
end

%% ============== FIGURE EIGHT ==============

function [wp, tp, velBC, accBC, info] = generate_figure_eight(start, R, V)
%GENERATE_FIGURE_EIGHT Figure-8 with approach/exit
%
% Structure:
%   1. APPROACH: Hover -> accelerate North
%   2. FIGURE-8: Two connected circles (CW then CCW)
%   3. EXIT: Decelerate to hover

% Approach/exit parameters
a_accel = 2.0;
approach_dist = V^2 / (2 * a_accel);
approach_dist = max(approach_dist, 1.5);
approach_time = V / a_accel;
approach_time = max(approach_time, 1.0);

% === APPROACH ===
wp_approach = [start, start + [approach_dist/2; 0; 0]];
vel_approach = [[0; 0; 0], [V/2; 0; 0]];
acc_approach = [[0; 0; 0], nan(3,1)];

% === FIGURE-8 ===
fig8_entry = start + [approach_dist; 0; 0];

% First circle (clockwise)
center1_N = fig8_entry(1);
center1_E = fig8_entry(2) + R;
theta1 = deg2rad([0, 90, 180, 270, 360]);
N1 = center1_N + R * sin(theta1);
E1 = center1_E - R * cos(theta1);

% Second circle (counter-clockwise)
center2_N = fig8_entry(1);
center2_E = fig8_entry(2) - R;
theta2 = deg2rad([0, -90, -180, -270, -360]);
N2 = center2_N + R * sin(theta2);
E2 = center2_E + R * cos(theta2);

% Combine (skip duplicate)
N_fig8 = [N1, N2(2:end)];
E_fig8 = [E1, E2(2:end)];
D_fig8 = fig8_entry(3) * ones(size(N_fig8));
wp_fig8 = [N_fig8; E_fig8; D_fig8];
n_fig8 = size(wp_fig8, 2);

% Velocities
vel_fig8 = nan(3, n_fig8);
for i = 1:5
    th = theta1(i);
    vel_fig8(:, i) = [V * cos(th); V * sin(th); 0];
end
for i = 2:5
    th = theta2(i);
    vel_fig8(:, 4+i) = [V * cos(th); -V * sin(th); 0];
end

acc_fig8 = nan(3, n_fig8);

% === EXIT ===
fig8_exit = wp_fig8(:, end);
wp_exit = [fig8_exit + [approach_dist/2; 0; 0], fig8_exit + [approach_dist; 0; 0]];
vel_exit = [[V/2; 0; 0], [0; 0; 0]];
acc_exit = [nan(3,1), [0; 0; 0]];

% === COMBINE ===
wp = [wp_approach, wp_fig8, wp_exit];
velBC = [vel_approach, vel_fig8, vel_exit];
accBC = [acc_approach, acc_fig8, acc_exit];

% === TIME ===
circumference = 2 * pi * R;
circle_time = circumference / V;
fig8_time = 2 * circle_time;
t_approach = [0, approach_time/2];
t_fig8 = approach_time + linspace(0, fig8_time, n_fig8);
t_exit = t_fig8(end) + [approach_time/2, approach_time];
tp = [t_approach, t_fig8, t_exit];

info.name = 'Figure Eight';
info.total_time = tp(end);
info.yaw_mode = 'tangent';
info.yaw_params = struct();
info.description = sprintf('Figure-8, R=%.1fm at %.1fm/s, with approach/exit', R, V);
end

%% ============== HELPER FUNCTION ==============


%% ============== RECTANGULAR CIRCUIT ==============

function [wp, tp, velBC, accBC, info] = generate_hippodrome(start, L, W, V)
%GENERATE_HIPPODROME Hippodrome (Stadium) shaped circuit
%
% Shape: Two parallel straights connected by two 180-degree turns (semicircles).
% Orientation: Long axis aligned with North-South.
% Direction: Counter-Clockwise (CCW).
% Start: Middle of the East straight, heading North.
%
% Inputs:
%   L - Total length (North-South bounding box)
%   W - Total width (East-West bounding box) = Turn Diameter
%   V - Velocity

% Derived dimensions
R = W / 2;          % Turn radius
L_straight = L - W; % Length of the straight component

if L_straight < 0
    warning('Length must be greater than Width for a hippodrome. Clamping to Circle.');
    L_straight = 0;
    R = L / 2;
    W = L;
end

% Key relative points in "Center Frame" (Center at 0,0)
% X = North, Y = East
% We fly CCW.
% Start is on East edge (Y=+R), at middle (X=0). Heading North (+V_x).

% Define points relative to geometric center
% 1. Start (Mid-East)
p1 = [0; R; 0];
v1 = [V; 0; 0];

% 2. North-East Corner (End of East Straight)
p2 = [L_straight/2; R; 0];
v2 = [V; 0; 0];

% 3. North Apex (Top of Turn)
p3 = [L_straight/2 + R; 0; 0];
v3 = [0; -V; 0];

% 4. North-West Corner (Start of West Straight)
p4 = [L_straight/2; -R; 0];
v4 = [-V; 0; 0];

% 5. South-West Corner (End of West Straight)
p5 = [-L_straight/2; -R; 0];
v5 = [-V; 0; 0];

% 6. South Apex (Bottom of Turn)
p6 = [-L_straight/2 - R; 0; 0];
v6 = [0; V; 0];

% 7. South-East Corner (Start of East Straight)
p7 = [-L_straight/2; R; 0];
v7 = [V; 0; 0];

% 8. Back to Start
p8 = p1;
v8 = v1;

% Assemble
pts_rel = [p1, p2, p3, p4, p5, p6, p7, p8];
vels    = [v1, v2, v3, v4, v5, v6, v7, v8];
accs    = nan(3, 8); % Free acceleration (let solver handle smooth transitions)

% Transform to Global Frame
% Start pos matches p1 relative to center
% So Center = Start - p1
center_pos = start - p1;

wp = center_pos + pts_rel;
velBC = vels;
accBC = accs;

% === TIME ALLOCATION ===
% Calculate distances along the path
% Straights: L_straight
% Turns: pi * R
%
% Segments:
% 1->2: L_straight/2
% 2->3: pi*R / 2
% 3->4: pi*R / 2
% 4->5: L_straight
% 5->6: pi*R / 2
% 6->7: pi*R / 2
% 7->8: L_straight/2

seg_lens = [
    L_straight/2;
    pi*R/2;
    pi*R/2;
    L_straight;
    pi*R/2;
    pi*R/2;
    L_straight/2
    ];

path_len = [0; cumsum(seg_lens)];
tp = path_len' / V;

info.name = 'Hippodrome Circuit';
info.total_time = tp(end);
info.yaw_mode = 'tangent';
info.yaw_params = struct();
info.description = sprintf('Hippodrome %gx%gm, V=%.1fm/s', L, W, V);
end

%% ============== HELPER FUNCTION ==============

function val = get_param(params, name, default)
if isfield(params, name)
    val = params.(name);
else
    val = default;
end
end
