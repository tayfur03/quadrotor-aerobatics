function [segment_times, info] = optimize_time_allocation(waypoints, v_max, a_max, velBC, params)
% OPTIMIZE_TIME_ALLOCATION 3-Stage Hybrid Time Allocation
%
% Stage 1: A priori feasibility bounds (turn angle + kinematic limits)
% Stage 2: Mellinger closed-form snap-optimal times with aggressiveness k_T
% Stage 3: Single-pass post-check & constrained refinement
%
% This replaces the iterative inflate-and-retry approach with a
% mathematically grounded optimization that produces "predator-like"
% trajectories: maximum speed on straights, controlled deceleration
% through turns, aggressive acceleration out of them.
%
% Inputs:
%   waypoints - [3 x N] waypoint positions
%   v_max     - Maximum allowable velocity [m/s]
%   a_max     - Maximum allowable acceleration [m/s^2]
%   velBC     - [3 x N] Velocity boundary conditions (NaN = free)
%   params    - (optional) Struct with:
%       - aggressiveness: k_T parameter [0.1=gentle, 1=balanced, 10=predator]
%                         (default: 2.0)
%       - samples_per_seg: Samples per segment for post-check (default: 50)
%       - v_safety: Velocity safety factor for overshoot (default: 1.3)
%
% Outputs:
%   segment_times - [1 x (N-1)] Optimized time per segment
%   info          - Struct with diagnostics

    if nargin < 5, params = struct(); end

    k_T            = get_param(params, 'aggressiveness', 2.0);
    samples_per_seg = get_param(params, 'samples_per_seg', 50);
    v_safety       = get_param(params, 'v_safety', 1.3);  % Polynomial overshoot factor

    n_wp = size(waypoints, 2);
    n_seg = n_wp - 1;

    if n_seg < 1
        error('Need at least 2 waypoints');
    end

    % =====================================================================
    % STAGE 1: A PRIORI FEASIBILITY BOUNDS
    % Compute hard lower bounds on T_i from physics BEFORE solving any QP.
    % This ensures the optimizer starts in a feasible region.
    % =====================================================================

    distances = zeros(1, n_seg);
    for i = 1:n_seg
        distances(i) = norm(waypoints(:, i+1) - waypoints(:, i));
    end

    % 1a. Velocity bound: T >= d / v_max * safety_factor
    %     The safety factor accounts for polynomial overshoot above avg speed
    T_min_vel = distances / v_max * v_safety;

    % 1b. Acceleration bound: trapezoidal profile assumption
    %     Accelerate at a_max for half the segment, decelerate for the other half
    %     d = 0.5 * a_max * (T/2)^2 * 2 = 0.25 * a_max * T^2
    %     => T >= 2 * sqrt(d / a_max)
    T_min_acc = 2 * sqrt(distances / a_max);

    % 1c. Turn angle bound: sharp turns need deceleration
    %     At waypoint i, the turn angle theta_i determines the centripetal
    %     demand. The drone must slow down so that v^2/R <= a_max.
    T_min_turn = zeros(1, n_seg);

    turn_angles = zeros(1, n_wp);
    for i = 2:n_wp-1
        v_in = waypoints(:, i) - waypoints(:, i-1);
        v_out = waypoints(:, i+1) - waypoints(:, i);
        len_in = norm(v_in);
        len_out = norm(v_out);

        if len_in > 1e-6 && len_out > 1e-6
            cos_theta = dot(v_in, v_out) / (len_in * len_out);
            cos_theta = max(-1, min(1, cos_theta));  % Clamp for acos
            turn_angles(i) = acos(cos_theta);  % [0, pi], 0=straight, pi=U-turn
        end
    end

    % For each segment, the limiting turn angle is at the waypoint it
    % enters. Sharp turns require more time on the incoming segment.
    v_cruise = v_max * 0.7;  % Nominal cruise speed for turn calculations
    for i = 1:n_seg
        % The turn at the END of this segment (=waypoint i+1)
        if i < n_seg
            theta = turn_angles(i+1);
            if theta > 0.1  % More than ~6 degrees
                % Required deceleration to make the turn safely
                % v_turn = sqrt(a_max * R), approximate R from geometry
                % Simpler: slow down proportionally to turn sharpness
                v_turn = v_cruise * max(0.2, cos(theta/2));
                decel_time = abs(v_cruise - v_turn) / a_max;
                T_min_turn(i) = max(T_min_turn(i), distances(i) / v_turn);
            end
        end

        % The turn at the START of this segment (=waypoint i)
        if i > 1
            theta = turn_angles(i);
            if theta > 0.1
                v_turn = v_cruise * max(0.2, cos(theta/2));
                accel_time = abs(v_cruise - v_turn) / a_max;
                T_min_turn(i) = max(T_min_turn(i), distances(i) / v_turn);
            end
        end
    end

    % Combined lower bound
    T_lower = max(T_min_vel, max(T_min_acc, T_min_turn));
    T_lower = max(T_lower, 0.5);  % Absolute minimum 0.5s

    % =====================================================================
    % STAGE 2: MELLINGER CLOSED-FORM SNAP-OPTIMAL TIMES
    % For a 7th-order polynomial, snap cost ~ c_i / T_i^7
    % Total cost J = sum(c_i/T_i^7) + k_T * sum(T_i)
    % Optimal: T_i* = (7*c_i / k_T)^(1/8)
    %
    % We approximate c_i from segment geometry (distance + turn angles).
    % =====================================================================

    % Estimate snap cost coefficients c_i
    % Straight segments with small deflection: c ~ d^4 (distance dominated)
    % Sharp turns: c increases with turn sharpness (more snap needed)
    c = zeros(1, n_seg);
    for i = 1:n_seg
        d = distances(i);

        % Base snap coefficient: proportional to d^4
        c(i) = d^4;

        % Amplify c for sharp turns (entering or exiting this segment)
        turn_factor = 1.0;
        if i > 1 && turn_angles(i) > 0.1
            turn_factor = turn_factor * (1 + 2 * (turn_angles(i) / pi)^2);
        end
        if i < n_seg && turn_angles(i+1) > 0.1
            turn_factor = turn_factor * (1 + 2 * (turn_angles(i+1) / pi)^2);
        end
        c(i) = c(i) * turn_factor;
    end

    % Closed-form optimal: T_i* = (7 * c_i / k_T)^(1/8)
    T_optimal = (7 * c / k_T) .^ (1/8);

    % Apply feasibility floor
    segment_times = max(T_optimal, T_lower);

    fprintf('Stage 1-2: T_lower=[%.1f–%.1fs], T_optimal=[%.1f–%.1fs], k_T=%.1f\n', ...
        min(T_lower), max(T_lower), min(T_optimal), max(T_optimal), k_T);

    % =====================================================================
    % STAGE 3: SINGLE-PASS POST-CHECK & REFINEMENT
    % Solve the actual min-snap problem once with the optimized times,
    % check v/a constraints, and do ONE targeted correction pass.
    % =====================================================================

    accBC = nan(3, n_wp);
    accBC(:, 1) = [0; 0; 0];
    accBC(:, end) = [0; 0; 0];

    timePoints = [0, cumsum(segment_times)];
    numSamples = n_seg * samples_per_seg;

    info.converged = false;
    info.iterations = 0;
    info.max_vel_ratio = NaN;
    info.max_acc_ratio = NaN;
    info.turn_angles_deg = turn_angles * 180/pi;
    info.T_lower = T_lower;
    info.T_optimal = T_optimal;

    % Two passes max: first solve, then one refinement if needed
    for pass = 1:2
        timePoints = [0, cumsum(segment_times)];

        try
            [~, vel, acc, ~, ~, ~, ~, tSamples] = minsnappolytraj(...
                waypoints, timePoints, numSamples, ...
                'VelocityBoundaryCondition', velBC, ...
                'AccelerationBoundaryCondition', accBC);
        catch ME
            warning('optimize_time_allocation: minsnappolytraj failed (pass %d): %s', pass, ME.message);
            % Use current times as best effort
            break;
        end

        % Check for NaN/Inf
        if any(isnan(vel(:))) || any(isinf(vel(:)))
            warning('optimize_time_allocation: NaN/Inf in solution (pass %d)', pass);
            % Increase all times by 20% and try once more
            if pass == 1
                segment_times = segment_times * 1.2;
                continue;
            end
            break;
        end

        % Per-segment violation analysis
        vel_norms = vecnorm(vel, 2, 1);
        acc_norms = vecnorm(acc, 2, 1);

        violations = ones(1, n_seg);
        seg_max_v = zeros(1, n_seg);
        seg_max_a = zeros(1, n_seg);

        for s = 1:n_seg
            if s < n_seg
                mask = tSamples >= timePoints(s) & tSamples < timePoints(s+1);
            else
                mask = tSamples >= timePoints(s) & tSamples <= timePoints(s+1);
            end

            if ~any(mask), continue; end

            seg_max_v(s) = max(vel_norms(mask));
            seg_max_a(s) = max(acc_norms(mask));

            scale_v = seg_max_v(s) / v_max;
            scale_a = sqrt(seg_max_a(s) / a_max);

            violations(s) = max(scale_v, scale_a);
        end

        info.max_vel_ratio = max(vel_norms) / v_max;
        info.max_acc_ratio = max(acc_norms) / a_max;
        info.iterations = pass;

        % Check convergence (5% tolerance)
        if max(violations) < 1.05
            info.converged = true;
            fprintf('Time allocation converged (pass %d): T=%.1fs, v_ratio=%.2f, a_ratio=%.2f\n', ...
                pass, sum(segment_times), info.max_vel_ratio, info.max_acc_ratio);
            break;
        end

        if pass == 1
            % Targeted refinement: only inflate violating segments
            n_fixed = 0;
            for s = 1:n_seg
                if violations(s) > 1.0
                    segment_times(s) = segment_times(s) * violations(s) * 1.05;
                    n_fixed = n_fixed + 1;
                end
            end
            fprintf('Stage 3 refinement: %d/%d segments inflated, re-solving...\n', n_fixed, n_seg);
        else
            fprintf('Time allocation (pass 2): T=%.1fs, v_ratio=%.2f, a_ratio=%.2f\n', ...
                sum(segment_times), info.max_vel_ratio, info.max_acc_ratio);
        end
    end

    % Final safety
    segment_times = max(segment_times, 0.1);

    info.total_time = sum(segment_times);
    info.distances = distances;
    info.seg_max_v = seg_max_v;
    info.seg_max_a = seg_max_a;

    fprintf('Final: %.1fs total, %.1f m/s max vel, %.1f m/s² max acc (aggressiveness=%.1f)\n', ...
        info.total_time, max(seg_max_v), max(seg_max_a), k_T);
end

function val = get_param(params, name, default)
    if isfield(params, name)
        val = params.(name);
    else
        val = default;
    end
end
