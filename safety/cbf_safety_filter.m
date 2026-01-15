classdef cbf_safety_filter < handle
%CBF_SAFETY_FILTER Control Barrier Function safety filter for UAV
%
% Implements safety constraints using Control Barrier Functions (CBFs).
% Modifies control inputs to ensure safety while minimally deviating
% from the nominal controller.
%
% Safety Constraints:
%   1. Terrain Clearance: h_terrain(x) = alt - terrain_height - h_min >= 0
%   2. Radar Exposure: h_radar(x) = P_max - P_det(x) >= 0 (optional)
%   3. Velocity Limits: h_vel(x) = v_max^2 - ||v||^2 >= 0
%
% CBF Condition:
%   For safety function h(x), require: dh/dt + alpha * h(x) >= 0
%   This ensures h(x) >= 0 is forward invariant.
%
% The filter solves a QP at each timestep:
%   min  ||u - u_nom||^2
%   s.t. CBF constraints
%
% Properties:
%   terrain      - terrain_map object
%   threat       - threat_map object (optional)
%   h_min        - Minimum terrain clearance [m]
%   P_max        - Maximum allowed detection probability
%   v_max        - Maximum velocity [m/s]
%   alpha        - CBF class-K function gain (larger = more aggressive)
%
% Methods:
%   filter(state, u_nom) - Apply safety filter to nominal control
%   get_barrier_values(state) - Get current barrier function values
%   is_safe(state) - Check if current state satisfies constraints
%
% Example:
%   cbf = cbf_safety_filter(terrain_map, threat_map);
%   cbf.h_min = 30;      % 30m clearance
%   cbf.alpha = 1.0;     % Moderate aggressiveness
%
%   % In control loop:
%   u_safe = cbf.filter(state, u_nominal);
%
% Author: Quadrotor Terrain Following Project

    properties
        terrain         % terrain_map object
        threat          % threat_map object (optional)
        h_min           % Minimum terrain clearance [m]
        P_max           % Maximum detection probability
        v_max           % Maximum velocity [m/s]
        a_max           % Maximum acceleration [m/s^2]
        alpha_terrain   % CBF gain for terrain constraint
        alpha_radar     % CBF gain for radar constraint
        alpha_vel       % CBF gain for velocity constraint
        enable_terrain  % Enable terrain CBF
        enable_radar    % Enable radar CBF
        enable_velocity % Enable velocity CBF
        relaxation      % Slack variable weight for soft constraints
    end

    properties (Access = private)
        last_filter_info  % Store last filter computation info
    end

    methods
        function obj = cbf_safety_filter(terrain_map_obj, threat_map_obj, params)
            %CBF_SAFETY_FILTER Constructor
            %   cbf = cbf_safety_filter(terrain)
            %   cbf = cbf_safety_filter(terrain, threat)
            %   cbf = cbf_safety_filter(terrain, threat, params)

            obj.terrain = terrain_map_obj;

            if nargin >= 2
                obj.threat = threat_map_obj;
            else
                obj.threat = [];
            end

            if nargin < 3
                params = struct();
            end

            % Default parameters
            obj.h_min = get_param(params, 'h_min', 20);
            obj.P_max = get_param(params, 'P_max', 0.3);
            obj.v_max = get_param(params, 'v_max', 20);
            obj.a_max = get_param(params, 'a_max', 10);
            obj.alpha_terrain = get_param(params, 'alpha_terrain', 1.0);
            obj.alpha_radar = get_param(params, 'alpha_radar', 0.5);
            obj.alpha_vel = get_param(params, 'alpha_vel', 2.0);
            obj.enable_terrain = get_param(params, 'enable_terrain', true);
            obj.enable_radar = get_param(params, 'enable_radar', false);
            obj.enable_velocity = get_param(params, 'enable_velocity', true);
            obj.relaxation = get_param(params, 'relaxation', 1000);
        end

        function [u_safe, info] = filter(obj, state, u_nom)
            %FILTER Apply CBF safety filter to nominal control
            %   [u_safe, info] = filter(state, u_nom)
            %
            %   state - Struct with:
            %           .x   - Position [N; E; D] in NED
            %           .v   - Velocity [vN; vE; vD] in NED
            %   u_nom - Nominal control (acceleration) [aN; aE; aD]
            %
            %   u_safe - Safety-filtered control
            %   info   - Struct with filter details

            x = state.x(:);
            v = state.v(:);
            u_nom = u_nom(:);

            % Get current barrier values
            [h_terrain, grad_h_terrain] = obj.terrain_barrier(x, v);
            [h_radar, grad_h_radar] = obj.radar_barrier(x, v);
            [h_vel, grad_h_vel] = obj.velocity_barrier(v);

            % Build QP
            % Decision variables: [u (3); slack (1)]
            % Cost: ||u - u_nom||^2 + relaxation * slack^2
            % Constraints: A * [u; slack] <= b (CBF conditions)

            H = blkdiag(2*eye(3), 2*obj.relaxation);
            f = [-2*u_nom; 0];

            A = [];
            b = [];

            % Terrain CBF constraint
            % h_dot + alpha * h >= 0
            % grad_h' * v + grad_h' * u * dt + alpha * h >= -slack
            % For position barrier: h(x) = alt - terrain - h_min
            % h_dot = -vD - d(terrain)/dt ≈ -vD (for slowly varying terrain)
            % We need h_dot + alpha*h >= 0, which gives constraint on vD_dot = aD

            if obj.enable_terrain && h_terrain < 100  % Only active when close
                % Lie derivative: L_f h = grad_h' * f(x) where f(x) = v
                % L_g h = grad_h (for single integrator on velocity)
                % CBF: L_f h + L_g h * u + alpha * h >= 0
                %      grad_h' * v + grad_h' * u + alpha * h >= 0
                %      -grad_h' * u <= grad_h' * v + alpha * h

                A_terrain = [-grad_h_terrain', -1];  % Include slack
                b_terrain = grad_h_terrain' * v + obj.alpha_terrain * h_terrain;

                A = [A; A_terrain];
                b = [b; b_terrain];
            end

            % Radar CBF constraint
            if obj.enable_radar && ~isempty(obj.threat) && h_radar < 0.5
                % Similar structure: -grad_h' * u <= grad_h' * v + alpha * h
                A_radar = [-grad_h_radar', -1];
                b_radar = grad_h_radar' * v + obj.alpha_radar * h_radar;

                A = [A; A_radar];
                b = [b; b_radar];
            end

            % Velocity CBF constraint
            if obj.enable_velocity && h_vel < obj.v_max^2 * 0.5
                % h_vel = v_max^2 - ||v||^2
                % h_dot = -2 * v' * v_dot = -2 * v' * u
                % CBF: -2*v'*u + alpha * h >= 0
                %      2*v'*u <= alpha * h

                A_vel = [2*v', -1];
                b_vel = obj.alpha_vel * h_vel;

                A = [A; A_vel];
                b = [b; b_vel];
            end

            % Acceleration limits (hard constraint, no slack)
            A_acc = [eye(3), zeros(3,1); -eye(3), zeros(3,1)];
            b_acc = [obj.a_max; obj.a_max; obj.a_max; obj.a_max; obj.a_max; obj.a_max];
            A = [A; A_acc];
            b = [b; b_acc];

            % Slack non-negativity
            A = [A; 0, 0, 0, -1];
            b = [b; 0];

            % Solve QP
            options = optimoptions('quadprog', 'Display', 'off');

            if isempty(A)
                % No constraints active, use nominal
                u_safe = u_nom;
                slack = 0;
            else
                try
                    x_opt = quadprog(H, f, A, b, [], [], [], [], [], options);
                    if isempty(x_opt)
                        % QP infeasible, use projected nominal
                        u_safe = obj.project_to_safe(u_nom, state);
                        slack = inf;
                    else
                        u_safe = x_opt(1:3);
                        slack = x_opt(4);
                    end
                catch
                    % Solver failed, use projected nominal
                    u_safe = obj.project_to_safe(u_nom, state);
                    slack = inf;
                end
            end

            % Build info struct
            info.h_terrain = h_terrain;
            info.h_radar = h_radar;
            info.h_vel = h_vel;
            info.u_nom = u_nom;
            info.u_safe = u_safe;
            info.modification = norm(u_safe - u_nom);
            info.slack = slack;
            info.active_terrain = obj.enable_terrain && h_terrain < 100;
            info.active_radar = obj.enable_radar && h_radar < 0.5;
            info.active_vel = obj.enable_velocity && h_vel < obj.v_max^2 * 0.5;

            obj.last_filter_info = info;
        end

        function [h, grad_h] = terrain_barrier(obj, x, v)
            %TERRAIN_BARRIER Compute terrain clearance barrier function
            %   h = altitude - terrain_height - h_min
            %   Positive when safe, negative when violated

            N = x(1);
            E = x(2);
            D = x(3);

            alt = -D;  % Altitude (positive up)
            terrain_h = obj.terrain.get_height(N, E);

            h = alt - terrain_h - obj.h_min;

            % Gradient with respect to position
            % dh/dN = -d(terrain_h)/dN
            % dh/dE = -d(terrain_h)/dE
            % dh/dD = -1 (since alt = -D)

            delta = obj.terrain.resolution / 2;
            dterrain_dN = (obj.terrain.get_height(N + delta, E) - ...
                           obj.terrain.get_height(N - delta, E)) / (2*delta);
            dterrain_dE = (obj.terrain.get_height(N, E + delta) - ...
                           obj.terrain.get_height(N, E - delta)) / (2*delta);

            grad_h = [-dterrain_dN; -dterrain_dE; -1];
        end

        function [h, grad_h] = radar_barrier(obj, x, v)
            %RADAR_BARRIER Compute radar exposure barrier function
            %   h = P_max - P_det(x)
            %   Positive when safe (low detection), negative when dangerous

            if isempty(obj.threat) || ~obj.threat.computed
                h = 1;  % Always safe if no threat map
                grad_h = zeros(3, 1);
                return;
            end

            N = x(1);
            E = x(2);
            alt = -x(3);

            P_det = obj.threat.get_risk(N, E, alt);
            h = obj.P_max - P_det;

            % Numerical gradient
            delta = 10;  % 10m for gradient computation
            dP_dN = (obj.threat.get_risk(N + delta, E, alt) - ...
                     obj.threat.get_risk(N - delta, E, alt)) / (2*delta);
            dP_dE = (obj.threat.get_risk(N, E + delta, alt) - ...
                     obj.threat.get_risk(N, E - delta, alt)) / (2*delta);
            dP_dalt = (obj.threat.get_risk(N, E, alt + delta) - ...
                       obj.threat.get_risk(N, E, alt - delta)) / (2*delta);

            % grad_h = -grad_P, and convert altitude gradient to D
            grad_h = [-dP_dN; -dP_dE; dP_dalt];  % Note: dalt/dD = -1
        end

        function [h, grad_h] = velocity_barrier(obj, v)
            %VELOCITY_BARRIER Compute velocity limit barrier function
            %   h = v_max^2 - ||v||^2

            h = obj.v_max^2 - norm(v)^2;
            grad_h = -2 * v;
        end

        function [h_terrain, h_radar, h_vel] = get_barrier_values(obj, state)
            %GET_BARRIER_VALUES Get current barrier function values

            x = state.x(:);
            v = state.v(:);

            h_terrain = obj.terrain_barrier(x, v);
            h_radar = obj.radar_barrier(x, v);
            h_vel = obj.velocity_barrier(v);
        end

        function safe = is_safe(obj, state)
            %IS_SAFE Check if current state satisfies all constraints

            [h_terrain, h_radar, h_vel] = obj.get_barrier_values(state);

            safe = true;
            if obj.enable_terrain && h_terrain < 0
                safe = false;
            end
            if obj.enable_radar && h_radar < 0
                safe = false;
            end
            if obj.enable_velocity && h_vel < 0
                safe = false;
            end
        end

        function u_safe = project_to_safe(obj, u_nom, state)
            %PROJECT_TO_SAFE Simple fallback: project to safe direction

            x = state.x(:);
            v = state.v(:);

            u_safe = u_nom;

            % If below minimum clearance, command climb
            [h_terrain, grad_h] = obj.terrain_barrier(x, v);
            if h_terrain < 0
                % Emergency climb: set vertical acceleration upward
                climb_acc = min(obj.a_max, -obj.alpha_terrain * h_terrain);
                u_safe(3) = -climb_acc;  % Negative D acceleration = climb
            end

            % Limit velocity if too fast
            if norm(v) > obj.v_max
                % Decelerate in velocity direction
                v_dir = v / norm(v);
                decel = min(obj.a_max, obj.alpha_vel * (norm(v) - obj.v_max));
                u_safe = u_safe - decel * v_dir;
            end

            % Clip to acceleration limits
            u_safe = max(-obj.a_max, min(obj.a_max, u_safe));
        end

        function fig = plot_diagnostics(obj, log, fig_handle)
            %PLOT_DIAGNOSTICS Plot CBF diagnostic information
            %   plot_diagnostics(log) - log contains t, h_terrain, h_radar, etc.

            if nargin < 3
                fig = figure('Name', 'CBF Diagnostics', 'Position', [100, 100, 1000, 800]);
            else
                fig = figure(fig_handle);
            end

            t = log.t;

            subplot(3, 2, 1);
            plot(t, log.h_terrain, 'b-', 'LineWidth', 1.5);
            hold on;
            plot(t, zeros(size(t)), 'r--', 'LineWidth', 1);
            xlabel('Time [s]');
            ylabel('h_{terrain} [m]');
            title('Terrain Clearance Barrier');
            legend('h(x)', 'Safe boundary');
            grid on;

            subplot(3, 2, 2);
            if isfield(log, 'h_radar')
                plot(t, log.h_radar, 'b-', 'LineWidth', 1.5);
                hold on;
                plot(t, zeros(size(t)), 'r--', 'LineWidth', 1);
            end
            xlabel('Time [s]');
            ylabel('h_{radar}');
            title('Radar Exposure Barrier');
            grid on;

            subplot(3, 2, 3);
            plot(t, log.h_vel, 'b-', 'LineWidth', 1.5);
            hold on;
            plot(t, zeros(size(t)), 'r--', 'LineWidth', 1);
            xlabel('Time [s]');
            ylabel('h_{vel} [m^2/s^2]');
            title('Velocity Limit Barrier');
            grid on;

            subplot(3, 2, 4);
            if isfield(log, 'modification')
                plot(t, log.modification, 'b-', 'LineWidth', 1.5);
            end
            xlabel('Time [s]');
            ylabel('||u_{safe} - u_{nom}||');
            title('Control Modification Magnitude');
            grid on;

            subplot(3, 2, 5);
            if isfield(log, 'slack')
                semilogy(t, log.slack + 1e-10, 'b-', 'LineWidth', 1.5);
            end
            xlabel('Time [s]');
            ylabel('Slack Variable');
            title('Constraint Relaxation');
            grid on;

            subplot(3, 2, 6);
            if isfield(log, 'alt') && isfield(log, 'terrain_h')
                plot(t, log.alt, 'b-', 'LineWidth', 1.5);
                hold on;
                plot(t, log.terrain_h, 'k-', 'LineWidth', 1.5);
                plot(t, log.terrain_h + obj.h_min, 'r--', 'LineWidth', 1);
            end
            xlabel('Time [s]');
            ylabel('Altitude [m]');
            title('Altitude vs Terrain');
            legend('UAV Alt', 'Terrain', 'Min Clearance');
            grid on;
        end
    end
end

function val = get_param(params, name, default)
    if isfield(params, name)
        val = params.(name);
    else
        val = default;
    end
end
