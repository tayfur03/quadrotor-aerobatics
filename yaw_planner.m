function [psi, psi_dot, psi_ddot] = yaw_planner(t, mode, state_ref, params)
%YAW_PLANNER Compute yaw reference for trajectory tracking.
%
% Yaw reference modes based on Tal & Karaman papers:
%   'constant'       - Fixed yaw angle (default: 0)
%   'constant_rate'  - Constant yaw rate rotation (for barrel roll, etc.)
%   'tangent'        - Align body X with horizontal velocity (2D)
%   'coordinated'    - Full 3D coordinated flight (body X aligned with velocity)
%   'knife_edge'     - Knife-edge flight mode (body Y vertical)
%   'poi'            - Point of interest tracking (look at a fixed point)
%
% Inputs:
%   t         - Current time [s]
%   mode      - Yaw mode string
%   state_ref - Reference state struct with fields:
%               .v  - Velocity vector [3x1]
%               .a  - Acceleration vector [3x1]
%   params    - Parameters struct, may contain:
%               .psi_constant    - Target yaw for 'constant' mode [rad]
%               .psi_rate        - Yaw rate for 'constant_rate' mode [rad/s]
%               .psi_initial     - Initial yaw for 'constant_rate' mode [rad]
%               .poi             - Point of interest [3x1] for 'poi' mode
%               .last_psi        - Previous yaw for continuity
%               .knife_edge_dir  - Direction for knife-edge: 'left' or 'right'
%
% Outputs:
%   psi       - Yaw angle [rad]
%   psi_dot   - Yaw rate [rad/s]
%   psi_ddot  - Yaw acceleration [rad/s^2]
%
% Reference: Tal et al. (2023) "Aerobatic Trajectory Generation for a VTOL
%            Fixed-Wing Aircraft Using Differential Flatness"

switch lower(mode)
    case "constant"
        % Fixed yaw angle
        if isfield(params, 'psi_constant')
            psi = params.psi_constant;
        else
            psi = 0;
        end
        psi_dot = 0;
        psi_ddot = 0;

    case "constant_rate"
        % Constant yaw rate - useful for barrel rolls and rolling maneuvers
        % From paper: "rolling/yawing motion where ψ changes at the same rate"
        if isfield(params, 'psi_rate')
            psi_rate = params.psi_rate;
        else
            psi_rate = pi/2;  % Default: 90 deg/s
        end

        if isfield(params, 'psi_initial')
            psi_initial = params.psi_initial;
        else
            psi_initial = 0;
        end

        psi = psi_initial + psi_rate * t;
        psi_dot = psi_rate;
        psi_ddot = 0;

    case {"tangent", "coordinated"}
        % Coordinated flight - align body X with velocity vector
        % From paper: "coordinated flight without yaw" and "3D coordinated flight"

        vx = state_ref.v(1);
        vy = state_ref.v(2);
        vz = state_ref.v(3);

        ax = state_ref.a(1);
        ay = state_ref.a(2);
        az = state_ref.a(3);

        v_sq = vx^2 + vy^2 + vz^2; % 3D Hız büyüklüğü

        % Hız eşiği: Çok yavaşsa (Hover) yön değiştirmeye çalışma
        if v_sq < 0.1
            if isfield(params, 'last_psi')
                psi = params.last_psi;
            else
                psi = 0;
            end
            psi_dot = 0;
            psi_ddot = 0;
        else
            %% KOORDİNELİ YAW HESABI (Doğal Takip)
            % Adım 1: İstenen Thrust Yönü (Body Z - b3)
            % a_ref = T_vec + g -> T_vec = a_ref - g
            % Drone Z ekseni Thrust'ın tersinedir.
            g_vec = [0;0;9.81];
            t_vec = [ax; ay; az] - g_vec;

            % Eğer serbest düşüşteyse (t_vec=0), Z eksenini aşağı kabul et
            if norm(t_vec) < 1e-3
                b3 = [0;0;1];
            else
                b3 = -t_vec / norm(t_vec);
            end

            % Adım 2: İstenen Gidiş Yönü (Velocity Unit Vector)
            vel_dir = [vx; vy; vz] / sqrt(v_sq);

            % Adım 3: Body Y Ekseni (b2)
            % Hız vektörü ile Thrust vektörünün oluşturduğu düzleme diktir.
            % Bu, drone'un "yatış" (bank) eksenidir.
            y_temp = cross(b3, vel_dir);

            % Singularity Check: Eğer hız vektörü ile Z ekseni paralelse (Tam dik dalış/tırmanış)
            if norm(y_temp) < 1e-3
                % Standart bir Y ekseni seç (Örn: Doğu)
                b2 = [0;1;0];
            else
                b2 = y_temp / norm(y_temp);
            end

            % Adım 4: Body X Ekseni (b1) - Doğal Burun Yönü
            % Z ve Y belliyse, X bunların cross product'ıdır.
            b1 = cross(b2, b3);

            % Adım 5: Rotasyon Matrisinden Yaw Çıkarma
            % R_des = [b1, b2, b3]
            % Yaw (Psi) = atan2(R(2,1), R(1,1))
            target_psi = atan2(b1(2), b1(1));

            % AÇI SÜREKLİLİĞİ (Unwrap)
            if isfield(params, 'last_psi')
                diff = target_psi - params.last_psi;
                while diff > pi,  diff = diff - 2*pi; end
                while diff < -pi, diff = diff + 2*pi; end
                psi = params.last_psi + diff;
            else
                psi = target_psi;
            end

            % TÜREVLER (Feedforward)
            v_hor_sq = vx^2 + vy^2;
            if v_hor_sq > 0.1
                psi_dot = (vx*ay - vy*ax) / v_hor_sq;

                % Analytical psi_ddot if Jerk is available
                if isfield(state_ref, 'j')
                    jx = state_ref.j(1);
                    jy = state_ref.j(2);

                    % N = vx*ay - vy*ax
                    % D = vx^2 + vy^2
                    % psi_ddot = (N_dot*D - N*D_dot) / D^2

                    N = vx*ay - vy*ax;
                    D = v_hor_sq;

                    % N_dot = vx*jy - vy*jx
                    N_dot = vx*jy - vy*jx;

                    % D_dot = 2*vx*ax + 2*vy*ay
                    D_dot = 2*vx*ax + 2*vy*ay;

                    psi_ddot = (N_dot*D - N*D_dot) / (D^2);
                else
                    psi_ddot = 0;
                end
            else
                psi_dot = 0;
                psi_ddot = 0;
            end
        end

    case "knife_edge"
        % Knife-edge flight mode - body Y axis vertical
        % From paper: "knife-edge flight" where ψ = π/2 or -π/2
        % The vehicle flies sideways with wingtip pointing up/down
        %
        % In knife-edge, the body Y-axis is aligned with the vertical (iz)
        % This means yaw is offset by ±90° from velocity direction

        vx = state_ref.v(1);
        vy = state_ref.v(2);
        v_horiz_sq = vx^2 + vy^2;

        % Get knife-edge direction (left or right bank)
        if isfield(params, 'knife_edge_dir') && strcmpi(params.knife_edge_dir, 'left')
            offset = -pi/2;  % Left knife-edge (roll left)
        else
            offset = pi/2;   % Right knife-edge (roll right)
        end

        if v_horiz_sq > 0.1
            % Yaw perpendicular to velocity direction
            vel_heading = atan2(vy, vx);
            target_psi = vel_heading + offset;
        else
            % At low speed, maintain current or default
            if isfield(params, 'last_psi')
                target_psi = params.last_psi;
            else
                target_psi = offset;
            end
        end

        % Angle continuity (unwrap)
        if isfield(params, 'last_psi')
            diff = target_psi - params.last_psi;
            while diff > pi,  diff = diff - 2*pi; end
            while diff < -pi, diff = diff + 2*pi; end
            psi = params.last_psi + diff;
        else
            psi = target_psi;
        end

        % Yaw rate derivative (simplified)
        ax = state_ref.a(1);
        ay = state_ref.a(2);
        if v_horiz_sq > 0.1
            psi_dot = (vx*ay - vy*ax) / v_horiz_sq;
        else
            psi_dot = 0;
        end
        psi_ddot = 0;

    case "poi"
        % Point of Interest tracking - look at a fixed point in space
        % Useful for cinematic shots or inspection tasks

        if ~isfield(params, 'poi') || ~isfield(state_ref, 'pos')
            psi = 0; psi_dot = 0; psi_ddot = 0;
            return;
        end

        poi = params.poi(:);
        pos = state_ref.pos(:);

        % Vector from drone to POI
        to_poi = poi - pos;
        dx = to_poi(1);
        dy = to_poi(2);

        dist_horiz_sq = dx^2 + dy^2;

        if dist_horiz_sq > 0.01  % More than 10cm horizontal distance
            target_psi = atan2(dy, dx);

            % Angle continuity
            if isfield(params, 'last_psi')
                diff = target_psi - params.last_psi;
                while diff > pi,  diff = diff - 2*pi; end
                while diff < -pi, diff = diff + 2*pi; end
                psi = params.last_psi + diff;
            else
                psi = target_psi;
            end

            % Yaw rate from velocity relative to POI
            if isfield(state_ref, 'v')
                vx = state_ref.v(1);
                vy = state_ref.v(2);
                % d/dt(atan2(dy,dx)) = (dx*vy_rel - dy*vx_rel) / (dx^2 + dy^2)
                psi_dot = (dx * (-vy) - dy * (-vx)) / dist_horiz_sq;
            else
                psi_dot = 0;
            end
        else
            if isfield(params, 'last_psi')
                psi = params.last_psi;
            else
                psi = 0;
            end
            psi_dot = 0;
        end
        psi_ddot = 0;

    case "velocity_aligned"
        % Full 3D velocity-aligned yaw (for high-speed aerobatics)
        % Body X always points in the direction of 3D velocity
        vx = state_ref.v(1);
        vy = state_ref.v(2);
        v_horiz_sq = vx^2 + vy^2;

        if v_horiz_sq > 0.1
            psi = atan2(vy, vx);
        else
            if isfield(params, 'last_psi')
                psi = params.last_psi;
            else
                psi = 0;
            end
        end
        psi_dot = 0;
        psi_ddot = 0;

    otherwise
        psi = 0; psi_dot = 0; psi_ddot = 0;
end
end