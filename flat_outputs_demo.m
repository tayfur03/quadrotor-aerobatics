function [x_ref, v_ref, a_ref, j_ref, s_ref, psi_ref, psi_dot_ref, psi_ddot_ref] = flat_outputs_demo(t, params)
%FLAT_OUTPUTS_DEMO  Smooth 3D trajectory (single param), with derivatives.
% NED frame: x = [N; E; D], D is positive downward.
%
% Uses yaw_planner() for yaw reference generation.
%
% Inputs:
%   t      - Current time [s]
%   params - Struct with fields:
%            .shape         - "eight", "heart", "roulette", or "barrel_roll"
%            .yaw           - Yaw mode: "constant", "tangent", "coordinated",
%                             "constant_rate", "knife_edge", "poi"
%            .psi_ref_const - (optional) Constant yaw angle [rad]
%            .yaw_params    - (optional) Additional yaw params for yaw_planner
%            .traj          - (optional) Trajectory overrides
%            .roulette      - (optional) Roulette curve parameters
%
% Outputs:
%   x_ref, v_ref, a_ref, j_ref, s_ref : [3x1] each (position to snap)
%   psi_ref, psi_dot_ref, psi_ddot_ref : yaw references

    % =========================
    % User-tunable trajectory params
    % =========================
    Tf   = params.Tf;              % total duration [s]
    A    = 5;                 % horizontal size [m] (Roulette modunda scale faktor olarak kullanilir)
    D0   = 0.5;               % base Down offset [m]
    Dz   = 1;                 % vertical oscillation amplitude [m]
    shape = params.shape ;    % "eight", "heart", "roulette", or "barrel_roll"
    yaw_mode = params.yaw ;   % yaw mode for yaw_planner
    cycles = 1.0;             % number of loops within Tf (can be 1, 1.5, 2, ...)

    % Optional trajectory overrides:
    % params.traj.Tf, params.traj.A, params.traj.D0, params.traj.Dz, params.traj.cycles
    if isfield(params, 'traj')
        if isfield(params.traj, 'Tf'), Tf = params.traj.Tf; end
        if isfield(params.traj, 'A'), A = params.traj.A; end
        if isfield(params.traj, 'D0'), D0 = params.traj.D0; end
        if isfield(params.traj, 'Dz'), Dz = params.traj.Dz; end
        if isfield(params.traj, 'cycles'), cycles = params.traj.cycles; end
    end

    % Optional roulette overrides (paper-accurate defaults).
    % params.roulette.scale, params.roulette.use_paper_time, params.roulette.r1..r6, params.roulette.k1..k3
    roulette = struct();
    if isfield(params, 'roulette')
        roulette = params.roulette;
    end
    use_paper_time = true;
    if isfield(roulette, 'use_paper_time')
        use_paper_time = logical(roulette.use_paper_time);
    end
    A_roulette = A;
    if isfield(roulette, 'scale')
        A_roulette = roulette.scale;
    elseif shape == "roulette" && use_paper_time
        A_roulette = 1.0;
    end
    % Persistent yaw for unwrap continuity across calls (used by yaw_planner)
    persistent last_psi_persist;
    if isempty(last_psi_persist)
        if isfield(params, 'psi_ref_const')
            last_psi_persist = params.psi_ref_const;
        elseif isfield(params, 'yaw_params') && isfield(params.yaw_params, 'psi_initial')
            last_psi_persist = params.yaw_params.psi_initial;
        else
            last_psi_persist = 0.0;
        end
    end
    if isfield(params, 'reset_yaw') && params.reset_yaw
        if isfield(params, 'psi_ref_const')
            last_psi_persist = params.psi_ref_const;
        elseif isfield(params, 'yaw_params') && isfield(params.yaw_params, 'psi_initial')
            last_psi_persist = params.yaw_params.psi_initial;
        else
            last_psi_persist = 0.0;
        end
    end
    if t <= 0
        if isfield(params, 'psi_ref_const')
            last_psi_persist = params.psi_ref_const;
        elseif isfield(params, 'yaw_params') && isfield(params.yaw_params, 'psi_initial')
            last_psi_persist = params.yaw_params.psi_initial;
        else
            last_psi_persist = 0.0;
        end
    end

    % =========================
    % Smooth time scaling s(t) in [0,1] with C^3 continuity
    % =========================
    tau = min(max(t / Tf, 0), 1);            % normalized time in [0,1]
    [s, s1, s2, s3, s4] = s_curve_c3(tau);   % s and derivatives wrt tau
    
    % Convert derivatives wrt time t:
    dtau = 1 / Tf;
    s_dot   = s1 * dtau;
    s_ddot  = s2 * dtau^2;
    s_jerk  = s3 * dtau^3;
    s_snap  = s4 * dtau^4;

    % =========================
    % Parametric angle theta(s)
    % =========================
    w = 2*pi*cycles;      % total angle over s in [0,1]
    if shape == "roulette" && use_paper_time
        th      = t;
        th_dot  = 1.0;
        th_ddot = 0.0;
        th_jerk = 0.0;
        th_snap = 0.0;
    else
        th      = w * s;
        th_dot  = w * s_dot;
        th_ddot = w * s_ddot;
        th_jerk = w * s_jerk;
        th_snap = w * s_snap;
    end

    % =========================
    % Define planar curve r(th) = [N; E] and its derivatives wrt theta
    % =========================
    switch shape
        case "eight"
            % Lemniscate-like figure-8
            N_th   = A * sin(th);
            E_th   = A * sin(th) * cos(th);

            dN_dth  = A * cos(th);
            dE_dth  = A * (cos(th)^2 - sin(th)^2);

            d2N_dth2 = -A * sin(th);
            d2E_dth2 = -4*A * sin(th) * cos(th);

            d3N_dth3 = -A * cos(th);
            d3E_dth3 = -4*A * (cos(th)^2 - sin(th)^2);

            d4N_dth4 = A * sin(th);
            d4E_dth4 = 16*A*cos(th)*sin(th);

        case "heart"
            % Classic heart curve
            N_th = A * (sin(th)^3);
            E_th = A * (13*cos(th) - 5*cos(2*th) - 2*cos(3*th) - cos(4*th)) / 16;

            % Derivatives wrt th 
            dN_dth  = A * 3*sin(th)^2*cos(th);
            dE_dth  = A * (-13*sin(th) + 10*sin(2*th) + 6*sin(3*th) + 4*sin(4*th)) / 16;

            d2N_dth2 = A * 3*(2*sin(th)*cos(th)^2 - sin(th)^3);
            d2E_dth2 = A * (-13*cos(th) + 20*cos(2*th) + 18*cos(3*th) + 16*cos(4*th)) / 16;

            d3N_dth3 = A * 3*(2*cos(th)^3 - 7*sin(th)^2*cos(th));
            d3E_dth3 = A * (13*sin(th) - 40*sin(2*th) - 54*sin(3*th) - 64*sin(4*th)) / 16;

            d4N_dth4 = A * 3*(6*cos(th)^2*sin(th) - 14*sin(th)*cos(th)^2);
            d4E_dth4 = A * (13*cos(th) - 80*cos(2*th) - 162*cos(3*th) - 256*cos(4*th)) / 16;

        case "roulette"
            % Tal & Karaman (2021) Roulette Curve
            % Makaledeki sabitler:
            r1=6; r2=1.8; r3=0.6; r4=-2.25; r5=-0.3; r6=-0.45;
            % Makaledeki frekanslar (k1=0.28, k2=2.8, k3=1.4)
            % Not: Bizim kodumuzda 'th' parametresi 'w' ile ölçekleniyor.
            % Makale orijinalinde doğrudan 't' kullanıyor.
            % 'A' parametresini genel ölçekleme (Scale) için kullanıyoruz.
            % Orjinal boyut için A=1 verilmelidir.

            k1=0.28; k2=2.8; k3=1.4;
            if isfield(roulette, 'r1'), r1 = roulette.r1; end
            if isfield(roulette, 'r2'), r2 = roulette.r2; end
            if isfield(roulette, 'r3'), r3 = roulette.r3; end
            if isfield(roulette, 'r4'), r4 = roulette.r4; end
            if isfield(roulette, 'r5'), r5 = roulette.r5; end
            if isfield(roulette, 'r6'), r6 = roulette.r6; end
            if isfield(roulette, 'k1'), k1 = roulette.k1; end
            if isfield(roulette, 'k2'), k2 = roulette.k2; end
            if isfield(roulette, 'k3'), k3 = roulette.k3; end

            % Parametrik hesaplamalar (th üzerinden)
            c1 = cos(k1*th); s1 = sin(k1*th);
            c2 = cos(k2*th); s2 = sin(k2*th);
            c3 = cos(k3*th); s3 = sin(k3*th);

            % Position
            N_th = A_roulette * (r1*c1 + r2*c2 + r3*s3);
            E_th = A_roulette * (r4*s1 + r5*s2 + r6*c3);

            % 1st Derivative (d/dth)
            dN_dth = A_roulette * (-r1*k1*s1 - r2*k2*s2 + r3*k3*c3);
            dE_dth = A_roulette * ( r4*k1*c1 + r5*k2*c2 - r6*k3*s3);

            % 2nd Derivative
            d2N_dth2 = A_roulette * (-r1*k1^2*c1 - r2*k2^2*c2 - r3*k3^2*s3);
            d2E_dth2 = A_roulette * (-r4*k1^2*s1 - r5*k2^2*s2 - r6*k3^2*c3);

            % 3rd Derivative
            d3N_dth3 = A_roulette * ( r1*k1^3*s1 + r2*k2^3*s2 - r3*k3^3*c3);
            d3E_dth3 = A_roulette * (-r4*k1^3*c1 - r5*k2^3*c2 + r6*k3^3*s3);

            % 4th Derivative
            d4N_dth4 = A_roulette * ( r1*k1^4*c1 + r2*k2^4*c2 + r3*k3^4*s3);
            d4E_dth4 = A_roulette * ( r4*k1^4*s1 + r5*k2^4*s2 + r6*k3^4*c3);
        case "barrel_roll"
            % Needs cruise flight before that manuever
            A = 1;
            N_th = 2*th;
            E_th = -1 + A*cos(th);
            D_th = - A*sin(th);

            % Derivatives wrt th 
            dN_dth = 5;
            dE_dth = -A*sin(th);
            dD_dth = -A*cos(th);

            d2N_dth2 = 0;
            d2E_dth2 = -A*cos(th);
            d2D_dth2 = A*sin(th);

            d3N_dth3 = 0;
            d3E_dth3 = A*sin(th);
            d3D_dth3 = A*cos(th);

            d4N_dth4 = 0;
            d4E_dth4 = A*cos(th);
            d4D_dth4 = -A*sin(th);

        otherwise
            error("Unknown shape. Use ""eight"", ""heart"" or ""roulette"".");
    end

    % =========================
    % Vertical (Down) profile
    % =========================
    if ~strcmp(shape, "barrel_roll")
        D_th = D0 + Dz * (1 - cos(th)) / 2;
        dD_dth   = Dz * (sin(th)) / 2;
        d2D_dth2 = Dz * (cos(th)) / 2;
        d3D_dth3 = -Dz * (sin(th)) / 2;
        d4D_dth4 = -Dz * (cos(th)) / 2;
    end

    % =========================
    % Chain rule: x(th(t)) -> Snap Calculation
    % =========================
    N      = N_th;
    E      = E_th;
    D      = D_th;

    N_dot  = dN_dth*th_dot;
    E_dot  = dE_dth*th_dot;
    D_dot  = dD_dth*th_dot;

    N_ddot = d2N_dth2*(th_dot^2) + dN_dth*th_ddot;
    E_ddot = d2E_dth2*(th_dot^2) + dE_dth*th_ddot;
    D_ddot = d2D_dth2*(th_dot^2) + dD_dth*th_ddot;

    N_jerk = d3N_dth3*(th_dot^3) + 3*d2N_dth2*th_dot*th_ddot + dN_dth*th_jerk;
    E_jerk = d3E_dth3*(th_dot^3) + 3*d2E_dth2*th_dot*th_ddot + dE_dth*th_jerk;
    D_jerk = d3D_dth3*(th_dot^3) + 3*d2D_dth2*th_dot*th_ddot + dD_dth*th_jerk;

    % SNAP (4. Türev) Zincir Kuralı Düzeltildi:
    % x(4) = x'''' th'^4 + 6 x''' th'² th'' + 3 x'' th''² + 4 x'' th' th''' + x' th(4)
    % Not: Orjinal kodda 3*x''*th''^2 terimi eksikti, buraya eklendi.
    
    N_snap = d4N_dth4*(th_dot^4) + 6*d3N_dth3*(th_dot^2)*th_ddot + ...
             3*d2N_dth2*(th_ddot^2) + 4*d2N_dth2*th_dot*th_jerk + dN_dth*th_snap;
             
    E_snap = d4E_dth4*(th_dot^4) + 6*d3E_dth3*(th_dot^2)*th_ddot + ...
             3*d2E_dth2*(th_ddot^2) + 4*d2E_dth2*th_dot*th_jerk + dE_dth*th_snap;
             
    D_snap = d4D_dth4*(th_dot^4) + 6*d3D_dth3*(th_dot^2)*th_ddot + ...
             3*d2D_dth2*(th_ddot^2) + 4*d2D_dth2*th_dot*th_jerk + dD_dth*th_snap;
    
    x_ref = [N; E; D];
    v_ref = [N_dot; E_dot; D_dot];
    a_ref = [N_ddot; E_ddot; D_ddot];
    j_ref = [N_jerk; E_jerk; D_jerk];
    s_ref = [N_snap; E_snap; D_snap];

    % =========================
    % Yaw reference using yaw_planner
    % =========================

    % Build state_ref struct for yaw_planner
    state_ref = struct();
    state_ref.pos = x_ref;
    state_ref.v = v_ref;
    state_ref.a = a_ref;

    % Build yaw_params struct for yaw_planner
    yaw_params = struct();
    yaw_params.last_psi = last_psi_persist;

    % Map legacy yaw mode names to yaw_planner modes
    switch lower(yaw_mode)
        case "align"
            % Legacy "align" mode maps to "tangent" (horizontal velocity alignment)
            actual_yaw_mode = "tangent";
        otherwise
            actual_yaw_mode = yaw_mode;
    end

    % Handle constant yaw: use psi_ref_const if provided
    if strcmpi(actual_yaw_mode, "constant")
        if isfield(params, 'psi_ref_const')
            yaw_params.psi_constant = params.psi_ref_const;
        else
            yaw_params.psi_constant = 0;
        end
    end

    % Handle constant_rate: use yaw_params if provided
    if strcmpi(actual_yaw_mode, "constant_rate")
        if isfield(params, 'yaw_params') && isfield(params.yaw_params, 'psi_rate')
            yaw_params.psi_rate = params.yaw_params.psi_rate;
        else
            yaw_params.psi_rate = pi/2;  % Default: 90 deg/s
        end
        if isfield(params, 'yaw_params') && isfield(params.yaw_params, 'psi_initial')
            yaw_params.psi_initial = params.yaw_params.psi_initial;
        elseif isfield(params, 'psi_ref_const')
            yaw_params.psi_initial = params.psi_ref_const;
        else
            yaw_params.psi_initial = 0;
        end
    end

    % Handle POI mode
    if strcmpi(actual_yaw_mode, "poi")
        if isfield(params, 'yaw_params') && isfield(params.yaw_params, 'poi')
            yaw_params.poi = params.yaw_params.poi;
        end
    end

    % Handle knife_edge mode
    if strcmpi(actual_yaw_mode, "knife_edge")
        if isfield(params, 'yaw_params') && isfield(params.yaw_params, 'knife_edge_dir')
            yaw_params.knife_edge_dir = params.yaw_params.knife_edge_dir;
        end
    end

    % Merge any additional yaw_params from input
    if isfield(params, 'yaw_params')
        fnames = fieldnames(params.yaw_params);
        for i = 1:length(fnames)
            if ~isfield(yaw_params, fnames{i})
                yaw_params.(fnames{i}) = params.yaw_params.(fnames{i});
            end
        end
    end

    % Call yaw_planner
    [psi_ref, psi_dot_ref, psi_ddot_ref] = yaw_planner(t, actual_yaw_mode, state_ref, yaw_params);

    % Update persistent for continuity
    last_psi_persist = psi_ref;
end

function [s, s1, s2, s3, s4] = s_curve_c3(tau)
% Standard C^3 "septic" polynomial:
    s  = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7;
    s1 = 140*tau^3 - 420*tau^4 + 420*tau^5 - 140*tau^6;
    s2 = 420*tau^2 - 1680*tau^3 + 2100*tau^4 - 840*tau^5;
    s3 = 840*tau   - 5040*tau^2 + 8400*tau^3 - 4200*tau^4;
    s4 = 840       - 10080*tau  + 25200*tau^2 - 16800*tau^3;
end
