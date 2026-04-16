function params = quad_params_indi()
%QUAD_PARAMS_INDI  Quadrotor parameters and control gains for INDI demo.
%
% NED frame: x=[N;E;D], D positive down.
% Body frame: bx forward, by right, bz down.

params.g = 9.81;
params.m = 1.0;                     % [kg]
params.J = diag([0.01, 0.01, 0.02]);% [kg m^2]

% Geometry / thrust coefficients (toy values, but consistent)
params.lx = 0.2;     % [m]
params.ly = 0.2;     % [m]
kT        = 1.0e-5;  % [N/(rad/s)^2]
kMz       = 2.0e-6;  % [Nm/(rad/s)^2]

params.kT  = kT;
params.kMz = kMz;

lx = params.lx;
ly = params.ly;

% Control effectiveness matrix G1: [mu_x; mu_y; mu_z; T] = G1 * omega.^2
params.G1 = [ ...
     ly*kT   -ly*kT   -ly*kT    ly*kT ;   % mu_x
     lx*kT    lx*kT   -lx*kT   -lx*kT ;   % mu_y
     -kMz      kMz     -kMz      kMz  ;    % mu_z
     kT       kT       kT       kT   ];   % T  (thrust magnitude)

% Hover rotor speed (approx)
params.omega_hover = sqrt((params.m*params.g) / (4*kT));

params.omega_min = 0;
params.omega_max = 2.5 * params.omega_hover;

params.psi_ref_const = 0.0;   % [rad] desired yaw (NED), can be nonzero later
params.yaw = "constant";

% % Outer-loop gains (Eq. 17)
% params.Kx = diag([4.0, 4.0, 6.0]);
% params.Kv = diag([4.0, 4.0, 4.5]);
% params.Ka = diag([0.8, 0.8, 0.8]);
%
% % Attitude gains
% params.Kxi    = diag([10.0, 10.0, 6.0]);    % angle error gains
% params.Komega = 0.8 * diag([4.0, 4.0, 3.0]);    % rate damping

% Outer-loop gains (Eq. 17)
% Gains are from Table-II
% [x,y,z]
params.Kx = diag([18.0, 18.0, 13.5]);
params.Kv = diag([7.8, 7.8, 5.9]);
params.Ka = diag([0.5, 0.5, 0.3]);

% Attitude gains
params.Kxi    = diag([175.0, 175.0, 82.0]);    % angle error gains
params.Komega = 0.8 * diag([19.0, 19.5, 19.2]);    % rate damping

% INDI filters (~30 Hz)
omega_c = 2*pi*30;
params.tau_f = 1 / omega_c;              % time constant

% Motor dynamics
params.tau_m = 0.02; %[s]

params.shape = "heart";

end