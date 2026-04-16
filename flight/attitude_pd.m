function [alpha_cmd, xi_e] = attitude_pd(q_cmd, q, omega, omega_ref, omega_dot_ref, params)
%ATTITUDE_PD  Quaternion attitude + rate feedback.
%
% q_cmd, q : [qw;qx;qy;qz], body->inertial
% omega    : body rates [p;q;r]
%
% alpha_cmd: desired angular acceleration in body frame

Kxi    = params.Kxi;
Komega = params.Komega;

% Error quaternion: q_e = q_cmd * conj(q)
q_e = quat_mul(q_cmd, quat_conj(q));
q_e = quat_normalize(q_e);

if q_e(1) < 0
    q_e = -q_e; % Ensure the quaternion error is positive
end

qw = q_e(1);
qv = q_e(2:4);

angle = 2*acos(max(min(qw,1),-1));
if angle < 1e-6
    xi_e = [0;0;0];
else
    axis = qv / sin(angle/2);
    xi_e = angle * axis;
end

% xi_e_max = deg2rad(15);
% xi_e = saturate(xi_e,-xi_e_max,xi_e_max);

alpha_cmd = omega_dot_ref + Kxi*xi_e + Komega*(omega_ref - omega);

end