function q_inc = compute_incremental_attitude_cmd(R_curr, tau_bz_c, psi_ref)
% Same attitude increment logic used in mission/INDI demos.

t_norm = norm(tau_bz_c);
if t_norm < 1e-6
    t_des_in = [0; 0; -1];
else
    t_des_in = tau_bz_c / t_norm;
end

bz_des_in = -t_des_in;
bz_des_body = R_curr' * bz_des_in;

current_z = [0; 0; 1];
cross_prod = cross(current_z, bz_des_body);
dot_prod = dot(current_z, bz_des_body);

if dot_prod < -0.9999
    q_tilt = [0; 1; 0; 0];
else
    s = sqrt(2 * (1 + dot_prod));
    q_tilt = [0.5 * s; cross_prod / s];
end
q_tilt = quat_normalize(q_tilt);

R_tilt = quat_to_R(q_tilt);
R_new_in = R_curr * R_tilt;
psi_curr = atan2(R_new_in(2, 1), R_new_in(1, 1));

psi_err = psi_ref - psi_curr;
while psi_err > pi,  psi_err = psi_err - 2*pi; end
while psi_err < -pi, psi_err = psi_err + 2*pi; end

q_yaw = [cos(psi_err/2); 0; 0; sin(psi_err/2)];
q_inc = quat_mul(q_tilt, q_yaw);
if q_inc(1) < 0
    q_inc = -q_inc;
end
end