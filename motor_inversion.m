function omega = motor_inversion(mu_cmd, T_cmd, params)
%MOTOR_INVERSION  Solve G1 * omega.^2 = [mu_cmd; T_cmd].

G1 = params.G1;
omega_min = params.omega_min;
omega_max = params.omega_max;

u_des = [mu_cmd(:); T_cmd];
omega_sq = pinv(G1) * u_des;

omega_sq = max(omega_sq, 0);
omega = sqrt(omega_sq);

omega = min(max(omega, omega_min), omega_max);

end