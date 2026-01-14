function [Omega_ref, omega_dot_ref] = flatness(R, tau, j_ref, s_ref, psi_dot_ref, psi_ddot_ref)
%FLATNESS_TAL_KARAMAN
%   Tal & Karaman 2021 makalesi Denklem (14) ve (15) implementasyonudur.
%
%   Girdiler:
%   R      : Gövde -> Inertial rotasyon matrisi (Mevcut durum veya Command)
%   tau    : Spesifik itki büyüklüğü (m/s^2) -> norm(a_cmd + g*e3)
%   j_ref  : Jerk referansı (Inertial) [3x1]
%   s_ref  : Snap referansı (Inertial) [3x1]
%
%   Çıktılar:
%   Omega_ref : Açısal hız referansı (Body frame) [3x1]
%   alpha_ref : Açısal ivme referansı (Body frame) [3x1] (Feedforward)
%   debug     : Debug struct with singularity diagnostics

% Sabitler
e3 = [0;0;1];
iz = e3;

% Gövde eksenleri
bx = R(:,1);
by = R(:,2);
bz = R(:,3); % Makalede itki yönü bz olarak tanımlı (Fig 1)

% Mevcut açısal hız (Omega) tahmini (Denklem 14 ve 15 türetimi için)
% Not: Makaledeki formülasyonda, referans türetilirken sistemin o anki
% R ve tau değerleri kullanılır.

% --- S Matrisinin Hesabı (Denklem 12) ---
% S, yaw rate ile body rates arasındaki dönüşümdür.
% r_psi vektörü projeksiyonu:
r_psi_sq = bx(1)^2 + bx(2)^2;

% [SOFT-EPSILON]: Avoid singularity at 90 deg pitch
% Reverted to small epsilon, handling singularity via DLS solver instead
eps_s = 1e-3; % Slightly increased from 1e-6 to 1e-3
r_psi_sq_safe = max(r_psi_sq, eps_s);

% S = [-bx(2), bx(1)] * [0 -bz(1) by(1); 0 -bz(2) by(2)] / r_sq
term1 = [-bx(2), bx(1)] / r_psi_sq_safe;
mat2  = [0, -bz(1), by(1);
    0, -bz(2), by(2)];
S = term1 * mat2;

% --- Denklem (14): Omega_ref Hesabı ---
% M matrisi (4x4): [tau * R * [iz]x , -bz; S, 0]

iz_skew = [0 -1 0; 1 0 0; 0 0 0]; % [iz]x

% M matrisi (4x4)
top_left  = tau * R * iz_skew; % 3x3
top_right = -bz;               % 3x1 (-bz as in jerk = tau R [iz]x Omega - tau_dot bz)
bot_left  = S;                 % 1x3
bot_right = 0;                 % 1x1

M = [top_left, top_right;
    bot_left, bot_right];

% Sağ taraf vektörü
rhs_14 = [j_ref; psi_dot_ref]; % 4x1

% Çözüm: [Omega_ref; tau_dot_ref]
% Use Damped Least Squares (DLS) for robustness near singularity
% Instead of inv(M), use (M'*M + lambda*I)^-1 * M'
lambda_dls = 1e-5; % Reduced damping factor for better tracking
M_sq = M' * M;
I_dls = eye(size(M_sq));
sol_14 = (M_sq + lambda_dls * I_dls) \ (M' * rhs_14);

Omega_ref   = sol_14(1:3);

% Safety saturation for Omega_ref
max_omega = 40.0; % rad/s (approx 2300 deg/s) - increased limit
if norm(Omega_ref) > max_omega
    Omega_ref = (Omega_ref / norm(Omega_ref)) * max_omega;
end

tau_dot_ref = sol_14(4);

% --- Denklem (15): Alpha_ref (Omega_dot) Hesabı ---
% Sağ taraf vektöründeki Coriolis terimleri

Omega_skew = [0 -Omega_ref(3) Omega_ref(2);
    Omega_ref(3) 0 -Omega_ref(1);
    -Omega_ref(2) Omega_ref(1) 0];

% Terim: R * (2*tau_dot + tau*[Omega]x) * [iz]x * Omega
% Note: Using iz_skew instead of transpose for correct physical derivation
term_coriolis_1 = R * (2*tau_dot_ref*eye(3) + tau*Omega_skew) * iz_skew * Omega_ref;

term_coriolis_2 = 0; % S_dot terms omitted

rhs_15 = [s_ref; psi_ddot_ref] - [term_coriolis_1; term_coriolis_2];

% M matrisi aynıdır
sol_15 = (M_sq + lambda_dls * I_dls) \ (M' * rhs_15);

omega_dot_ref = sol_15(1:3); % Feedforward acceleration

% --- Debug Output ---
% if nargout >= 3
%     debug.r_psi_sq = r_psi_sq;
%     debug.pitch = asin(-R(3,1));  % Pitch angle from R
%     debug.cond_M = cond(M);       % Condition number of M matrix
%     debug.Omega_norm = norm(Omega_ref);
%     debug.Omega_ref = Omega_ref;
%     debug.tau = tau;
% end

end
