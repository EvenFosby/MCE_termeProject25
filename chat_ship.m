%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSV wave-induced motion + simplified closed-loop DP using MSS toolbox
% - LF ship model: 6-DOF linear mass + hydrostatics + linear damping
% - WF motion: waveMotionRAO (uses vessel.motionRAO)
% - Total measured outputs: y = LF + WF (per Fossen 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
rng(1);

%--- Reset any persistent RAO tables inside waveMotionRAO -----------------
clear waveMotionRAO;

%--- Load vessel (must contain vessel.motionRAO.*) ------------------------
load supply   % provides 'vessel'
assert(isfield(vessel,'motionRAO'), 'supply.mat must include vessel.motionRAO');

%--- Initial LF states ----------------------------------------------------
eta  = zeros(6,1);   % [x y z phi theta psi]
nu   = zeros(6,1);   % [u v w p q r]
x    = [eta; nu];

%--- Mass/Inertia ---------------------------------------------------------
m   = vessel.main.m;
MRB = vessel.MRB;              % rigid-body
MA  = vessel.A(:,:,1);         % added mass (at low freq)
M   = MRB + MA;                % total inertia

%--- Hydrostatics ---------------------------------------------------------
rho = 1025; g = 9.81;
Awp   = vessel.main.Lwl * vessel.main.B * 0.8;  % waterplane approx
GM_T  = vessel.main.GM_T;
GM_L  = vessel.main.GM_L;

% Linear restoring (G) ~ diag(k_surge, k_sway, k_heave, k_roll, k_pitch, k_yaw)
G = diag([ M(1,1)*0.05^2, ...
           M(2,2)*0.05^2, ...
           rho*g*Awp, ...
           m*g*GM_T, ...
           m*g*GM_L, ...
           M(6,6)*0.05^2 ]);

% Linear damping (modal ζ)
zeta = [1 1 0.20 0.03 0.05 1];
D = diag([ 2*zeta(1)*sqrt(M(1,1)*G(1,1)), ...
           2*zeta(2)*sqrt(M(2,2)*G(2,2)), ...
           2*zeta(3)*sqrt(M(3,3)*G(3,3)), ...
           2*zeta(4)*sqrt(M(4,4)*G(4,4)), ...
           2*zeta(5)*sqrt(M(5,5)*G(5,5)), ...
           2*zeta(6)*sqrt(M(6,6)*G(6,6)) ]);

%--- Convenience handles ---------------------------------------------------
J6 = @(eta_) eulerang(eta_(4), eta_(5), eta_(6)); % MSS: 6x6 kinematic map

% Continuous-time state derivatives:
% eta_dot = J(eta)*nu
% nu_dot  = -M\G * eta - M\D * nu + M\tau
Minv = M \ eye(6);
ship_dynamics = @(eta_,nu_,tau_) [ J6(eta_) * nu_ ;
                                  -Minv*(G*eta_ + D*nu_) + Minv*tau_ ];

%--- DP Controller (operate on LF states ONLY) ----------------------------
eta_ref = zeros(6,1);                 % station-keeping at origin
S = diag([1 1 0 0 0 1]);              % control {x,y,psi}; others open
alpha_p = 1.0; alpha_d = 1.0;         % gain scalars (tweakable)
Kp = alpha_p * G;
Kd = alpha_d * D;

use_integral = false;                 % set true to enable integral action
Ki = 0.02 * diag(diag(G));            % small integral if used
eta_int = zeros(6,1);

dp_tau = @(eta_,nu_) ...
    S * ( -Kp*(eta_ - eta_ref) - Kd*nu_ - (use_integral * Ki * eta_int) );

%--- Sea state & directional spectrum -------------------------------------
Hs    = 1.5;                 % significant wave height [m]
Tz    = 10;                 % zero-crossing period [s]
T0    = Tz / 0.710;         % modal (peak) period [s] (Fossen 2021, Eq. 10.61)
w0    = 2*pi / T0;          % modal (peak) frequency [rad/s]
gamma = 3.3;                % JONSWAP peak enhancement
beta  = deg2rad(145);       % wave direction rel. bow [rad]

maxFreq          = 3.0;     % RAO computation cap [rad/s]
numFreqIntervals = 60;      % > 50 recommended
numDirections    = 24;      % > 15 recommended
spreadingFlag    = false;   % no directional spreading (set true to enable)

% Trim RAO tables to [0, maxFreq]  --- *** motionRAO *** ------------------
if vessel.motionRAO.w(end) > maxFreq
    w_index = find(vessel.motionRAO.w > maxFreq, 1) - 1;
    vessel.motionRAO.w = vessel.motionRAO.w(1:w_index); % freq vector
    for DOF = 1:length(vessel.motionRAO.amp)
        vessel.motionRAO.amp{DOF}   = vessel.motionRAO.amp{DOF}(1:w_index, :, :);
        vessel.motionRAO.phase{DOF} = vessel.motionRAO.phase{DOF}(1:w_index, :, :);
    end
end
omegaMax = vessel.motionRAO.w(end);

% Build directional spectrum + complex RAO amplitude array
spectrumParam = [Hs, w0, gamma];
[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum( ...
    'JONSWAP', spectrumParam, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%--- Simulation setup ------------------------------------------------------
dt              = 0.1;
T_final         = 200;
T_initTransient = 20;

t = 0:dt:(T_final + T_initTransient);
N = numel(t);

% Logs
x_log           = zeros(N, 12);
tau_log         = zeros(N, 6);
eta_wf_log      = zeros(N, 6);
nu_wf_log       = zeros(N, 6);
nudot_wf_log    = zeros(N, 6);
zeta_log        = zeros(N, 1);    % wave elevation
y_eta_log       = zeros(N, 6);    % measured total position/orientation
y_nu_log        = zeros(N, 6);    % measured total velocities
y_nudot_log     = zeros(N, 6);    % measured total accelerations (WF only + LF approx)

% Ship speed and heading for wave kinematics in RAO evaluation
% For DP stationkeeping, set U ≈ 0 and psi = eta(6)
U_ship  = 0;  % m/s

%--- Main loop ------------------------------------------------------------
for k = 1:N
    tk  = t(k);
    eta = x(1:6);
    nu  = x(7:12);
    psi = eta(6);

    % WF motion from RAO at time tk
    [eta_wf, nu_wf, nudot_wf, zeta] = waveMotionRAO( ...
        tk, S_M, Amp, Omega, mu, vessel, U_ship, psi, beta, numFreqIntervals);

    % Ensure column vectors
    eta_wf   = eta_wf(:);
    nu_wf    = nu_wf(:);
    nudot_wf = nudot_wf(:);

    % DP control on LF states only (recommended; avoids wave chasing)
    tau = dp_tau(eta, nu);

    % Optional integral (LF only on controlled DOFs)
    if use_integral
        e_eta = eta_ref - eta;
        eta_int = eta_int + dt * (S * e_eta);
    end

    % Integrate LF dynamics (Euler explicit; small dt is fine here)
    xdot = ship_dynamics(eta, nu, tau);
    x    = x + dt * xdot;

    % "Measured" total output = LF + WF (per Fossen 2021)
    y_eta   = eta + eta_wf;
    y_nu    = nu  + nu_wf;
    y_nudot =                nudot_wf; % LF accel not kept explicitly here

    % Log
    x_log(k,:)        = [eta.' nu.'];
    tau_log(k,:)      = tau.';
    eta_wf_log(k,:)   = eta_wf.';
    nu_wf_log(k,:)    = nu_wf.';
    nudot_wf_log(k,:) = nudot_wf.';
    zeta_log(k)       = zeta;
    y_eta_log(k,:)    = y_eta.';
    y_nu_log(k,:)     = y_nu.';
    y_nudot_log(k,:)  = y_nudot.';
end

%--- Trim initial transient ----------------------------------------------
k0  = max(1, floor(T_initTransient/dt) + 1);
tt  = t(k0:end) - t(k0);
y_eta_log   = y_eta_log(k0:end,:);
y_nu_log    = y_nu_log(k0:end,:);
eta_wf_log  = eta_wf_log(k0:end,:);
nu_wf_log   = nu_wf_log(k0:end,:);
zeta_log    = zeta_log(k0:end);
tau_log     = tau_log(k0:end,:);

%--- Plots ----------------------------------------------------------------
figure(1); clf;

% Spectrum plot
subplot(2,1,1); hold on; grid on;
if spreadingFlag
    plot(Omega, S_M(:, floor(length(mu)/2)), 'LineWidth', 2);
    plot(Omega, S_M(:, floor(length(mu)/4)), 'LineWidth', 2);
    plot(Omega, S_M(:, length(mu)), 'LineWidth', 2);
    plot([w0 w0], [min(S_M,[],'all') max(S_M,[],'all')], 'LineWidth', 2);
    legend('\mu = 0°', '\mu = 45°', '\mu = 90°', sprintf('\\omega_0 = %.3f rad/s', w0));
else
    plot(Omega, S_M(:,1), 'LineWidth', 2);
    plot([w0 w0], [min(S_M,[],'all') max(S_M,[],'all')], 'LineWidth', 2);
    legend('S(\Omega)', sprintf('\\omega_0 = %.3f rad/s', w0));
end
xlabel('\Omega (rad/s)'); ylabel('m^2 s');
title('Directional wave spectrum');

% Wave elevation
subplot(2,1,2); hold on; grid on;
plot(tt, zeta_log, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('m');
title(sprintf('Wave Elevation: \\beta_{wave} = %.1f° , H_s = %.1f m', rad2deg(beta), Hs));

% 6-DOF total WF+LF positions (η)
figure(2); clf;
DOFtxt = {'x (m)','y (m)','z (m)','\phi (deg)','\theta (deg)','\psi (deg)'};
Tscale = [1 1 1 180/pi 180/pi 180/pi];
for i = 1:6
    subplot(6,1,i); grid on;
    plot(tt, Tscale(i)*y_eta_log(:,i), 'LineWidth', 1.6);
    ylabel(DOFtxt{i});
end
xlabel('Time (s)');
sgtitle('Measured positions & Euler angles (LF + WF)');

% 6-DOF measured velocities
figure(3); clf;
DOFtxt_v = {'u (m/s)','v (m/s)','w (m/s)','p (deg/s)','q (deg/s)','r (deg/s)'};
for i = 1:6
    subplot(6,1,i); grid on;
    plot(tt, Tscale(i)*y_nu_log(:,i), 'LineWidth', 1.6);
    ylabel(DOFtxt_v{i});
end
xlabel('Time (s)');
sgtitle('Measured velocities (LF + WF)');

% Control effort
figure(4); clf; grid on;
plot(tt, tau_log, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\tau (N / N·m)');
title('DP control effort (LF only)');
legend('X','Y','Z','K','M','N');
