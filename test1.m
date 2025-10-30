%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the motion of a supply vessel with a closed-loop 
% DP controller using equivalent added mass and damping matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
rng(1);

clear waveMotionRAO;

% Load vessel
load supply;

% Initial LF states
eta_0   = [0 0 0 0 0 0]';       % eta = [x y z phi theta psi]
nu_0    = [0 0 0 0 0 0]';       % nu = [u v w p q r]
x       = [eta_0; nu_0];

% Vessel mass 
m = vessel.main.m;
MRB = vessel.MRB; % Rigid body mass

% Hydrostatics
rho = 1025; g = 9.81;
Awp = vessel.main.Lwl * vessel.main.B * 0.8; % Waterplane displacement
GM_T = vessel.main.GM_T;
GM_L = vessel.main.GM_L;

% Linear restoring matrix
G = diag([MRB(1,1)*0.05^2, ...
          MRB(2,2)*0.05^2, ...
          rho*g*Awp, ...
          m*g*GM_T, ...
          m*g*GM_L, ...
          MRB(6,6)*0.05^2]);

%% Sea state and wave spectrum 
% Sea state
Hs      = 2.5;               % Significant wave height [m]
gamma   = 3.3; 
beta    = deg2rad(145);     % Wave direction relative to bow [rad]

Tz = 10 ;            % Zero-crossing period [s]
T0 = Tz / 0.710;    % Wave spectrum modal (peak) period [s] (Fossen 2021, Eq. 10.61)
w0 = 2*pi / T0;     % Wave spectrum modal (peak) frequency [rad/s]

spectrumParam = [Hs, w0, gamma];

maxFreq = 2*pi; % 3.0;                  % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 100;          % Number of wave frequency intervals (>50)
numDirections = 24;             % Number of wave directions (>15)

spreadingFlag = true;

% Reshape vessel data to use 0 to maxFreq
if vessel.forceRAO.w(end) > maxFreq
    w_index = find(vessel.forceRAO.w > maxFreq, 1) - 1;
    vessel.forceRAO.w = vessel.forceRAO.w(1:w_index); % frequency vector
    for DOF = 1:length(vessel.forceRAO.amp)
        vessel.forceRAO.amp{DOF} = vessel.forceRAO.amp{DOF}(1:w_index, :, :);
        vessel.forceRAO.phase{DOF} = vessel.forceRAO.phase{DOF}(1:w_index, :, :);
    end
end

omegaMax = vessel.forceRAO.w(end);

[S_M, omega, Amp, ~, ~, mu] = waveDirectionalSpectrum('JONSWAP', ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% Compute A_eq and B_eq using computeManeuveringModel
U_ship = 5;  % Ship forward speed (m/s)

% Compute encounter frequency (approximate shift for following/head seas)
omega_p = w0 - (w0^2 / g) * U_ship * cos(beta);

% Compute equivalent added mass and damping matrices
plotFlag = 0;  % Set to 1 to see plots of A(ω) and B(ω)
vessel = computeManeuveringModel(vessel, omega_p, plotFlag);

% Extract the diagonal elements of A_eq and B_eq for the first velocity
MA = diag([vessel.A_eq(1,1,1,1), ...
           vessel.A_eq(2,2,1,1), ...
           vessel.A_eq(3,3,1,1), ...
           vessel.A_eq(4,4,1,1), ...
           vessel.A_eq(5,5,1,1), ...
           vessel.A_eq(6,6,1,1)]);

D = diag([vessel.B_eq(1,1,1,1), ...
          vessel.B_eq(2,2,1,1), ...
          vessel.B_eq(3,3,1,1), ...
          vessel.B_eq(4,4,1,1), ...
          vessel.B_eq(5,5,1,1), ...
          vessel.B_eq(6,6,1,1)]);

% System inertia matrix
M = MRB + MA;

% Linear damping with additional tuning
zeta = [1 1 0.20 0.03 0.05 1];
D_tuned = diag([2*zeta(1)*sqrt(M(1,1)*G(1,1)), ...
                2*zeta(2)*sqrt(M(2,2)*G(2,2)), ...
                2*zeta(3)*sqrt(M(3,3)*G(3,3)), ...
                2*zeta(4)*sqrt(M(4,4)*G(4,4)), ...
                2*zeta(5)*sqrt(M(5,5)*G(5,5)), ...
                2*zeta(6)*sqrt(M(6,6)*G(6,6))]);

% Use the spectrum-weighted damping D, but can also blend with tuned damping if needed
% Option 1: Use only B_eq from computeManeuveringModel
% D_final = D;

% Option 2: Blend with zeta-tuned damping (recommended for better control performance)
alpha_blend = 0.5;  % Blending factor: 0 = all B_eq, 1 = all D_tuned
D_final = (1 - alpha_blend) * D + alpha_blend * D_tuned;

% Kinematic mapping
J6 = @(eta) eulerang(eta(4), eta(5), eta(6));

% Continuous-time state derivatives:
% eta_dot = J(eta)*nu
% nu_dot  = -M\G * eta - M\D * nu + M\tau
Minv = M \ eye(6);
ship_dynamics = @(x, tau) [J6(x(1:6)) * x(7:12); 
                           -Minv*(D_final*x(7:12) + G*x(1:6)) + Minv*tau]; 

% DP controller
eta_ref = zeros(6,1);
S = diag([1 1 0 0 0 1]);
alpha_p = 1.0; alpha_d = 1.0;         % gain scalars (tweakable)
Kp = alpha_p * G;
Kd = alpha_d * D_final;

use_integral = false; 
Ki = 0.02 * diag(diag(G));
eta_int = zeros(6,1);

tau_dp = @(eta, nu) S*(-Kp*(eta - eta_ref) - Kd*nu - (use_integral*Ki*eta_int));

%% Ship motion simulation
dt = 0.1;
T_final = 200;
T_initTransient = 20; 

t = (0:dt:T_final+T_initTransient-1);
N = numel(t);

% Low pass filtering psi
Tpsi = 12;
alpha = 1 - exp(-dt/Tpsi);
psi_lp = 0;

% Simulation state log
x_log           = zeros(N, 12);
tau_log         = zeros(N, 6);
eta_wf_log      = zeros(N, 6);
nu_wf_log       = zeros(N, 6);
nudot_wf_log    = zeros(N, 6);
zeta_log        = zeros(N, 1);    % wave elevation
y_eta_log       = zeros(N, 6);    % measured total position/orientation
y_nu_log        = zeros(N, 6);    % measured total velocities
y_nudot_log     = zeros(N, 6);    % measured total accelerations (WF only + LF approx)
psi_lp_log      = zeros(N,1);

y_eta = zeros(6,1);

%% Main Loop
for k = 1:N
    tk = t(k);
    eta = x(1:6); 
    nu = x(7:12); 

    % Low pass filtered heading
    psi_meas = ssa(y_eta(6));
    psi_err  = ssa(psi_meas - psi_lp);
    psi_lp   = ssa(psi_lp + alpha*psi_err);
    
    % 6-DOF wave-frequency (WF) motion
    [eta_wf, nu_wf, nudot_wf, zeta] = waveMotionRAO(tk, S_M, Amp, ...
        omega, mu, vessel, U_ship, 0, beta, numFreqIntervals);

    % Add DP control on LF states
    tau = tau_dp(eta, nu);

    if use_integral
        e_eta = eta_ref - eta;
        eta_int = eta_int + dt * (S * e_eta);
    end

    % Update states using rk4
    x = rk4(ship_dynamics, dt, x, tau);

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

    psi_lp_log(k,:) = psi_lp;
end

% Time-series
startIndex  = max(1, floor(T_initTransient/dt) + 1);
tt  = t(startIndex:end) - t(startIndex);
y_eta_log   = y_eta_log(startIndex:end,:);
y_nu_log    = y_nu_log(startIndex:end,:);
eta_wf_log  = eta_wf_log(startIndex:end,:);
nu_wf_log   = nu_wf_log(startIndex:end,:);
zeta_log    = zeta_log(startIndex:end);
tau_log     = tau_log(startIndex:end,:);
psi_lp_log  = psi_lp_log(startIndex:end,:);

%% Display equivalent matrices
fprintf('\n========================================\n');
fprintf('Equivalent Added Mass Matrix (MA diagonal):\n');
fprintf('A_eq(1,1) = %.2e kg (surge)\n', MA(1,1));
fprintf('A_eq(2,2) = %.2e kg (sway)\n', MA(2,2));
fprintf('A_eq(3,3) = %.2e kg (heave)\n', MA(3,3));
fprintf('A_eq(4,4) = %.2e kg·m² (roll)\n', MA(4,4));
fprintf('A_eq(5,5) = %.2e kg·m² (pitch)\n', MA(5,5));
fprintf('A_eq(6,6) = %.2e kg·m² (yaw)\n', MA(6,6));

fprintf('\nEquivalent Damping Matrix (D diagonal):\n');
fprintf('B_eq(1,1) = %.2e N·s/m (surge)\n', D(1,1));
fprintf('B_eq(2,2) = %.2e N·s/m (sway)\n', D(2,2));
fprintf('B_eq(3,3) = %.2e N·s/m (heave)\n', D(3,3));
fprintf('B_eq(4,4) = %.2e N·m·s/rad (roll)\n', D(4,4));
fprintf('B_eq(5,5) = %.2e N·m·s/rad (pitch)\n', D(5,5));
fprintf('B_eq(6,6) = %.2e N·m·s/rad (yaw)\n', D(6,6));
fprintf('========================================\n\n');

%% Plots
figure(1); clf;

% Spectrum plot
subplot(2,1,1); hold on; grid on;
if spreadingFlag
    plot(omega, S_M(:, floor(length(mu)/2)), 'LineWidth', 2);
    plot(omega, S_M(:, floor(length(mu)/4)), 'LineWidth', 2);
    plot(omega, S_M(:, length(mu)), 'LineWidth', 2);
    plot([w0 w0], [min(S_M,[],'all') max(S_M,[],'all')], 'LineWidth', 2);
    legend('\mu = 0°', '\mu = 45°', '\mu = 90°', sprintf('\\omega_0 = %.3f rad/s', w0));
else
    plot(omega, S_M(:,1), 'LineWidth', 2);
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
    subplot(7,1,i); grid on;
    plot(tt, Tscale(i)*y_eta_log(:,i), 'LineWidth', 1.6);
    ylabel(DOFtxt{i});
end
subplot(7,1,7)
plot(tt, Tscale(6)*psi_lp_log(:), 'LineWidth', 1.6)
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