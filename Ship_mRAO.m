%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the motion of a supply vessel under closed-loop 
% DP-control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

rng(1);

% Load vessel 
load supply;

% Simulation flags
spreadingFlag = true; 
plotFlag = 0;

% Initial LF states
eta_0   = [0 0 0 0 0 0]';       % eta = [x y z phi theta psi]
nu_0    = [0 0 0 0 0 0]';       % nu = [u v w p q r]
x       = [eta_0; nu_0];

% vessel speed and heading
psi_ship = 0;
U_ship = 0;

% Simulation parametes
h = 0.1; % sampeling time
T_final = 200; 
T_initTransient = 50;

t = (0:h:T_final+T_initTransient-1);
N = numel(t);

% Sea state and wave spectrum
Hs = 2.5; % Significant wave height [m]
omega_p = 0.2;
gamma = 3.3; % Peakedness factor
beta_wave = deg2rad(145); % Wave direction relative to bow [rad]

spectrumType = 'JONSWAP';
spectrumParam = [Hs, omega_p, gamma];

maxFreq = 3.0;              % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 100;     % Number of wave frequency intervals (>50)
numDirections = 24;         % Number of wave directions (>15)

% Reshape vessel data o use 0 to maxFreq
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

%% Compute A_eq and B_eq
g = 9.81;
omega_p = omega_p - (omega_p^2 / g) * U_ship * cos(beta_wave); % Shifted encounter frequency
vessel = computeManeuveringModel(vessel, omega_p, plotFlag);

%% Vessel manuvering model

MRB = vessel.MRB;

% Extract the diagonal elements of A_eq and B_eq for the first velocity
MA = diag([vessel.A_eq(1,1,1,1), ...
           vessel.A_eq(2,2,1,1), ...
           vessel.A_eq(3,3,1,1), ...
           vessel.A_eq(4,4,1,1), ...
           vessel.A_eq(5,5,1,1), ...
           vessel.A_eq(6,6,1,1)]);

M = MRB + MA;

D = diag([vessel.B_eq(1,1,1,1), ...
          vessel.B_eq(2,2,1,1), ...
          vessel.B_eq(3,3,1,1), ...
          vessel.B_eq(4,4,1,1), ...
          vessel.B_eq(5,5,1,1), ...
          vessel.B_eq(6,6,1,1)]);


% Linear restroing matrix
m = vessel.main.m; % Vessel mass
rho = 1025; g = 9.81; 
Awp = vessel.main.Lwl * vessel.main.B * 0.8; % Waterplane displacement
GM_T = vessel.main.GM_T;
GM_L = vessel.main.GM_L;

R33 = rho*g*Awp; 
R44 = m*g*GM_T; 
R55 = m * g * GM_L; % Calculate the restoring force for the roll motion

G = diag([0, 0, R33, R44, R55, 0]);

% Kinematic mapping
J = @(eta) eulerang(eta(4), eta(5), eta(6));

% Continuous-time state derivatives:
% eta_dot = J(eta)*nu
% nu_dot  = -M\G * eta - M\D * nu + M\tau
Minv = M \ eye(6);
ship_dynamics = @(x, tau) [J(x(1:6)) * x(7:12); 
                           -Minv*(D*x(7:12) + G*x(1:6)) + Minv*tau]; 

% DP controller
eta_ref = zeros(6,1);
S = diag([1 1 0 0 0 1]);
alpha_p = 1.0; alpha_d = 1.0;         % gain scalars (tweakable)
Kp = alpha_p * G;
Kd = alpha_d * D;

use_integral = false; 
Ki = 0.02 * diag(diag(G));
eta_int = zeros(6,1);

tau_dp = @(eta, nu) S*(-Kp*(eta - eta_ref) - Kd*nu - (use_integral*Ki*eta_int));

% Low pass filterting psi
Tpsi = 12;
alpha = 1 - exp(-h/Tpsi);
psi_lp = 0;

% Simualtion state log
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
        omega, mu, vessel, U_ship, psi_ship, beta_wave, numFreqIntervals);

    % Add DP control on LF states
    tau = tau_dp(eta, nu);

    if use_integral
        e_eta = eta_ref - eta;
        eta_int = eta_int + dt * (S * e_eta);
    end

    % Update states using rk4
    x = rk4(ship_dynamics, h, x, tau);

    % "Measured" total output = LF + WF (per Fossen 2021)
    y_eta   = eta + eta_wf;
    y_nu    = nu  + nu_wf;
    y_nudot = nudot_wf; % LF accel not kept explicitly here

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
startIndex  = max(1, floor(T_initTransient/h) + 1);
tt  = t(startIndex:end) - t(startIndex);
y_eta_log   = y_eta_log(startIndex:end,:);
y_nu_log    = y_nu_log(startIndex:end,:);
eta_wf_log  = eta_wf_log(startIndex:end,:);
nu_wf_log   = nu_wf_log(startIndex:end,:);
zeta_log    = zeta_log(startIndex:end);
tau_log     = tau_log(startIndex:end,:);
psi_lp_log  = psi_lp_log(startIndex:end,:);

%% Plots
figure(1); clf;

% Spectrum plot
subplot(2,1,1); hold on; grid on;
if spreadingFlag
    plot(omega, S_M(:, floor(length(mu)/2)), 'LineWidth', 2);
    plot(omega, S_M(:, floor(length(mu)/4)), 'LineWidth', 2);
    plot(omega, S_M(:, length(mu)), 'LineWidth', 2);
    plot([omega_p omega_p], [min(S_M,[],'all') max(S_M,[],'all')], 'LineWidth', 2);
    legend('\mu = 0°', '\mu = 45°', '\mu = 90°', sprintf('\\omega_0 = %.3f rad/s', omega_p));
else
    plot(omega, S_M(:,1), 'LineWidth', 2);
    plot([omega_p omega_p], [min(S_M,[],'all') max(S_M,[],'all')], 'LineWidth', 2);
    legend('S(\Omega)', sprintf('\\omega_0 = %.3f rad/s', omeaga_p));
end
xlabel('\Omega (rad/s)'); ylabel('m^2 s');
title('Directional wave spectrum');

% Wave elevation
subplot(2,1,2); hold on; grid on;
plot(tt, zeta_log, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('m');
title(sprintf('Wave Elevation: \\beta_{wave} = %.1f° , H_s = %.1f m', rad2deg(beta_wave), Hs));

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
