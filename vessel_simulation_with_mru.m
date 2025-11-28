%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Vessel Motion Simulation with MRU Noise
% This script simulates vessel motion under wave disturbances and applies
% realistic MRU sensor noise to generate both true and measured motion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% PART 1: VESSEL MOTION SIMULATION (from Ship_motion_RAO.m)
fprintf('=== PART 1: VESSEL MOTION SIMULATION ===\n');

rng(1);

% Load vessel
load supply;

% Simulation Flags
spreadingFlag = true;
plotFlag = false;
forceRaoFlag = true;
useIntegralAction = false;
recordFlag = true;

% Simulation parameters
h = 0.1;
T_final = 300;
T_initTransient = 60;

t = 0:h:T_final+T_initTransient-1;
N = numel(t);

% Control objective
psi = 0;
U = 0;

%% Sea state and wave spectrum
Hs = 2.5;
Tz = 6;

T0 = Tz / 0.710;
w0 = 2*pi / T0;

gamma = 3.3;
beta = deg2rad(140);

spectrumType = 'JONSWAP';
spectrumParam = [Hs, w0, gamma];

maxFreq = 3.0;
numFreqIntervals = 60;
numDirections = 48;

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

[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum(spectrumType, ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% Vessel model
MRB = vessel.MRB;

% Compute A_eq and B_eq
g = 9.81;
omega_p = w0 - (w0^2 / g) * U * cos(beta); % Shifted encounter frequency
vessel = computeManeuveringModel(vessel, omega_p, plotFlag);

% Extract the diagonal elements of A_eq and B_eq for the first velocity
MA = diag([vessel.A_eq(1,1,1,1), ...
           vessel.A_eq(2,2,1,1), ...
           vessel.A_eq(3,3,1,1), ...
           vessel.A_eq(4,4,1,1), ...
           vessel.A_eq(5,5,1,1), ...
           vessel.A_eq(6,6,1,1)]);

M = MRB + MA;
Minv = M \ eye(6);

D = diag([vessel.B_eq(1,1,1,1), ...
          vessel.B_eq(2,2,1,1), ...
          vessel.B_eq(3,3,1,1), ...
          vessel.B_eq(4,4,1,1), ...
          vessel.B_eq(5,5,1,1), ...
          vessel.B_eq(6,6,1,1)]);

% Linear restoring matrix
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

ship_dynamics = @(x, tau) [J(x(1:6)) * x(7:12);
                           -Minv*(D*x(7:12) + G*x(1:6)) + Minv*tau];

%% DP controller
S = [1 1 0 0 0 1]';

% Computing PID-gains using Algorithm 15.2 from (Fossen, 2021)
omega_b1 = 0.05; omega_b2 = 0.05; omega_b6 = 0.03;
omega_b = [omega_b1, omega_b2, 0, 0, 0, omega_b6];
Omega_b = diag(omega_b);

zeta_pid1 = 0.80; zeta_pid2 = 0.80; zeta_pid6 = 1;
zeta_pid = [zeta_pid1, zeta_pid2, 0, 0, 0, zeta_pid6];
Zeta_pid = diag(zeta_pid);

omega_n = zeros(1,6);
for i = 1:length(omega_n)
    omega_n(i) = omega_b(i) / ( sqrt(1 - 2*zeta_pid(i)^2 + sqrt(4*zeta_pid(i)^4 - 4*zeta_pid(i)^2 + 2) ) );
end
Omega_n = diag(omega_n);

% Assuming roll, pitch and yaw is small => J_Theta(eta) = I
Kp = M*Omega_n^2;

Kd = 2.*M*Zeta_pid*Omega_n;

Ki = 0.10*Kp*Omega_n;

% PID controller
tau_pid = @(eta, nu, eta_int) S.*(eulerang(eta(4), eta(5), eta(6))'*(-Kp*eta ...
    - Kd*eulerang(eta(4), eta(5), eta(6))*nu - (useIntegralAction*Ki*eta_int)));

% Heading lowpass filter
T_psi = 12;
alpha = h/(T_psi + h);
psi_lp = 0;

%% Main simulation loop
fprintf('Running vessel motion simulation...\n');

% Initial vessel states
eta_0   = [0 0 0 0 0 0]';       % eta = [x y z phi theta psi]
nu_0    = [0 0 0 0 0 0]';       % nu = [u v w p q r]
x       = [eta_0; nu_0];

eta_int = [0 0 0 0 0 0]';

% desired vessel states
eta_d   = [0 0 0 0 0 0]';
nu_d    = [0 0 0 0 0 0]';
x_d     = [eta_0; nu_0];

% Preallocate data log
x_log           = zeros(N, 12);
tau_log         = zeros(N, 6);
tau_ctrl_log    = zeros(N, 6);
eta_wf_log      = zeros(N,6);
nu_wf_log       = zeros(N, 6);
nudot_wf_log    = zeros(N, 6);
zeta_log        = zeros(N, 1);    % wave elevation
y_eta_log       = zeros(N, 6);    % measured total position/orientation
y_nu_log        = zeros(N, 6);    % measured total velocities
y_nudot_log     = zeros(N, 6);    % measured total accelerations (WF only + LF approx)
psi_lp_log      = zeros(N,1);

simdata_fRAO = zeros(N, 7);
simdata_mRAO = zeros(N, 19);

% Recording arrays (only used when recordFlag = true)
eta_recorded = zeros(N, 6);
nu_recorded = zeros(N, 6);
t_recorded = zeros(N, 1);
record_counter = 0;

for k = 1:N
    tk = t(k);
    eta = x(1:6);
    nu = x(7:12);

    % Lowpass filtering heading
    psi_err = eta(6) - psi_lp;
    psi_lp = psi_lp + alpha*psi_err;
    psi_lp = ssa(psi_lp);

    % Add DP control on LF states
    tau_control = tau_pid(eta, nu, eta_int);

    if useIntegralAction
        e_eta = eta_d - eta;
        eta_int = eta_int + h * (S .* e_eta);
    end

    if forceRaoFlag
        % Force RAO
        [tau_wave, zeta_wave] = waveForceRAO_v1(tk, S_M, Amp, Omega, mu, ...
        vessel, U, psi, beta, numFreqIntervals);
        simdata_fRAO(k, :) = [tau_wave', zeta_wave']; % Log force data

        eta_wf = zeros(6,1);
        nu_wf = zeros(6,1);
        nudot_wf = zeros(6,1);
    else
        % Motion RAO
        [eta_wf, nu_wf, nudot_wf, zeta_wave] = waveMotionRAO_v1(tk, ...
            S_M, Amp, Omega, mu, vessel, U, psi, beta, numFreqIntervals);
        simdata_mRAO(k, :) = [eta_wf', nu_wf', nudot_wf', zeta_wave]; % Log motion data

        tau_wave = zeros(6,1);
    end

    tau = tau_control + tau_wave;

    % Update states using rk4
    x = rk4(ship_dynamics, h, x, tau);

    % "Measured" total output = LF + WF (per Fossen 2021)
    y_eta   = eta + eta_wf;
    y_nu    = nu  + nu_wf;
    y_nudot = nudot_wf; % LF accel not kept explicitly here

    % Log everything
    x_log(k,:)        = [eta.' nu.'];
    tau_log(k,:)      = tau.';
    tau_ctrl_log(k,:) = tau_control.';
    eta_wf_log(k,:)   = eta_wf.';
    nu_wf_log(k,:)    = nu_wf.';
    nudot_wf_log(k,:) = nudot_wf.';
    zeta_log(k)       = zeta_wave;
    y_eta_log(k,:)    = y_eta.';
    y_nu_log(k,:)     = y_nu.';
    y_nudot_log(k,:)  = y_nudot.';
    psi_lp_log(k,:)   = psi_lp;

    % Record eta and nu if recordFlag is true
    if recordFlag
        if forceRaoFlag
            record_counter = record_counter + 1;
            eta_recorded(record_counter, :) = eta.';
            nu_recorded(record_counter, :) = nu.';
            t_recorded(record_counter) = tk;
        else
            record_counter = record_counter + 1;
            eta_recorded(record_counter, :) = y_eta.';
            nu_recorded(record_counter, :) = y_nu.';
            t_recorded(record_counter) = tk;
        end

    end

end

% Trim recorded data to actual size and save if recording was enabled
if recordFlag
    eta_recorded = eta_recorded(1:record_counter, :);
    nu_recorded = nu_recorded(1:record_counter, :);
    t_recorded = t_recorded(1:record_counter);

    % Save recorded data to file
    motion_data.time = t_recorded;
    motion_data.eta = eta_recorded;
    motion_data.nu = nu_recorded;
    motion_data.description = 'Vessel eta (position/orientation) and nu (velocities) time series';
    motion_data.eta_labels = {'x (m)', 'y (m)', 'z (m)', 'phi (rad)', 'theta (rad)', 'psi (rad)'};
    motion_data.nu_labels = {'u (m/s)', 'v (m/s)', 'w (m/s)', 'p (rad/s)', 'q (rad/s)', 'r (rad/s)'};

    save('vessel_motion_data.mat', 'motion_data');
    fprintf('Recorded %d samples of eta and nu to vessel_motion_data.mat\n', record_counter);
end

% After your main loop, ignore transient
idx0 = max(1, floor(T_initTransient/h) + 1);
zeta = simdata_mRAO(idx0:end,19);

% Spectrum validation
fprintf('std(zeta)=%.3f m  (expected ~ %.3f m)\n', std(zeta), Hs/4);

dOmega = Omega(2)-Omega(1);                % rad/s
if spreadingFlag
    dmu = 2*pi/numDirections;              % radians
    m0  = sum(S_M(:))*dOmega*dmu;          % m^2
else
    m0  = sum(S_M(:,1))*dOmega;            % m^2
end
fprintf('m0 from S_M = %.3f m^2,  std(zeta)^2 = %.3f m^2\n', m0, var(zeta));

%% PART 2: MRU NOISE GENERATION (from mru_noise_sim.m)
fprintf('\n=== PART 2: MRU NOISE GENERATION ===\n');

% Extract data from recorded motion
eta_true = eta_recorded;
nu_true = nu_recorded;
t_mru = t_recorded;
N_mru = length(t_mru);
h_mru = t_mru(2) - t_mru(1);

fprintf('  Processing %d samples, dt = %.4f s, duration = %.1f s\n', N_mru, h_mru, t_mru(end));

%% MRU Noise Specifications
% Position noise standard deviations (RMS)
sigma_eta = [
    0.10;      % Surge (x) - 10 cm
    0.10;      % Sway (y) - 10 cm
    0.05;      % Heave (z) - 5 cm
    deg2rad(0.03);  % Roll (phi) - 0.03 degrees RMS
    deg2rad(0.03);  % Pitch (theta) - 0.03 degrees RMS
    deg2rad(0.10);  % Yaw (psi) - 0.10 degrees RMS
];

% Velocity noise standard deviations (RMS)
sigma_nu = [
    0.05;      % Surge velocity (u) - 5 cm/s
    0.05;      % Sway velocity (v) - 5 cm/s
    0.03;      % Heave velocity (w) - 3 cm/s
    deg2rad(0.05);  % Roll rate (p) - 0.05 deg/s
    deg2rad(0.05);  % Pitch rate (q) - 0.05 deg/s
    deg2rad(0.08);  % Yaw rate (r) - 0.08 deg/s
];

% Optional: Small sensor bias (can drift slowly over time)
bias_eta = zeros(6, 1);
bias_nu = zeros(6, 1);

% Optional: Bias drift rates
drift_rate_eta = sigma_eta * 0.001;
drift_rate_nu = sigma_nu * 0.001;

fprintf('\nMRU Noise Model (Kongsberg Maritime state-of-the-art):\n');
fprintf('Position noise (RMS):\n');
fprintf('  Surge/Sway: %.1f cm, Heave: %.1f cm\n', sigma_eta(1)*100, sigma_eta(3)*100);
fprintf('  Roll/Pitch: %.3f deg, Yaw: %.3f deg\n', rad2deg(sigma_eta(4)), rad2deg(sigma_eta(6)));
fprintf('Velocity noise (RMS):\n');
fprintf('  Linear velocities: %.1f cm/s\n', sigma_nu(1)*100);
fprintf('  Angular velocities: %.3f deg/s\n', rad2deg(sigma_nu(4)));

%% Generate MRU measurements with noise
% Preallocate measurement arrays
eta_mru = zeros(N_mru, 6);
nu_mru = zeros(N_mru, 6);

% Preallocate bias arrays
bias_eta_log = zeros(N_mru, 6);
bias_nu_log = zeros(N_mru, 6);

% Set random seed for reproducibility
rng(42);

% Generate measurements for each time step
for k = 1:N_mru
    if k > 1
        bias_eta = bias_eta + sqrt(h_mru) * drift_rate_eta .* randn(6, 1);
        bias_nu = bias_nu + sqrt(h_mru) * drift_rate_nu .* randn(6, 1);
    end

    % Add white Gaussian noise to position measurements
    noise_eta = sigma_eta .* randn(6, 1);
    eta_mru(k, :) = eta_true(k, :)' + noise_eta + bias_eta;

    % Add white Gaussian noise to velocity measurements
    noise_nu = sigma_nu .* randn(6, 1);
    nu_mru(k, :) = nu_true(k, :)' + noise_nu + bias_nu;

    % Log bias for analysis
    bias_eta_log(k, :) = bias_eta';
    bias_nu_log(k, :) = bias_nu';
end

%% Calculate measurement errors
eta_error = eta_mru - eta_true;
nu_error = nu_mru - nu_true;

% Calculate RMS errors
eta_rms = sqrt(mean(eta_error.^2, 1));
nu_rms = sqrt(mean(nu_error.^2, 1));

fprintf('\nMeasurement Error Statistics (RMS):\n');
fprintf('Position errors:\n');
fprintf('  Surge: %.2f cm, Sway: %.2f cm, Heave: %.2f cm\n', ...
    eta_rms(1)*100, eta_rms(2)*100, eta_rms(3)*100);
fprintf('  Roll: %.3f deg, Pitch: %.3f deg, Yaw: %.3f deg\n', ...
    rad2deg(eta_rms(4)), rad2deg(eta_rms(5)), rad2deg(eta_rms(6)));
fprintf('Velocity errors:\n');
fprintf('  u: %.2f cm/s, v: %.2f cm/s, w: %.2f cm/s\n', ...
    nu_rms(1)*100, nu_rms(2)*100, nu_rms(3)*100);
fprintf('  p: %.3f deg/s, q: %.3f deg/s, r: %.3f deg/s\n', ...
    rad2deg(nu_rms(4)), rad2deg(nu_rms(5)), rad2deg(nu_rms(6)));

%% Save MRU data for use in motion compensation simulations
mru_data.time = t_mru;
mru_data.eta_true = eta_true;
mru_data.nu_true = nu_true;
mru_data.eta_mru = eta_mru;
mru_data.nu_mru = nu_mru;
mru_data.eta_error = eta_error;
mru_data.nu_error = nu_error;
mru_data.sigma_eta = sigma_eta;
mru_data.sigma_nu = sigma_nu;
mru_data.description = 'MRU measurements with simulated noise';

save('vessel_mru_data.mat', 'mru_data');
fprintf('\nMRU data saved to vessel_mru_data.mat\n');

% %% Signal processing using Kalman Filter
% f_m = 1/h_mru;
% f_s = 100; 
% 
% Z = f_s / f_m;
% if (mod(Z, 1) ~= 0 || Z < 1)
%     error('f_s must be specified such that Z = f_s/f_m is an integer >= 1. Current Z = %.4f', Z);
% end
% 
% h_s = 1/f_s;
% 
% % Model state space equations
% J = eulerang(eta(4), eta(5), eta(6));
% A = [zeros(6), J, zeros(6);
%     zeros(6), -Minv*D, Minv*J;
%     zeros(6), zeros(6), zeros(6)];
% 
% B = [zeros(6); Minv*B_out; zeros(6)];
% 
% C = [eye(6), zeros(6), zeros(6)];
% 
% E = [zeros(6); Minv; zeors(6)];
% 
% 
% % Discrete matrices
% Ad = eye(6) + h_s * A;
% Bd = h_s * B;
% Cd = C;
% Ed = h_s * E;



%% Visualization
fprintf('\n=== GENERATING PLOTS ===\n');

%% === VESSEL MOTION RAO PLOTS (from Ship_motion_RAO.m) ===

% Time-series (discard initial transient)
startIndex = max(1, floor(T_initTransient / h) + 1);
tt_vessel = t(startIndex:end) - t(startIndex);

% Unpack motion data
eta_WF    = simdata_mRAO(startIndex:end, 1:6);      % positions
nu_WF     = simdata_mRAO(startIndex:end, 7:12);     % velocities
nudot_WF  = simdata_mRAO(startIndex:end, 13:18);    % accelerations
waveElevM = simdata_mRAO(startIndex:end, 19);       % wave elevation (from motion RAO)

% Unpack force data
tau_wave  = simdata_fRAO(startIndex:end, 1:6);  % 6 DOF forces/moments
waveElevF = simdata_fRAO(startIndex:end, 7);    % wave elevation (from force RAO)

% Extract LF motion from x_log
eta_LF = x_log(startIndex:end, 1:6);      % Low-frequency positions
nu_LF  = x_log(startIndex:end, 7:12);     % Low-frequency velocities
tau_ctrl = tau_ctrl_log(startIndex:end, :);

% ---- Figure 101: Wave spectrum and wave elevation (Motion RAO) ----
figure(101); clf;
set(gcf, 'Name', 'Motion RAO: Spectrum and Wave Elevation', 'NumberTitle', 'off');

subplot(2,1,1); hold on;
if spreadingFlag
    % Show a few representative directions
    midIdx  = max(1, floor(length(mu)/2));
    qtrIdx  = max(1, floor(length(mu)/4));
    endIdx  = length(mu);

    plot(Omega, S_M(:, midIdx), 'LineWidth', 2);
    plot(Omega, S_M(:, qtrIdx), 'LineWidth', 2);
    plot(Omega, S_M(:, endIdx), 'LineWidth', 2);

    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 1.5);
    legend(sprintf('\\mu = %.0f°', rad2deg(mu(midIdx))), ...
           sprintf('\\mu = %.0f°', rad2deg(mu(qtrIdx))), ...
           sprintf('\\mu = %.0f°', rad2deg(mu(endIdx))), ...
           sprintf('\\omega_0 = %.3g rad/s', w0), ...
           'Location','best');
else
    plot(Omega, S_M(:,1), 'LineWidth', 2);
    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 1.5);
    legend('S(\Omega)', sprintf('\\omega_0 = %.3g rad/s', w0), 'Location','best');
end
xlabel('\Omega (rad/s)'); ylabel('m^2 s');
title([spectrumType, ' spectrum']); grid on; hold off;

subplot(2,1,2);
plot(tt_vessel, waveElevM, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('m'); grid on;
title(sprintf('Wave Elevation for \\beta = %.0f° and H_s = %.1f m', rad2deg(beta), Hs));

% ---- Figure 102: Wave-frequency positions (6-DOF) ----
figure(102); clf;
set(gcf, 'Name', 'Motion RAO: Wave-Frequency Positions', 'NumberTitle', 'off');
DOF_txt = {'x-position (m)', 'y-position (m)', 'z-position (m)', ...
           'Roll angle (deg)', 'Pitch angle (deg)', 'Yaw angle (deg)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, T_scale(k)*eta_WF(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Wave-frequency positions (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

% ---- Figure 103: Wave-frequency velocities (6-DOF) ----
figure(103); clf;
set(gcf, 'Name', 'Motion RAO: Wave-Frequency Velocities', 'NumberTitle', 'off');
DOF_txt_v = {'Surge vel (m/s)','Sway vel (m/s)','Heave vel (m/s)', ...
             'Roll rate (deg/s)','Pitch rate (deg/s)','Yaw rate (deg/s)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, T_scale(k)*nu_WF(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_v{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Wave-frequency velocities (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

% ---- Figure 104: Wave-frequency accelerations (6-DOF) ----
figure(104); clf;
set(gcf, 'Name', 'Motion RAO: Wave-Frequency Accelerations', 'NumberTitle', 'off');
DOF_txt_a = {'Surge acc (m/s^2)','Sway acc (m/s^2)','Heave acc (m/s^2)', ...
             'Roll accel (deg/s^2)','Pitch accel (deg/s^2)','Yaw accel (deg/s^2)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, T_scale(k)*nudot_WF(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_a{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Wave-frequency accelerations (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

% ---- Figure 201: Force RAO wave spectrum and elevation ----
figure(201); clf;
set(gcf, 'Name', 'Force RAO: Spectrum and Wave Elevation', 'NumberTitle', 'off');

subplot(2,1,1); hold on;
if spreadingFlag
    midIdx  = max(1, floor(length(mu)/2));
    qtrIdx  = max(1, floor(length(mu)/4));
    endIdx  = length(mu);

    plot(Omega, S_M(:, midIdx), 'LineWidth', 2);
    plot(Omega, S_M(:, qtrIdx), 'LineWidth', 2);
    plot(Omega, S_M(:, endIdx), 'LineWidth', 2);

    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 1.5);
    legend(sprintf('\\mu = %.0f°', rad2deg(mu(midIdx))), ...
           sprintf('\\mu = %.0f°', rad2deg(mu(qtrIdx))), ...
           sprintf('\\mu = %.0f°', rad2deg(mu(endIdx))), ...
           sprintf('\\omega_0 = %.3g rad/s', w0), ...
           'Location','best');
else
    plot(Omega, S_M(:,1), 'LineWidth', 2);
    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 1.5);
    legend('S(\Omega)', sprintf('\\omega_0 = %.3g rad/s', w0), 'Location','best');
end
xlabel('\Omega (rad/s)'); ylabel('m^2 s');
title([spectrumType, ' spectrum']); grid on; hold off;

subplot(2,1,2);
plot(tt_vessel, waveElevF, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('m'); grid on;
title(sprintf('Wave Elevation for \\beta = %.0f° and H_s = %.1f m', rad2deg(beta), Hs));

% ---- Figure 202: 6-DOF generalized 1st-order wave forces ----
figure(202); clf;
set(gcf, 'Name', 'Force RAO: Wave Forces', 'NumberTitle', 'off');
DOF_txt_tau = {'Surge (N)','Sway (N)','Heave (N)','Roll (N·m)','Pitch (N·m)','Yaw (N·m)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, tau_wave(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_tau{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Generalized 1st-order Wave Forces (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

% ---- Figure 301: LF Positions (6-DOF) ----
figure(301); clf;
set(gcf, 'Name', 'Low-Frequency Ship Motion: Positions', 'NumberTitle', 'off');
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, T_scale(k)*eta_LF(:,k), 'LineWidth', 1.8, 'Color', [0 0.4470 0.7410]);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt{k});
end
if exist('sgtitle','file')
    sgtitle(sprintf('Low-Frequency Ship Motion (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));
end

% ---- Figure 302: LF Velocities (6-DOF) ----
figure(302); clf;
set(gcf, 'Name', 'Low-Frequency Ship Motion: Velocities', 'NumberTitle', 'off');
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, T_scale(k)*nu_LF(:,k), 'LineWidth', 1.8, 'Color', [0 0.4470 0.7410]);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_v{k});
end
if exist('sgtitle','file')
    sgtitle(sprintf('Low-Frequency Ship Velocities (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));
end

% ---- Figure 303: Control forces and moments ----
figure(303); clf;
set(gcf, 'Name', 'DP Control Effort', 'NumberTitle', 'off');
for k = 1:6
    subplot(6,1,k);
    plot(tt_vessel, tau_ctrl(:,k), 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980]);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_tau{k});
end
if exist('sgtitle','file')
    sgtitle(sprintf('DP Control Effort (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));
end

% ---- Figure 304: XY Position Plot ----
figure(304); clf;
set(gcf, 'Name', 'Vessel XY Trajectory', 'NumberTitle', 'off');
plot(eta_LF(:,1), eta_LF(:,2), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
hold on;
plot(eta_LF(1,1), eta_LF(1,2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot(eta_LF(end,1), eta_LF(end,2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(0, 0, 'kx', 'MarkerSize', 15, 'LineWidth', 3);
grid on; axis equal;
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Vessel Trajectory', 'Start', 'End', 'Setpoint', 'Location', 'best');
title(sprintf('Vessel XY Position (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));

% ---- Figure 305: Control Effort Statistics ----
figure(305); clf;
set(gcf, 'Name', 'Control Effort Statistics', 'NumberTitle', 'off');

% Calculate RMS values for each DOF
tau_rms = sqrt(mean(tau_ctrl.^2, 1));
tau_max = max(abs(tau_ctrl), [], 1);
tau_mean = mean(abs(tau_ctrl), 1);

subplot(3,1,1);
bar(tau_rms);
set(gca, 'XTickLabel', {'Surge','Sway','Heave','Roll','Pitch','Yaw'});
ylabel('RMS'); title('RMS Control Effort'); grid on;

subplot(3,1,2);
bar(tau_max);
set(gca, 'XTickLabel', {'Surge','Sway','Heave','Roll','Pitch','Yaw'});
ylabel('Max Magnitude'); title('Maximum Control Effort'); grid on;

subplot(3,1,3);
bar(tau_mean);
set(gca, 'XTickLabel', {'Surge','Sway','Heave','Roll','Pitch','Yaw'});
ylabel('Mean Magnitude'); title('Mean Absolute Control Effort'); grid on;

if exist('sgtitle','file')
    sgtitle('Control Effort Statistics');
end

% ---- Figure 306: Position Error Statistics ----
figure(306); clf;
set(gcf, 'Name', 'Position Error Statistics', 'NumberTitle', 'off');

% Calculate position errors (assuming setpoint is at origin)
pos_error = eta_LF;  % Since eta_d = [0 0 0 0 0 0]
pos_error(:,4:6) = pos_error(:,4:6) * 180/pi;  % Convert angles to degrees

pos_rms = sqrt(mean(pos_error.^2, 1));
pos_max = max(abs(pos_error), [], 1);
pos_mean = mean(abs(pos_error), 1);

subplot(3,1,1);
bar(pos_rms);
set(gca, 'XTickLabel', {'X (m)','Y (m)','Z (m)','Roll (°)','Pitch (°)','Yaw (°)'});
ylabel('RMS'); title('RMS Position Error'); grid on;

subplot(3,1,2);
bar(pos_max);
set(gca, 'XTickLabel', {'X (m)','Y (m)','Z (m)','Roll (°)','Pitch (°)','Yaw (°)'});
ylabel('Max'); title('Maximum Position Error'); grid on;

subplot(3,1,3);
bar(pos_mean);
set(gca, 'XTickLabel', {'X (m)','Y (m)','Z (m)','Roll (°)','Pitch (°)','Yaw (°)'});
ylabel('Mean'); title('Mean Absolute Position Error'); grid on;

if exist('sgtitle','file')
    sgtitle('Position Error Statistics');
end

% Print summary statistics
fprintf('\n=== CONTROL PERFORMANCE SUMMARY ===\n');
fprintf('Position Errors (RMS):\n');
fprintf('  X: %.3f m,  Y: %.3f m,  Z: %.3f m\n', pos_rms(1), pos_rms(2), pos_rms(3));
fprintf('  Roll: %.3f°,  Pitch: %.3f°,  Yaw: %.3f°\n', pos_rms(4), pos_rms(5), pos_rms(6));
fprintf('\nControl Effort (RMS):\n');
fprintf('  Surge: %.1f N,  Sway: %.1f N,  Heave: %.1f N\n', tau_rms(1), tau_rms(2), tau_rms(3));
fprintf('  Roll: %.1f N·m,  Pitch: %.1f N·m,  Yaw: %.1f N·m\n', tau_rms(4), tau_rms(5), tau_rms(6));

%% === MRU NOISE PLOTS (from mru_noise_sim.m) ===

% Create time vector relative to start
tt = t_mru - t_mru(1);

% Labels for plots
eta_labels = {'Surge [m]', 'Sway [m]', 'Heave [m]', 'Roll [deg]', 'Pitch [deg]', 'Yaw [deg]'};
nu_labels = {'u [m/s]', 'v [m/s]', 'w [m/s]', 'p [deg/s]', 'q [deg/s]', 'r [deg/s]'};
scale_eta = [1, 1, 1, 180/pi, 180/pi, 180/pi];
scale_nu = [1, 1, 1, 180/pi, 180/pi, 180/pi];

% Figure 1: Position measurements with noise
figure(1); clf;
set(gcf, 'Name', 'MRU Position Measurements', 'NumberTitle', 'off');
for i = 1:6
    subplot(3, 2, i);
    hold on; grid on;
    plot(tt, scale_eta(i) * eta_true(:, i), 'b-', 'LineWidth', 1.5, 'DisplayName', 'True');
    plot(tt, scale_eta(i) * eta_mru(:, i), 'r--', 'LineWidth', 1.0, 'DisplayName', 'MRU Measured');
    xlabel('Time [s]');
    ylabel(eta_labels{i});
    legend('Location', 'best');
    title(sprintf('%s (RMS error: %.3f)', eta_labels{i}, eta_rms(i)*scale_eta(i)));
end
if exist('sgtitle', 'file')
    sgtitle('MRU Position Measurements vs True Motion');
end

% Figure 2: Velocity measurements with noise
figure(2); clf;
set(gcf, 'Name', 'MRU Velocity Measurements', 'NumberTitle', 'off');
for i = 1:6
    subplot(3, 2, i);
    hold on; grid on;
    plot(tt, scale_nu(i) * nu_true(:, i), 'b-', 'LineWidth', 1.5, 'DisplayName', 'True');
    plot(tt, scale_nu(i) * nu_mru(:, i), 'r--', 'LineWidth', 1.0, 'DisplayName', 'MRU Measured');
    xlabel('Time [s]');
    ylabel(nu_labels{i});
    legend('Location', 'best');
    title(sprintf('%s (RMS error: %.3f)', nu_labels{i}, nu_rms(i)*scale_nu(i)));
end
if exist('sgtitle', 'file')
    sgtitle('MRU Velocity Measurements vs True Motion');
end

% Figure 3: Measurement errors
figure(3); clf;
set(gcf, 'Name', 'MRU Measurement Errors', 'NumberTitle', 'off');
for i = 1:6
    subplot(3, 2, i);
    hold on; grid on;
    plot(tt, scale_eta(i) * eta_error(:, i), 'r-', 'LineWidth', 1.0);
    plot(tt([1 end]), scale_eta(i) * [1 1] * sigma_eta(i), 'k--', 'LineWidth', 1.5, 'DisplayName', '1-\sigma');
    plot(tt([1 end]), -scale_eta(i) * [1 1] * sigma_eta(i), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Time [s]');
    ylabel(['Error ' eta_labels{i}]);
    legend('Location', 'best');
    title(sprintf('Position Error (RMS: %.3f)', eta_rms(i)*scale_eta(i)));
    ylim([-4*sigma_eta(i)*scale_eta(i), 4*sigma_eta(i)*scale_eta(i)]);
end
if exist('sgtitle', 'file')
    sgtitle('MRU Position Measurement Errors');
end

% Figure 4: Velocity measurement errors
figure(4); clf;
set(gcf, 'Name', 'MRU Velocity Errors', 'NumberTitle', 'off');
for i = 1:6
    subplot(3, 2, i);
    hold on; grid on;
    plot(tt, scale_nu(i) * nu_error(:, i), 'r-', 'LineWidth', 1.0);
    plot(tt([1 end]), scale_nu(i) * [1 1] * sigma_nu(i), 'k--', 'LineWidth', 1.5, 'DisplayName', '1-\sigma');
    plot(tt([1 end]), -scale_nu(i) * [1 1] * sigma_nu(i), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Time [s]');
    ylabel(['Error ' nu_labels{i}]);
    legend('Location', 'best');
    title(sprintf('Velocity Error (RMS: %.3f)', nu_rms(i)*scale_nu(i)));
    ylim([-4*sigma_nu(i)*scale_nu(i), 4*sigma_nu(i)*scale_nu(i)]);
end
if exist('sgtitle', 'file')
    sgtitle('MRU Velocity Measurement Errors');
end

% Figure 5: Focus on heave (most important for MCE)
figure(5); clf;
set(gcf, 'Name', 'Heave Motion Detail', 'NumberTitle', 'off');
subplot(2, 1, 1);
hold on; grid on;
plot(tt, eta_true(:, 3), 'b-', 'LineWidth', 2, 'DisplayName', 'True Heave');
plot(tt, eta_mru(:, 3), 'r--', 'LineWidth', 1.5, 'DisplayName', 'MRU Measured');
xlabel('Time [s]');
ylabel('Heave [m]');
legend('Location', 'best');
title(sprintf('Heave Position (MRU accuracy: %.1f cm RMS)', eta_rms(3)*100));

subplot(2, 1, 2);
hold on; grid on;
plot(tt, eta_error(:, 3)*100, 'r-', 'LineWidth', 1.0, 'DisplayName', 'Measurement Error');
plot(tt([1 end]), [1 1] * sigma_eta(3)*100, 'k--', 'LineWidth', 1.5, 'DisplayName', '1-\sigma bound');
plot(tt([1 end]), -[1 1] * sigma_eta(3)*100, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time [s]');
ylabel('Heave Error [cm]');
legend('Location', 'best');
title(sprintf('Heave Measurement Error (RMS: %.2f cm)', eta_rms(3)*100));
ylim([-4*sigma_eta(3)*100, 4*sigma_eta(3)*100]);

if exist('sgtitle', 'file')
    sgtitle('Heave Motion');
end

% Figure 6: Error histograms to verify Gaussian distribution
figure(6); clf;
set(gcf, 'Name', 'MRU Error Distributions', 'NumberTitle', 'off');
for i = 1:6
    subplot(3, 2, i);
    histogram(scale_eta(i) * eta_error(:, i), 30, 'Normalization', 'pdf');
    hold on; grid on;

    % Overlay theoretical Gaussian (scaled to match units)
    sigma_scaled = scale_eta(i) * sigma_eta(i);  % Convert sigma to same units as plotted data
    x_range = linspace(-4*sigma_scaled, 4*sigma_scaled, 100);
    pdf_gaussian = (1/(sigma_scaled*sqrt(2*pi))) * exp(-0.5*(x_range/sigma_scaled).^2);
    plot(x_range, pdf_gaussian, 'r-', 'LineWidth', 2);

    xlabel(['Error ' eta_labels{i}]);
    ylabel('Probability Density');
    title(eta_labels{i});
    legend('Measured', 'Theoretical', 'Location', 'best');
end
if exist('sgtitle', 'file')
    sgtitle('Position Error Distributions');
end

fprintf('\n=== SIMULATION COMPLETE ===\n');
fprintf('Generated files:\n');
fprintf('  - vessel_motion_data.mat (true vessel motion)\n');
fprintf('  - vessel_mru_data.mat (MRU measurements with noise)\n');
