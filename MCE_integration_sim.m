%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DP vessel with motion compensated platfrom for MCE operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

rng(1);

% Load vessel
load supply;

% Simulation Flags
spreadingFlag = true;
plotFlag = false;
forceRaoFlag = true;
useIntegralAction = false;

% Simulation parameters
h = 0.1;
T_final = 300;
T_initTransient = 0;

t = 0:h:T_final+T_initTransient-1;
N = numel(t);

% DP control objectiv
eta_d = [0; 0; 0; 0; 0; 0];
nu_d = [0; 0; 0; 0; 0; 0];
x_d = [eta_d; nu_d];

psi = eta_d(6);
U = sqrt(nu_d(1)^2 + nu_d(2)^2);

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
numDirections = 24;

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

ship_dynamics = @(x, tau) [J(x(1:6)) * x(7:12); 
                           -Minv*(D*x(7:12) + G*x(1:6)) + Minv*tau]; 

%% DP controller 
S = [1 1 0 0 0 1]';

% Kp = diag([1 1 0 0 0 1]);
% Ki = diag([1 1 0 0 0 1]);
% Kd = diag([1 1 0 0 0 1]);

% Computing PID-gains using Algorithem 15.2 from (Fossen, 2021)
omega_b1 = 0.08; omega_b2 = 0.08; omega_b6 = 0.12;
omega_b = [omega_b1, omega_b2, 0, 0, 0, omega_b6];
Omega_b = diag(omega_b);

zeta_pid1 = 1; zeta_pid2 = 1; zeta_pid6 = 1;
zeta_pid = [zeta_pid1, zeta_pid2, 0, 0, 0, zeta_pid6];
Zeta_pid = diag(zeta_pid);

omega_n = zeros(1,6);
for i = 1:length(omega_n)
    omega_n(i) = omega_b(i) / ( sqrt(1 - 2*zeta_pid(i)^2 + sqrt(4*zeta_pid(i)^4 - 4*zeta_pid(i)^2 + 2) ) );
end
Omega_n = diag(omega_n);

% Assuming roll, pitch and yaw is small => J_Theta(eta) = I
% Kp = M*Omega_n^2;
% Kd = 2.*M*Zeta_pid*Omega_n; % - D; 
% Ki = 0.10*Kp*Omega_n;

% PID controller
tau_pid = @(eta, nu, eta_int, Kp, Ki, Kd) S.*(eulerang(eta(4), eta(5), eta(6))'*(-Kp*eta ...
    - Kd*eulerang(eta(4), eta(5), eta(6))*nu - (useIntegralAction*Ki*eta_int)));

% Heading lowpass filter
T_psi = 12;
alpha = h/(T_psi + h);
psi_lp = 0;

%% Motion compensated platform configuration

% Motion compensation mode
mc_mode = 'stabilize';

% platform base location  {b}
r_0b = [-30; 0; -3];

% Joints limits
q1_lim = deg2rad([-15, 15]); % Roll limits
q2_lim = deg2rad([-10, 10]); % Pitch limits
d3_lim = [3, 6];  % Cylinder extension limits

% Fixed rotation of the manipulator w.r.t. platform base
%R_





%% Main simulation loop
% Initial vessel state
eta = zeros(6,1);       % eta = [x y z phi theta psi]
nu = zeros(6,1);        % nu = [u v w p q r]
x = [eta; nu];

eta_int = zeros(6,1);

% Preallocate log data
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

for k = 1:N
    tk = t(k);
    eta = x(1:6);
    nu = x(7:12);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DP Control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Lowpass filtering heading
    psi_err = eta(6) - psi_lp;
    psi_lp = psi_lp + alpha*psi_err;
    psi_lp = ssa(psi_lp);
    
    % Add DP control on LF states
    Jinv = eulerang(eta(4), eta(5), eta(6)) \ eye(6);
    M_star = Jinv' * M * Jinv;

    Kp = M_star*Omega_n^2;
    Kd = 2.*M_star*Zeta_pid*Omega_n; % - D; 
    Ki = 0.10*Kp*Omega_n;

    tau_control = tau_pid(eta, nu, eta_int, Kp, Ki, Kd);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Manipulator motion compensation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Rotation matrix from body to NED
    R_nb = Rzyx(eta(4), eta(5), eta(6));

    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update vessel states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
end

% After your main loop, ignore transient
idx0 = max(1, floor(T_initTransient/h) + 1);
zeta = simdata_mRAO(idx0:end,19);

% 1) Should be ~ Hs/4
fprintf('std(zeta)=%.3f m  (expected ~ %.3f m)\n', std(zeta), Hs/4);

% 2) Spectrum zeroth moment should match variance
dOmega = Omega(2)-Omega(1);                % rad/s
if spreadingFlag
    dmu = 2*pi/numDirections;              % *** radians ***
    m0  = sum(S_M(:))*dOmega*dmu;          % m^2
else
    m0  = sum(S_M(:,1))*dOmega;            % m^2
end
fprintf('m0 from S_M = %.3f m^2,  std(zeta)^2 = %.3f m^2\n', m0, var(zeta));

%% === MOTION RAO PLOTS ===
figure(101); clf;

% Time-series (discard initial transient)
startIndex = max(1, floor(T_initTransient / h) + 1);
tt = t(startIndex:end) - t(startIndex);

% Unpack motion data
eta_WF    = simdata_mRAO(startIndex:end, 1:6);      % positions
nu_WF     = simdata_mRAO(startIndex:end, 7:12);     % velocities
nudot_WF  = simdata_mRAO(startIndex:end, 13:18);    % accelerations
waveElevM = simdata_mRAO(startIndex:end, 19);       % wave elevation (from motion RAO)

% ---- Wave spectrum ----
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

% ---- Wave elevation ----
subplot(2,1,2);
plot(tt, waveElevM, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('m'); grid on;
title(sprintf('Wave Elevation for \\beta = %.0f° and H_s = %.1f m', rad2deg(beta), Hs));

% ---- Positions (6-DOF) ----
figure(102); clf;
DOF_txt = {'x-position (m)', 'y-position (m)', 'z-position (m)', ...
           'Roll angle (deg)', 'Pitch angle (deg)', 'Yaw angle (deg)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];
for k = 1:6
    subplot(6,1,k);
    plot(tt, T_scale(k)*eta_WF(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Wave-frequency positions (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

% ---- Velocities (6-DOF) ----
figure(103); clf;
DOF_txt_v = {'Surge vel (m/s)','Sway vel (m/s)','Heave vel (m/s)', ...
             'Roll rate (deg/s)','Pitch rate (deg/s)','Yaw rate (deg/s)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt, T_scale(k)*nu_WF(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_v{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Wave-frequency velocities (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

% ---- Accelerations (6-DOF) ----
figure(104); clf;
DOF_txt_a = {'Surge acc (m/s^2)','Sway acc (m/s^2)','Heave acc (m/s^2)', ...
             'Roll accel (deg/s^2)','Pitch accel (deg/s^2)','Yaw accel (deg/s^2)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt, T_scale(k)*nudot_WF(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_a{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Wave-frequency accelerations (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

%% === FORCE RAO PLOTS ===
figure(201); clf;

% Time-series (discard initial transient)
startIndex = max(1, floor(T_initTransient / h) + 1);
tt = t(startIndex:end) - t(startIndex);

% Unpack force data
tau_wave  = simdata_fRAO(startIndex:end, 1:6);  % 6 DOF forces/moments
waveElevF = simdata_fRAO(startIndex:end, 7);    % wave elevation (from force RAO)

% ---- Wave spectrum ----
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

% ---- Wave elevation ----
subplot(2,1,2);
plot(tt, waveElevF, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('m'); grid on;
title(sprintf('Wave Elevation for \\beta = %.0f° and H_s = %.1f m', rad2deg(beta), Hs));

% ---- 6-DOF generalized 1st-order wave forces ----
figure(202); clf;
DOF_txt_tau = {'Surge (N)','Sway (N)','Heave (N)','Roll (N·m)','Pitch (N·m)','Yaw (N·m)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt, tau_wave(:,k), 'LineWidth', 1.8);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_tau{k});
end
if exist('sgtitle','file'), sgtitle(sprintf('Generalized 1st-order Wave Forces (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs)); end

%% === SHIP MOTION (Low-Frequency) PLOTS ===
figure(301); clf;

% Time-series (discard initial transient)
startIndex = max(1, floor(T_initTransient / h) + 1);
tt = t(startIndex:end) - t(startIndex);

% Extract LF motion from x_log
eta_LF = x_log(startIndex:end, 1:6);      % Low-frequency positions
nu_LF  = x_log(startIndex:end, 7:12);     % Low-frequency velocities

% ---- LF Positions (6-DOF) ----
DOF_txt = {'x-position (m)', 'y-position (m)', 'z-position (m)', ...
           'Roll angle (deg)', 'Pitch angle (deg)', 'Yaw angle (deg)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];

for k = 1:6
    subplot(6,1,k);
    plot(tt, T_scale(k)*eta_LF(:,k), 'LineWidth', 1.8, 'Color', [0 0.4470 0.7410]);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt{k});
end
if exist('sgtitle','file')
    sgtitle(sprintf('Low-Frequency Ship Motion (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));
end

% ---- LF Velocities (6-DOF) ----
figure(302); clf;
DOF_txt_v = {'Surge vel (m/s)','Sway vel (m/s)','Heave vel (m/s)', ...
             'Roll rate (deg/s)','Pitch rate (deg/s)','Yaw rate (deg/s)'};
for k = 1:6
    subplot(6,1,k);
    plot(tt, T_scale(k)*nu_LF(:,k), 'LineWidth', 1.8, 'Color', [0 0.4470 0.7410]);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_v{k});
end
if exist('sgtitle','file')
    sgtitle(sprintf('Low-Frequency Ship Velocities (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));
end

%% === CONTROL EFFORT PLOTS ===
figure(303); clf;

% Extract control forces/moments from tau_log
tau_ctrl = tau_ctrl_log(startIndex:end, :);

% ---- Control forces and moments ----
DOF_txt_tau = {'Surge Force (N)','Sway Force (N)','Heave Force (N)', ...
               'Roll Moment (N·m)','Pitch Moment (N·m)','Yaw Moment (N·m)'};

for k = 1:6
    subplot(6,1,k);
    plot(tt, tau_ctrl(:,k), 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980]);
    grid on; xlabel('Time (s)'); ylabel(DOF_txt_tau{k});
end
if exist('sgtitle','file')
    sgtitle(sprintf('DP Control Effort (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));
end

%% === XY POSITION PLOT ===
figure(304); clf;
plot(eta_LF(:,1), eta_LF(:,2), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
hold on;
plot(eta_LF(1,1), eta_LF(1,2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot(eta_LF(end,1), eta_LF(end,2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(0, 0, 'kx', 'MarkerSize', 15, 'LineWidth', 3);
grid on; axis equal;
xlabel('X Position (m)'); ylabel('Y Position (m)');
legend('Vessel Trajectory', 'Start', 'End', 'Setpoint', 'Location', 'best');
title(sprintf('Vessel XY Position (\\beta = %.0f°, H_s = %.1f m)', rad2deg(beta), Hs));

%% === CONTROL EFFORT STATISTICS ===
figure(305); clf;

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

%% === POSITION ERROR STATISTICS ===
figure(306); clf;

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

