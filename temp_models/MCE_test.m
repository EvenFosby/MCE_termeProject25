%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DP vessel with Motion Compensated Platform (MCP)
% Integrated simulation with gimbal-based motion compensation
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
MA = vessel.A(:,:,1); % Added mass
M = MRB + MA; % System inertia matrix

% Hydrostatics
rho = 1025; g = 9.81;
Awp = vessel.main.Lwl * vessel.main.B * 0.8; % Waterplane displacement
GM_T = vessel.main.GM_T;
GM_L = vessel.main.GM_L;

% Linear restoring matrix
G = diag([M(1,1)*0.05^2, ...
          M(2,2)*0.05^2, ...
          rho*g*Awp, ...
          m*g*GM_T, ...
          m*g*GM_L, ...
          M(6,6)*0.05^2]);

% Linear damping
zeta = [1 1 0.20 0.03 0.05 1];
D = diag([2*zeta(1)*sqrt(M(1,1)*G(1,1)), ...
          2*zeta(2)*sqrt(M(2,2)*G(2,2)), ...
          2*zeta(3)*sqrt(M(3,3)*G(3,3)), ...
          2*zeta(4)*sqrt(M(4,4)*G(4,4)), ...
          2*zeta(5)*sqrt(M(5,5)*G(5,5)), ...
          2*zeta(6)*sqrt(M(6,6)*G(6,6))]);

% Kinematic mapping
J6 = @(eta) eulerang(eta(4), eta(5), eta(6));

% Continuous-time state derivatives:
% eta_dot = J(eta)*nu
% nu_dot  = -M\G * eta - M\D * nu + M\tau
Minv = M \ eye(6);
ship_dynamics = @(x, tau) [J6(x(1:6)) * x(7:12); 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MOTION COMPENSATED PLATFORM CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Platform mounting location on vessel (in body frame {b})
% [x_b, y_b, z_b]: positive = forward, starboard, down
r_0b = [30; 0; -8];  % meters - adjust based on your vessel geometry

% Target position in NED frame (what you want to reach/stabilize)
p_target_n = [50; 20; -5];  % [North, East, Down] meters

% Motion compensation control gains
K_rp = 0.90;      % Roll/Pitch compensation gain (0.8-0.95 typical)
K_heave = 0.85;   % Heave compensation gain (0.7-0.9 typical)
K_pos = 0.3;      % Position tracking gain
K_damp = 0.5;     % Velocity damping

% Joint limits
q1_lim = deg2rad([-25, 25]);    % Roll compensation range
q2_lim = deg2rad([-25, 25]);    % Pitch compensation range
d3_lim = [8, 30];                % Extension range [m]

% Nominal extension (when vessel is level)
d3_nominal = 18;  % [m]

% Initial manipulator configuration
q_manip = [0; 0; d3_nominal];       % [q1, q2, d3]
qdot_manip = zeros(3,1);

% Motion compensation mode
mc_mode = 'stabilize';  % 'stabilize': keep end-effector inertially fixed
                         % 'track': follow target position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sea state and wave spectrum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sea state
Hs      = 2.5;               % Significant wave height [m]
gamma   = 3.3; 
beta    = deg2rad(145);     % Wave direction relative to bow [rad]

Tz = 10;            % Zero-crossing period [s]
T0 = Tz / 0.710;    % Wave spectrum modal (peak) period [s] (Fossen 2021, Eq. 10.61)
w0 = 2*pi / T0;     % Wave spectrum modal (peak) frequency [rad/s]

spectrumParam = [Hs, w0, gamma];

maxFreq = 2*pi; % 3.0;                  % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 60;          % Number of wave frequency intervals (>50)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ship motion simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.1;
T_final = 200;
T_initTransient = 20; 

t = (0:dt:T_final+T_initTransient-1);
N = numel(t);

U_ship = 0;

% Low pass filtering psi
Tpsi = 12;
alpha = 1 - exp(-dt/Tpsi);
psi_lp = 0;

% Simulation state log - Vessel
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

% Simulation state log - MCP
q_manip_log     = zeros(N, 3);
p_e_n_log       = zeros(N, 3);
v_e_n_log       = zeros(N, 3);
R_e_n_log       = zeros(3, 3, N);
mc_error_log    = zeros(N, 3);
vessel_rp_log   = zeros(N, 2);  % Vessel roll and pitch for comparison

y_eta = zeros(6,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Starting DP + MCP Simulation...\n');
fprintf('Sea State: Hs = %.1f m, Tp = %.1f s, Wave direction = %.0f deg\n', Hs, T0, rad2deg(beta));
fprintf('MCP Mode: %s\n\n', mc_mode);

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
        omega, mu, vessel, U_ship, psi_lp, beta, numFreqIntervals);

    % Add DP control on LF states
    tau = tau_dp(eta, nu);

    if use_integral
        e_eta = eta_ref - eta;
        eta_int = eta_int + dt * (S * e_eta);
    end

    % Update vessel states using rk4
    x = rk4(ship_dynamics, dt, x, tau);

    % "Measured" total output = LF + WF (per Fossen 2021)
    y_eta   = eta + eta_wf;
    y_nu    = nu  + nu_wf;
    y_nudot = nudot_wf; % LF accel not kept explicitly here
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MOTION COMPENSATED PLATFORM CONTROL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Total vessel motion (what affects the platform)
    eta_total = y_eta;
    nu_total = y_nu;
    
    % Extract vessel roll, pitch, heave
    phi_vessel = eta_total(4);      % Roll
    theta_vessel = eta_total(5);    % Pitch
    z_vessel = eta_total(3);        % Heave (down positive in NED)
    
    p_vessel = nu_total(4);         % Roll rate
    q_vessel = nu_total(5);         % Pitch rate
    w_vessel = nu_total(3);         % Heave rate
    
    % Compute current end-effector state
    [p_e_n, v_e_n, R_e_n] = motionCompensatedPlatform(eta_total, nu_total, ...
        q_manip, qdot_manip, r_0b);
    
    % Motion compensation strategy
    if strcmp(mc_mode, 'stabilize')
        % Goal: Keep end-effector orientation level (gravity-aligned)
        %       and height constant in NED frame
        
        % Gimbal angles should counter vessel roll/pitch
        q1_desired = -K_rp * phi_vessel;      % Counter roll
        q2_desired = -K_rp * theta_vessel;    % Counter pitch
        
        % Extension compensates heave (keep constant NED height)
        % When vessel heaves down (z increases), extend more
        z_offset = z_vessel - eta_0(3);  % Deviation from initial depth
        d3_desired = d3_nominal + K_heave * z_offset;
        
        % Joint velocity commands (simple PD control)
        % Feedforward + feedback compensation
        q1dot_cmd = -5.0 * (q_manip(1) - q1_desired) - K_damp * qdot_manip(1) - K_rp * p_vessel;
        q2dot_cmd = -5.0 * (q_manip(2) - q2_desired) - K_damp * qdot_manip(2) - K_rp * q_vessel;
        d3dot_cmd = -2.0 * (q_manip(3) - d3_desired) - K_damp * qdot_manip(3) - K_heave * w_vessel;
        
    elseif strcmp(mc_mode, 'track')
        % Goal: Track target position while compensating motion
        
        % Position error in NED
        e_pos = p_target_n - p_e_n;
        
        % Desired velocity (tracking + motion compensation)
        v_desired_n = K_pos * e_pos;
        
        % Compute joint velocities using inverse kinematics
        qdot_cmd = inverseVelocityKinematics(eta_total, nu_total, q_manip, ...
            v_desired_n, r_0b, K_rp, K_heave);
        
        q1dot_cmd = qdot_cmd(1);
        q2dot_cmd = qdot_cmd(2);
        d3dot_cmd = qdot_cmd(3);
    end
    
    % Apply rate limits
    qdot_max = [deg2rad(15); deg2rad(15); 3.0];  % [rad/s, rad/s, m/s]
    q1dot_cmd = max(-qdot_max(1), min(qdot_max(1), q1dot_cmd));
    q2dot_cmd = max(-qdot_max(2), min(qdot_max(2), q2dot_cmd));
    d3dot_cmd = max(-qdot_max(3), min(qdot_max(3), d3dot_cmd));
    
    qdot_cmd = [q1dot_cmd; q2dot_cmd; d3dot_cmd];
    
    % Integrate joint positions
    q_manip_new = q_manip + dt * qdot_cmd;
    
    % Apply position limits with saturation
    q_manip_new(1) = max(q1_lim(1), min(q1_lim(2), q_manip_new(1)));
    q_manip_new(2) = max(q2_lim(1), min(q2_lim(2), q_manip_new(2)));
    q_manip_new(3) = max(d3_lim(1), min(d3_lim(2), q_manip_new(3)));
    
    % Update manipulator state
    qdot_manip = (q_manip_new - q_manip) / dt;
    q_manip = q_manip_new;
    
    % Compute tracking error
    if strcmp(mc_mode, 'track')
        mc_error = p_target_n - p_e_n;
    else
        % For stabilization, error is deviation from initial position
        if k == 1
            p_e_n_initial = p_e_n;
        end
        mc_error = p_e_n_initial - p_e_n;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Logging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Vessel states
    x_log(k,:)        = [eta.' nu.'];
    tau_log(k,:)      = tau.';
    eta_wf_log(k,:)   = eta_wf.';
    nu_wf_log(k,:)    = nu_wf.';
    nudot_wf_log(k,:) = nudot_wf.';
    zeta_log(k)       = zeta;
    y_eta_log(k,:)    = y_eta.';
    y_nu_log(k,:)     = y_nu.';
    y_nudot_log(k,:)  = y_nudot.';
    psi_lp_log(k,:)   = psi_lp;
    
    % MCP states
    q_manip_log(k,:)   = q_manip';
    p_e_n_log(k,:)     = p_e_n';
    v_e_n_log(k,:)     = v_e_n';
    R_e_n_log(:,:,k)   = R_e_n;
    mc_error_log(k,:)  = mc_error';
    vessel_rp_log(k,:) = [phi_vessel, theta_vessel]';
    
    % Progress indicator
    if mod(k, floor(N/10)) == 0
        fprintf('Progress: %.0f%%\n', 100*k/N);
    end
end

fprintf('Simulation complete!\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-process: Trim transient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startIndex  = max(1, floor(T_initTransient/dt) + 1);
tt  = t(startIndex:end) - t(startIndex);

% Vessel data
y_eta_log   = y_eta_log(startIndex:end,:);
y_nu_log    = y_nu_log(startIndex:end,:);
eta_wf_log  = eta_wf_log(startIndex:end,:);
nu_wf_log   = nu_wf_log(startIndex:end,:);
zeta_log    = zeta_log(startIndex:end);
tau_log     = tau_log(startIndex:end,:);
psi_lp_log  = psi_lp_log(startIndex:end,:);

% MCP data
q_manip_log    = q_manip_log(startIndex:end,:);
p_e_n_log      = p_e_n_log(startIndex:end,:);
v_e_n_log      = v_e_n_log(startIndex:end,:);
mc_error_log   = mc_error_log(startIndex:end,:);
vessel_rp_log  = vessel_rp_log(startIndex:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots - Vessel Motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
ylabel('\psi_{LP} (deg)');
sgtitle('Measured Vessel positions & Euler angles (LF + WF)');

% 6-DOF measured velocities
figure(3); clf;
DOFtxt_v = {'u (m/s)','v (m/s)','w (m/s)','p (deg/s)','q (deg/s)','r (deg/s)'};
for i = 1:6
    subplot(6,1,i); grid on;
    plot(tt, Tscale(i)*y_nu_log(:,i), 'LineWidth', 1.6);
    ylabel(DOFtxt_v{i});
end
xlabel('Time (s)');
sgtitle('Measured Vessel velocities (LF + WF)');

% Control effort
figure(4); clf; grid on;
plot(tt, tau_log, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('\tau (N / N·m)');
title('DP control effort (LF only)');
legend('X','Y','Z','K','M','N');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots - Motion Compensated Platform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Joint positions with vessel motion comparison
figure(5); clf;

subplot(3,1,1); hold on; grid on;
plot(tt, rad2deg(q_manip_log(:,1)), 'b', 'LineWidth', 1.8);
plot(tt, rad2deg(-vessel_rp_log(:,1)), 'r--', 'LineWidth', 1.2);
ylabel('q_1 (deg)'); 
legend('Gimbal Roll', '-Vessel Roll', 'Location', 'best');
title('Roll Compensation');

subplot(3,1,2); hold on; grid on;
plot(tt, rad2deg(q_manip_log(:,2)), 'b', 'LineWidth', 1.8);
plot(tt, rad2deg(-vessel_rp_log(:,2)), 'r--', 'LineWidth', 1.2);
ylabel('q_2 (deg)');
legend('Gimbal Pitch', '-Vessel Pitch', 'Location', 'best');
title('Pitch Compensation');

subplot(3,1,3); hold on; grid on;
plot(tt, q_manip_log(:,3), 'b', 'LineWidth', 1.8);
plot(tt, d3_nominal*ones(size(tt)), 'r--', 'LineWidth', 1.2);
ylabel('d_3 (m)');
xlabel('Time (s)');
legend('Extension', 'Nominal', 'Location', 'best');
title('Heave Compensation');

sgtitle('Motion Compensation Joint States');

% End-effector position stability
figure(6); clf;

subplot(3,1,1); hold on; grid on;
plot(tt, p_e_n_log(:,1), 'LineWidth', 1.6);
if strcmp(mc_mode, 'track')
    plot(tt, p_target_n(1)*ones(size(tt)), 'r--', 'LineWidth', 1.2);
    legend('Actual', 'Target', 'Location', 'best');
end
ylabel('North (m)'); grid on;
title('End-Effector North Position');

subplot(3,1,2); hold on; grid on;
plot(tt, p_e_n_log(:,2), 'LineWidth', 1.6);
if strcmp(mc_mode, 'track')
    plot(tt, p_target_n(2)*ones(size(tt)), 'r--', 'LineWidth', 1.2);
    legend('Actual', 'Target', 'Location', 'best');
end
ylabel('East (m)'); grid on;
title('End-Effector East Position');

subplot(3,1,3); hold on; grid on;
plot(tt, p_e_n_log(:,3), 'LineWidth', 1.6);
if strcmp(mc_mode, 'track')
    plot(tt, p_target_n(3)*ones(size(tt)), 'r--', 'LineWidth', 1.2);
    legend('Actual', 'Target', 'Location', 'best');
end
ylabel('Down (m)'); xlabel('Time (s)'); grid on;
title('End-Effector Down Position');

sgtitle('End-Effector Position in NED Frame');

% Position error statistics
figure(7); clf;

subplot(2,1,1); hold on; grid on;
plot(tt, mc_error_log(:,1), 'LineWidth', 1.5);
plot(tt, mc_error_log(:,2), 'LineWidth', 1.5);
plot(tt, mc_error_log(:,3), 'LineWidth', 1.5);
ylabel('Error (m)');
legend('North', 'East', 'Down', 'Location', 'best');
title('Position Error Components');
grid on;

subplot(2,1,2); hold on; grid on;
error_norm = sqrt(sum(mc_error_log.^2, 2));
plot(tt, error_norm, 'LineWidth', 1.8);
ylabel('|Error| (m)'); xlabel('Time (s)');
title('Total Position Error Magnitude');
grid on;

sgtitle(sprintf('Motion Compensation Performance (Mode: %s)', mc_mode));

% 3D trajectory plot
figure(8); clf; hold on; grid on;
plot3(p_e_n_log(:,1), p_e_n_log(:,2), p_e_n_log(:,3), 'b-', 'LineWidth', 2);
if strcmp(mc_mode, 'track')
    plot3(p_target_n(1), p_target_n(2), p_target_n(3), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
    legend('End-effector trajectory', 'Target', 'Location', 'best');
else
    plot3(p_e_n_log(1,1), p_e_n_log(1,2), p_e_n_log(1,3), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    legend('End-effector trajectory', 'Initial position', 'Location', 'best');
end
xlabel('North (m)'); ylabel('East (m)'); zlabel('Down (m)');
title('3D End-Effector Motion in NED Frame');
axis equal; view(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performance Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('=== Motion Compensation Performance ===\n');
fprintf('Position Error RMS:  [%.3f, %.3f, %.3f] m\n', rms(mc_error_log));
fprintf('Position Error Max:  [%.3f, %.3f, %.3f] m\n', max(abs(mc_error_log)));
fprintf('Position Error Mean: [%.3f, %.3f, %.3f] m\n', mean(mc_error_log));
fprintf('Total Error RMS: %.3f m\n', rms(error_norm));
fprintf('Total Error Max: %.3f m\n', max(error_norm));
fprintf('\n');

fprintf('=== Joint Motion Statistics ===\n');
fprintf('Roll joint (q1):   Mean = %.2f°, Std = %.2f°, Range = [%.2f, %.2f]°\n', ...
    rad2deg(mean(q_manip_log(:,1))), rad2deg(std(q_manip_log(:,1))), ...
    rad2deg(min(q_manip_log(:,1))), rad2deg(max(q_manip_log(:,1))));
fprintf('Pitch joint (q2):  Mean = %.2f°, Std = %.2f°, Range = [%.2f, %.2f]°\n', ...
    rad2deg(mean(q_manip_log(:,2))), rad2deg(std(q_manip_log(:,2))), ...
    rad2deg(min(q_manip_log(:,2))), rad2deg(max(q_manip_log(:,2))));
fprintf('Extension (d3):    Mean = %.2f m, Std = %.2f m, Range = [%.2f, %.2f] m\n', ...
    mean(q_manip_log(:,3)), std(q_manip_log(:,3)), ...
    min(q_manip_log(:,3)), max(q_manip_log(:,3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_e_n, v_e_n, R_e_n] = motionCompensatedPlatform(eta_vessel, nu_vessel, ...
    q_manip, qdot_manip, r_0b)
% Forward kinematics for gimbal + prismatic manipulator

    % Extract vessel motion
    x_n = eta_vessel(1:3);
    phi = eta_vessel(4);
    theta = eta_vessel(5);
    psi = eta_vessel(6);
    
    omega_b = nu_vessel(4:6);
    v_b = nu_vessel(1:3);
    
    % Rotation matrix from body to NED
    R_nb = Rzyx(phi, theta, psi);
    
    % Extract manipulator configuration
    q1 = q_manip(1);
    q2 = q_manip(2);
    d3 = q_manip(3);
    
    % Forward kinematics: R(x, q1) -> R(y, q2) -> P(z, d3)
    R_x_q1 = [1   0         0;
              0   cos(q1)  -sin(q1);
              0   sin(q1)   cos(q1)];
    
    R_y_q2 = [cos(q2)   0   sin(q2);
              0         1   0;
             -sin(q2)   0   cos(q2)];
    
    % Combined rotation from base {0} to end-effector
    R_e0 = R_x_q1 * R_y_q2;
    
    % End-effector position relative to platform base {0}
    p_e0 = R_e0 * [0; 0; d3];
    
    % Position in vessel body frame
    p_eb = r_0b + p_e0;
    
    % Position in NED frame
    p_e_n = x_n + R_nb * p_eb;
    
    % Rotation matrix: NED to end-effector
    R_e_n = R_nb * R_e0;
    
    % Velocity kinematics
    J_manip = manipulatorJacobian(q1, q2, d3);
    v_e0_in0 = J_manip * qdot_manip;
    v_e0_inb = R_e0 * v_e0_in0;
    
    % Vessel contribution to end-effector velocity
    v_vessel_contribution = v_b + cross(omega_b, p_eb);
    
    % Total end-effector velocity in body frame
    v_eb = v_vessel_contribution + v_e0_inb;
    
    % Transform to NED frame
    v_e_n = R_nb * v_eb;
end

function J = manipulatorJacobian(q1, q2, d3)
% Geometric Jacobian for gimbal + prismatic

    c1 = cos(q1); s1 = sin(q1);
    c2 = cos(q2); s2 = sin(q2);
    
    % ∂p/∂q1 (rotation about x-axis)
    dp_dq1 = [0;
              -d3*c1*c2;
              -d3*s1*c2];
    
    % ∂p/∂q2 (rotation about y-axis after x rotation)
    dp_dq2 = [d3*c2;
              d3*s1*s2;
              -d3*c1*s2];
    
    % ∂p/∂d3 (translation along rotated z-axis)
    dp_dd3 = [s2;
              -s1*c2;
              c1*c2];
    
    J = [dp_dq1, dp_dq2, dp_dd3];
end

function R = Rzyx(phi, theta, psi)
% Rotation matrix from body to NED using Z-Y-X Euler angles

    cphi = cos(phi); sphi = sin(phi);
    cth = cos(theta); sth = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);
    
    R = [cth*cpsi, -cphi*spsi + sphi*sth*cpsi,  sphi*spsi + cphi*sth*cpsi;
         cth*spsi,  cphi*cpsi + sphi*sth*spsi, -sphi*cpsi + cphi*sth*spsi;
         -sth,      sphi*cth,                   cphi*cth];
end

function qdot = inverseVelocityKinematics(eta_vessel, nu_vessel, q_manip, ...
    v_desired_n, r_0b, K_rp, K_heave)
% Compute joint velocities for desired end-effector velocity
% Includes motion compensation for vessel disturbances

    q1 = q_manip(1);
    q2 = q_manip(2);
    d3 = q_manip(3);
    
    % Vessel motion
    phi = eta_vessel(4);
    theta = eta_vessel(5);
    psi = eta_vessel(6);
    omega_b = nu_vessel(4:6);
    v_b = nu_vessel(1:3);
    
    % Rotation matrices
    R_nb = Rzyx(phi, theta, psi);
    R_x_q1 = [1 0 0; 0 cos(q1) -sin(q1); 0 sin(q1) cos(q1)];
    R_y_q2 = [cos(q2) 0 sin(q2); 0 1 0; -sin(q2) 0 cos(q2)];
    R_e0 = R_x_q1 * R_y_q2;
    
    % End-effector position in body frame
    p_e0 = R_e0 * [0; 0; d3];
    p_eb = r_0b + p_e0;
    
    % Manipulator Jacobian
    J_manip_0 = manipulatorJacobian(q1, q2, d3);
    J_manip_b = R_e0 * J_manip_0;
    J_manip_n = R_nb * J_manip_b;
    
    % Vessel velocity contribution
    v_vessel_n = R_nb * (v_b + cross(omega_b, p_eb));
    
    % Add motion compensation term
    v_comp_n = -[K_rp * v_vessel_n(1:2); K_heave * v_vessel_n(3)];
    
    % Required velocity
    v_rel_n = v_desired_n + v_comp_n - v_vessel_n;
    
    % Damped least-squares inverse
    lambda = 0.01;
    J_dls = J_manip_n' / (J_manip_n * J_manip_n' + lambda^2 * eye(3));
    
    qdot = J_dls * v_rel_n;
end