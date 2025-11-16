%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DP Vessel with Motion Compensated Platform for MCE Operations
% Complete Implementation with Improved Stability
%
% Features:
% - Dynamic Positioning (DP) control of supply vessel
% - RRP manipulator with motion compensation at stern
% - Improved control gains for stable end-effector performance
% - Platform base at stern: [-35m, 0, 4m] in body frame
% - SNAME convention: x forward, y starboard, z down
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

rng(1);

% Load vessel data
load supply;

%% Simulation Configuration

% Simulation Flags
spreadingFlag = true;       % Use directional wave spreading
plotFlag = false;           % Internal plotting flags
forceRaoFlag = false;       % Use motion RAO (not force RAO)
useIntegralAction = false;  % DP controller integral action

% Simulation parameters
h = 0.1;                    % Time step [s]
T_final = 300;              % Simulation duration [s]
T_initTransient = 0;        % Initial transient to discard [s]

t = 0:h:T_final+T_initTransient-1;
N = numel(t);

% Vessel control objective
psi = 0;                    % Desired heading [rad]
U = 0;                      % Desired forward speed [m/s]

%% Sea State and Wave Spectrum

Hs = 2.5;                   % Significant wave height [m]
Tz = 6;                     % Zero-crossing period [s]

T0 = Tz / 0.710;            % Peak period [s]
w0 = 2*pi / T0;             % Peak frequency [rad/s]

gamma = 3.3;                % JONSWAP peak enhancement factor
beta = deg2rad(140);        % Wave direction relative to bow [rad]

spectrumType = 'JONSWAP';
spectrumParam = [Hs, w0, gamma];

maxFreq = 3.0;              % Maximum frequency for RAO [rad/s]
numFreqIntervals = 60;      % Number of frequency intervals
numDirections = 24;         % Number of wave directions

% Trim vessel RAO data to maxFreq
if vessel.forceRAO.w(end) > maxFreq
    w_index = find(vessel.forceRAO.w > maxFreq, 1) - 1;
    vessel.forceRAO.w = vessel.forceRAO.w(1:w_index);
    for DOF = 1:length(vessel.forceRAO.amp)
        vessel.forceRAO.amp{DOF} = vessel.forceRAO.amp{DOF}(1:w_index, :, :);
        vessel.forceRAO.phase{DOF} = vessel.forceRAO.phase{DOF}(1:w_index, :, :);
    end
end

omegaMax = vessel.forceRAO.w(end);

% Generate directional wave spectrum
[S_M, Omega, Amp, ~, ~, mu] = waveDirectionalSpectrum(spectrumType, ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% Vessel Dynamic Model

MRB = vessel.MRB;           % Rigid-body mass matrix

% Compute frequency-dependent added mass and damping
g = 9.81;
omega_p = w0 - (w0^2 / g) * U * cos(beta);  % Encounter frequency
vessel = computeManeuveringModel(vessel, omega_p, plotFlag);

% Extract diagonal elements for linear model
MA = diag([vessel.A_eq(1,1,1,1), ...
           vessel.A_eq(2,2,1,1), ...
           vessel.A_eq(3,3,1,1), ...
           vessel.A_eq(4,4,1,1), ...
           vessel.A_eq(5,5,1,1), ...
           vessel.A_eq(6,6,1,1)]);

M = MRB + MA;               % Total mass matrix
Minv = M \ eye(6);          % Inverse mass matrix

D = diag([vessel.B_eq(1,1,1,1), ...
          vessel.B_eq(2,2,1,1), ...
          vessel.B_eq(3,3,1,1), ...
          vessel.B_eq(4,4,1,1), ...
          vessel.B_eq(5,5,1,1), ...
          vessel.B_eq(6,6,1,1)]);

% Linear restoring matrix
m = vessel.main.m;          % Vessel mass [kg]
rho = 1025;                 % Seawater density [kg/m^3]
g = 9.81;                   % Gravity [m/s^2]
Awp = vessel.main.Lwl * vessel.main.B * 0.8;  % Waterplane area [m^2]
GM_T = vessel.main.GM_T;    % Transverse metacentric height [m]
GM_L = vessel.main.GM_L;    % Longitudinal metacentric height [m]

R33 = rho*g*Awp;            % Heave restoring
R44 = m*g*GM_T;             % Roll restoring
R55 = m*g*GM_L;             % Pitch restoring

G = diag([0, 0, R33, R44, R55, 0]);

% Kinematic transformation matrix
J = @(eta) eulerang(eta(4), eta(5), eta(6));

% Continuous-time state-space model
% eta_dot = J(eta)*nu
% nu_dot  = -M^(-1)*D*nu - M^(-1)*G*eta + M^(-1)*tau
ship_dynamics = @(x, tau) [J(x(1:6)) * x(7:12); 
                           -Minv*(D*x(7:12) + G*x(1:6)) + Minv*tau]; 

%% DP Controller Design

S = [1 1 0 0 0 1]';         % Controlled DOFs [x, y, -, -, -, psi]

% PID gains using Algorithm 15.2 from Fossen (2021)
omega_b1 = 0.05; omega_b2 = 0.05; omega_b6 = 0.03;  % Bandwidth [rad/s]
omega_b = [omega_b1, omega_b2, 0, 0, 0, omega_b6];
Omega_b = diag(omega_b);

zeta_pid1 = 0.80; zeta_pid2 = 0.80; zeta_pid6 = 1;  % Damping ratios
zeta_pid = [zeta_pid1, zeta_pid2, 0, 0, 0, zeta_pid6];
Zeta_pid = diag(zeta_pid);

% Compute natural frequencies
omega_n = zeros(1,6);
for i = 1:length(omega_n)
    omega_n(i) = omega_b(i) / ( sqrt(1 - 2*zeta_pid(i)^2 + ...
        sqrt(4*zeta_pid(i)^4 - 4*zeta_pid(i)^2 + 2) ) );
end
Omega_n = diag(omega_n);

% Controller gains
Kp = M*Omega_n^2;           % Proportional gain
Kd = 2*M*Zeta_pid*Omega_n;  % Derivative gain
Ki = 0.10*Kp*Omega_n;       % Integral gain

% PID control law
tau_pid = @(eta, nu, eta_int) S.*(eulerang(eta(4), eta(5), eta(6))'*(-Kp*eta ...
    - Kd*eulerang(eta(4), eta(5), eta(6))*nu - (useIntegralAction*Ki*eta_int)));

% Heading low-pass filter
T_psi = 12;                 % Filter time constant [s]
alpha = h/(T_psi + h);      % Filter coefficient
psi_lp = 0;                 % Filtered heading

%% Motion Compensated Platform Configuration

% Motion compensation mode
mc_mode = 'stabilize';      % Options: 'stabilize' or 'track'

% Platform base location in body frame {b}
% SNAME convention: x forward, y starboard, z down
% Ship rotates around CO (center of rotation) at origin
r_0b = [-35; 0; 4];         % [x, y, z] in body frame [m]
                            % x = -35m → 35m AFT of CO (at stern)
                            % y = 0m   → On centerline
                            % z = 4m   → 4m below reference plane (on deck)
                            % Extension is UPWARD (negative z direction)

% Target position in NED frame (for 'track' mode)
p_target_n = [0; 0; -10];   % [N, E, D] [m] - 10m above sea level

% Joint limits
q1_lim = deg2rad([-15, 15]);    % Roll gimbal limits [rad]
q2_lim = deg2rad([-10, 10]);    % Pitch gimbal limits [rad]
d3_lim = [3, 6];                % Prismatic extension limits [m]

% Nominal cylinder extension (when vessel is level)
d3_nominal = 4.5;               % [m]

% Initial manipulator configuration
q_manip = [0; 0; d3_nominal];   % [q1, q2, d3]
qdot_manip = [0; 0; 0];         % Initial joint velocities

% Manipulator actuator dynamics (lowpass filters)
% Each actuator is modeled as a first-order system: T*qdot_dot + qdot = qdot_cmd
T_actuator = [0.3; 0.3; 0.5];   % Time constants [s] for [roll, pitch, extension]
                                 % Roll/pitch: faster hydraulic gimbals (~0.3s)
                                 % Extension: slower hydraulic cylinder (~0.5s)

%% Improved Motion Compensation Control Gains

% Feedforward gains (direct compensation of vessel motion)
% Values close to 1.0 provide aggressive compensation
K_rp_ff = 0.98;             % Roll/pitch feedforward gain (0.95-1.0)
K_heave_ff = 0.95;          % Heave feedforward gain (0.90-0.98)

% Feedback position gains (correct positioning errors)
% Higher values = stiffer position control
K_rp_fb = 10.0;             % Roll/pitch position feedback gain
K_heave_fb = 6.0;           % Heave position feedback gain

% Velocity feedback gains (damping to reduce oscillations)
K_rp_vel = 1.5;             % Roll/pitch velocity damping
K_heave_vel = 1.0;          % Heave velocity damping

% Rate limits (maximum joint velocities)
qdot_max = [deg2rad(25);    % Roll rate limit [rad/s]
            deg2rad(25);    % Pitch rate limit [rad/s]
            3.5];           % Extension rate limit [m/s]

% Low-pass filter for commanded velocities
use_velocity_filter = true; % Enable velocity filtering
T_filter = 0.5;             % Filter time constant [s]
alpha_filter = h / (T_filter + h);
qdot_cmd_filtered = [0; 0; 0];

% Position tracking gain (for 'track' mode)
K_pos = 0.5;

% Print configuration
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║     DP VESSEL WITH MOTION COMPENSATED PLATFORM SIMULATION     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('=== VESSEL CONFIGURATION ===\n');
fprintf('Length: %.1f m, Beam: %.1f m, Draft: %.1f m\n', ...
    vessel.main.Lwl, vessel.main.B, vessel.main.T);
fprintf('Mass: %.1f kg\n', vessel.main.m);
fprintf('\n');
fprintf('=== SEA STATE ===\n');
fprintf('Significant wave height: Hs = %.1f m\n', Hs);
fprintf('Peak period: Tp = %.1f s\n', T0);
fprintf('Wave direction: β = %.0f° (relative to bow)\n', rad2deg(beta));
fprintf('Spectrum type: %s (γ = %.1f)\n', spectrumType, gamma);
fprintf('\n');
fprintf('=== MOTION COMPENSATION PLATFORM ===\n');
fprintf('Mode: %s\n', mc_mode);
fprintf('Platform base (body frame): [%.1f, %.1f, %.1f] m\n', r_0b(1), r_0b(2), r_0b(3));
fprintf('Location: %.1f m aft of CO, %.1f m below reference (deck level)\n', abs(r_0b(1)), r_0b(3));
fprintf('Manipulator type: RRP (Roll-Roll-Prismatic gimbal)\n');
fprintf('Nominal extension: %.1f m UPWARD from deck\n', d3_nominal);
fprintf('Actuator time constants: [%.2f, %.2f, %.2f] s (roll, pitch, extension)\n', ...
    T_actuator(1), T_actuator(2), T_actuator(3));
fprintf('\n');
fprintf('=== CONTROL GAINS ===\n');
fprintf('Feedforward:  Roll/Pitch = %.2f, Heave = %.2f\n', K_rp_ff, K_heave_ff);
fprintf('Feedback:     Roll/Pitch = %.1f, Heave = %.1f\n', K_rp_fb, K_heave_fb);
fprintf('Damping:      Roll/Pitch = %.2f, Heave = %.2f\n', K_rp_vel, K_heave_vel);
fprintf('Rate limits:  [%.1f, %.1f, %.1f] deg/s, deg/s, m/s\n', ...
    rad2deg(qdot_max(1)), rad2deg(qdot_max(2)), qdot_max(3));
fprintf('\n');
fprintf('=== SIMULATION ===\n');
fprintf('Duration: %.1f s, Time step: %.2f s\n', T_final, h);
fprintf('Starting simulation...\n\n');

%% Initialize Simulation

% Initial vessel states
eta_0   = [0 0 0 0 0 0]';   % Initial position/orientation
nu_0    = [0 0 0 0 0 0]';   % Initial velocities
x       = [eta_0; nu_0];    % Combined state vector

eta_int = [0 0 0 0 0 0]';   % Integral error (for DP controller)

% Desired vessel states (station-keeping at origin)
eta_d   = [0 0 0 0 0 0]';

% Preallocate data logs - Vessel
x_log           = zeros(N, 12);
tau_log         = zeros(N, 6);
tau_ctrl_log    = zeros(N, 6);
eta_wf_log      = zeros(N, 6);
nu_wf_log       = zeros(N, 6);
nudot_wf_log    = zeros(N, 6);
zeta_log        = zeros(N, 1);
y_eta_log       = zeros(N, 6);
y_nu_log        = zeros(N, 6);
y_nudot_log     = zeros(N, 6);
psi_lp_log      = zeros(N, 1);

% Preallocate data logs - Motion Compensated Platform
q_manip_log      = zeros(N, 3);
qdot_manip_log   = zeros(N, 3);
qdot_cmd_log     = zeros(N, 3);
p_e_n_log        = zeros(N, 3);
v_e_n_log        = zeros(N, 3);
R_e_n_log        = zeros(3, 3, N);
mc_error_log     = zeros(N, 3);
vessel_rp_log    = zeros(N, 2);
vessel_rpvel_log = zeros(N, 2);
limit_hit_log    = zeros(N, 3);

simdata_fRAO = zeros(N, 7);
simdata_mRAO = zeros(N, 19);

y_eta = zeros(6,1);

tic;  % Start timing

%% Main Simulation Loop

for k = 1:N
    tk = t(k);
    eta = x(1:6);
    nu = x(7:12);

    % === Low-pass filter heading ===
    psi_err = eta(6) - psi_lp;
    psi_lp = psi_lp + alpha*psi_err;
    psi_lp = ssa(psi_lp);
    
    % === DP control (Low-Frequency) ===
    tau_control = tau_pid(eta, nu, eta_int);

    if useIntegralAction
        e_eta = eta_d - eta;
        eta_int = eta_int + h * (S .* e_eta);
    end

    % === Wave excitation (Wave-Frequency) ===
    if forceRaoFlag
        % Force RAO approach
        [tau_wave, zeta_wave] = waveForceRAO_v1(tk, S_M, Amp, Omega, mu, ...
            vessel, U, psi, beta, numFreqIntervals);
        simdata_fRAO(k, :) = [tau_wave', zeta_wave'];

        eta_wf = zeros(6,1);
        nu_wf = zeros(6,1);
        nudot_wf = zeros(6,1);
    else
        % Motion RAO approach (superposition)
        [eta_wf, nu_wf, nudot_wf, zeta_wave] = waveMotionRAO_v1(tk, ...
            S_M, Amp, Omega, mu, vessel, U, psi, beta, numFreqIntervals);
        simdata_mRAO(k, :) = [eta_wf', nu_wf', nudot_wf', zeta_wave'];
        
        tau_wave = zeros(6,1);
    end
    
    % Total control forces
    tau = tau_control + tau_wave;

    % === Integrate vessel dynamics ===
    x = rk4(ship_dynamics, h, x, tau);

    % === Measured outputs (LF + WF) ===
    y_eta   = eta + eta_wf;     % Total position/orientation
    y_nu    = nu  + nu_wf;      % Total velocities
    y_nudot = nudot_wf;         % Wave-frequency accelerations

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MOTION COMPENSATED PLATFORM CONTROL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Total vessel motion (measured)
    eta_total = y_eta;
    nu_total = y_nu;
    
    % Extract vessel roll, pitch, heave and their rates
    phi_vessel = eta_total(4);      % Roll [rad]
    theta_vessel = eta_total(5);    % Pitch [rad]
    z_vessel = eta_total(3);        % Heave [m] (down positive in NED)
    
    p_vessel = nu_total(4);         % Roll rate [rad/s]
    q_vessel = nu_total(5);         % Pitch rate [rad/s]
    w_vessel = nu_total(3);         % Heave rate [m/s]
    
    % === Forward kinematics: compute current end-effector state ===
    [p_e_n, v_e_n, R_e_n] = motionCompensatedPlatform(eta_total, nu_total, ...
        q_manip, qdot_manip, r_0b);
    
    % === Motion compensation control law ===
    if strcmp(mc_mode, 'stabilize')
        % GOAL: Keep end-effector orientation level and height constant in NED
        
        % FEEDFORWARD: Desired gimbal angles to counter vessel roll/pitch
        q1_desired = -K_rp_ff * phi_vessel;
        q2_desired = -K_rp_ff * theta_vessel;
        
        % FEEDFORWARD: Extension to maintain constant NED height
        % When vessel heaves down (z increases), RETRACT cylinder (reduce d3)
        z_offset = z_vessel - eta_0(3);
        d3_desired = d3_nominal - K_heave_ff * z_offset;
        
        % Position errors
        e_q1 = q1_desired - q_manip(1);
        e_q2 = q2_desired - q_manip(2);
        e_d3 = d3_desired - q_manip(3);
        
        % CONTROL LAW: Feedforward + Feedback + Damping
        % This structure provides:
        % 1. Fast reaction to vessel motion (feedforward)
        % 2. Correction of steady-state errors (feedback)
        % 3. Oscillation suppression (damping)
        
        q1dot_cmd_raw = K_rp_ff * (-p_vessel) + ...        % FF: cancel vessel roll rate
                        K_rp_fb * e_q1 + ...               % FB: correct position error
                        K_rp_vel * (-qdot_manip(1));       % Damp: reduce oscillations
        
        q2dot_cmd_raw = K_rp_ff * (-q_vessel) + ...        % FF: cancel vessel pitch rate
                        K_rp_fb * e_q2 + ...               % FB: correct position error
                        K_rp_vel * (-qdot_manip(2));       % Damp: reduce oscillations
        
        d3dot_cmd_raw = K_heave_ff * (-w_vessel) + ...     % FF: cancel vessel heave rate
                        K_heave_fb * e_d3 + ...            % FB: correct position error
                        K_heave_vel * (-qdot_manip(3));    % Damp: reduce oscillations
        
    elseif strcmp(mc_mode, 'track')
        % GOAL: Track target position while compensating vessel motion
        
        % Position error in NED frame
        e_pos = p_target_n - p_e_n;
        
        % Desired velocity for position tracking
        v_desired_n = K_pos * e_pos;
        
        % Inverse kinematics with motion compensation
        qdot_cmd_ik = inverseVelocityKinematics(eta_total, nu_total, q_manip, ...
            v_desired_n, r_0b, K_rp_ff, K_heave_ff);
        
        q1dot_cmd_raw = qdot_cmd_ik(1);
        q2dot_cmd_raw = qdot_cmd_ik(2);
        d3dot_cmd_raw = qdot_cmd_ik(3);
    end
    
    % Combine into vector
    qdot_cmd_raw = [q1dot_cmd_raw; q2dot_cmd_raw; d3dot_cmd_raw];
    
    % === Low-pass filter for smooth commands ===
    if use_velocity_filter
        qdot_cmd_filt = qdot_cmd_filtered + alpha_filter * (qdot_cmd_raw - qdot_cmd_filtered);
    else
        qdot_cmd_filt = qdot_cmd_raw;
    end
    
    % === Apply rate limits ===
    qdot_cmd = zeros(3,1);
    limit_hit = zeros(3,1);
    for i = 1:3
        if qdot_cmd_filt(i) > qdot_max(i)
            qdot_cmd(i) = qdot_max(i);
            limit_hit(i) = 1;
        elseif qdot_cmd_filt(i) < -qdot_max(i)
            qdot_cmd(i) = -qdot_max(i);
            limit_hit(i) = -1;
        else
            qdot_cmd(i) = qdot_cmd_filt(i);
        end
    end

    % === Actuator dynamics: First-order lowpass filter ===
    % Model: T*qdot_dot + qdot = qdot_cmd
    % Discrete: qdot(k+1) = qdot(k) + (h/T)*(qdot_cmd(k) - qdot(k))
    qdot_manip_new = zeros(3,1);
    for i = 1:3
        alpha_actuator = h / (T_actuator(i) + h);
        qdot_manip_new(i) = qdot_manip(i) + alpha_actuator * (qdot_cmd(i) - qdot_manip(i));
    end

    % === Integrate joint positions using filtered velocities ===
    q_manip_new = q_manip + h * qdot_manip_new;
    
    % === Apply position limits with saturation ===
    if q_manip_new(1) > q1_lim(2)
        q_manip_new(1) = q1_lim(2);
        limit_hit(1) = 2;
    elseif q_manip_new(1) < q1_lim(1)
        q_manip_new(1) = q1_lim(1);
        limit_hit(1) = -2;
    end
    
    if q_manip_new(2) > q2_lim(2)
        q_manip_new(2) = q2_lim(2);
        limit_hit(2) = 2;
    elseif q_manip_new(2) < q2_lim(1)
        q_manip_new(2) = q2_lim(1);
        limit_hit(2) = -2;
    end
    
    if q_manip_new(3) > d3_lim(2)
        q_manip_new(3) = d3_lim(2);
        limit_hit(3) = 2;
    elseif q_manip_new(3) < d3_lim(1)
        q_manip_new(3) = d3_lim(1);
        limit_hit(3) = -2;
    end
    
    % === Update manipulator state ===
    qdot_manip = qdot_manip_new;  % Use filtered velocity from actuator dynamics
    q_manip = q_manip_new;
    qdot_cmd_filtered = qdot_cmd_filt;
    
    % === Compute tracking error ===
    if strcmp(mc_mode, 'track')
        mc_error = p_target_n - p_e_n;
    else
        % For stabilization mode, error is deviation from initial position
        if k == 1
            p_e_n_initial = p_e_n;
        end
        mc_error = p_e_n_initial - p_e_n;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DATA LOGGING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Vessel states
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
    
    % Motion compensated platform states
    q_manip_log(k,:)      = q_manip';
    qdot_manip_log(k,:)   = qdot_manip';
    qdot_cmd_log(k,:)     = qdot_cmd';
    p_e_n_log(k,:)        = p_e_n';
    v_e_n_log(k,:)        = v_e_n';
    R_e_n_log(:,:,k)      = R_e_n;
    mc_error_log(k,:)     = mc_error';
    vessel_rp_log(k,:)    = [phi_vessel, theta_vessel]';
    vessel_rpvel_log(k,:) = [p_vessel, q_vessel]';
    limit_hit_log(k,:)    = limit_hit';
    
    % Progress indicator
    if mod(k, floor(N/10)) == 0
        fprintf('Progress: %3.0f%% (%.1f s elapsed)\n', 100*k/N, toc);
    end
end

elapsed_time = toc;
fprintf('\nSimulation complete! Total time: %.2f s\n', elapsed_time);
fprintf('Average iteration time: %.3f ms\n\n', 1000*elapsed_time/N);

%% Post-Processing: Remove Initial Transient

idx0 = max(1, floor(T_initTransient/h) + 1);
tt = t(idx0:end) - t(idx0);
zeta = simdata_mRAO(idx0:end,19);

% Wave statistics
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('WAVE STATISTICS\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('std(ζ) = %.3f m  (expected ~ %.3f m = Hs/4)\n', std(zeta), Hs/4);

dOmega = Omega(2)-Omega(1);
if spreadingFlag
    dmu = 2*pi/numDirections;
    m0  = sum(S_M(:))*dOmega*dmu;
else
    m0  = sum(S_M(:,1))*dOmega;
end
fprintf('m₀ from spectrum = %.3f m²,  var(ζ) = %.3f m²\n', m0, var(zeta));
fprintf('\n');

% Trim all data logs
startIndex = idx0;
y_eta_log       = y_eta_log(startIndex:end,:);
y_nu_log        = y_nu_log(startIndex:end,:);
eta_wf_log      = eta_wf_log(startIndex:end,:);
nu_wf_log       = nu_wf_log(startIndex:end,:);
tau_ctrl_log    = tau_ctrl_log(startIndex:end,:);
x_log           = x_log(startIndex:end,:);
psi_lp_log      = psi_lp_log(startIndex:end,:);

q_manip_log      = q_manip_log(startIndex:end,:);
qdot_manip_log   = qdot_manip_log(startIndex:end,:);
qdot_cmd_log     = qdot_cmd_log(startIndex:end,:);
p_e_n_log        = p_e_n_log(startIndex:end,:);
v_e_n_log        = v_e_n_log(startIndex:end,:);
mc_error_log     = mc_error_log(startIndex:end,:);
vessel_rp_log    = vessel_rp_log(startIndex:end,:);
vessel_rpvel_log = vessel_rpvel_log(startIndex:end,:);
limit_hit_log    = limit_hit_log(startIndex:end,:);

%% Performance Statistics

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('VESSEL DP PERFORMANCE\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

eta_LF = x_log(:, 1:6);
pos_error = eta_LF;
pos_error(:,4:6) = pos_error(:,4:6) * 180/pi;
pos_rms = sqrt(mean(pos_error.^2, 1));
pos_max = max(abs(pos_error), [], 1);

fprintf('Position Errors:\n');
fprintf('  RMS:  X = %6.3f m,  Y = %6.3f m,  Z = %6.3f m\n', ...
    pos_rms(1), pos_rms(2), pos_rms(3));
fprintf('        φ = %6.3f°,  θ = %6.3f°,  ψ = %6.3f°\n', ...
    pos_rms(4), pos_rms(5), pos_rms(6));
fprintf('  Max:  X = %6.3f m,  Y = %6.3f m,  Z = %6.3f m\n', ...
    pos_max(1), pos_max(2), pos_max(3));
fprintf('        φ = %6.3f°,  θ = %6.3f°,  ψ = %6.3f°\n', ...
    pos_max(4), pos_max(5), pos_max(6));

tau_rms = sqrt(mean(tau_ctrl_log.^2, 1));
tau_max = max(abs(tau_ctrl_log), [], 1);
fprintf('\nControl Effort:\n');
fprintf('  RMS:  Fx = %7.1f N,  Fy = %7.1f N,  Fz = %7.1f N\n', ...
    tau_rms(1), tau_rms(2), tau_rms(3));
fprintf('        Mx = %7.1f Nm, My = %7.1f Nm, Mz = %7.1f Nm\n', ...
    tau_rms(4), tau_rms(5), tau_rms(6));
fprintf('  Max:  Fx = %7.1f N,  Fy = %7.1f N,  Fz = %7.1f N\n', ...
    tau_max(1), tau_max(2), tau_max(3));
fprintf('        Mx = %7.1f Nm, My = %7.1f Nm, Mz = %7.1f Nm\n', ...
    tau_max(4), tau_max(5), tau_max(6));

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('MOTION COMPENSATION PERFORMANCE\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

error_norm = sqrt(sum(mc_error_log.^2, 2));
fprintf('End-Effector Position Error:\n');
fprintf('  RMS:   North = %6.3f m,  East = %6.3f m,  Down = %6.3f m\n', ...
    rms(mc_error_log));
fprintf('         Total = %6.3f m\n', rms(error_norm));
fprintf('  Max:   North = %6.3f m,  East = %6.3f m,  Down = %6.3f m\n', ...
    max(abs(mc_error_log)));
fprintf('         Total = %6.3f m\n', max(error_norm));
fprintf('  Mean:  North = %6.3f m,  East = %6.3f m,  Down = %6.3f m\n', ...
    mean(mc_error_log));
fprintf('  Std:   North = %6.3f m,  East = %6.3f m,  Down = %6.3f m\n', ...
    std(mc_error_log));

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('JOINT MOTION STATISTICS\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

fprintf('Roll Joint (q₁):\n');
fprintf('  Mean = %6.2f°,  Std = %5.2f°\n', ...
    rad2deg(mean(q_manip_log(:,1))), rad2deg(std(q_manip_log(:,1))));
fprintf('  Range: [%6.2f°, %6.2f°]\n', ...
    rad2deg(min(q_manip_log(:,1))), rad2deg(max(q_manip_log(:,1))));

fprintf('Pitch Joint (q₂):\n');
fprintf('  Mean = %6.2f°,  Std = %5.2f°\n', ...
    rad2deg(mean(q_manip_log(:,2))), rad2deg(std(q_manip_log(:,2))));
fprintf('  Range: [%6.2f°, %6.2f°]\n', ...
    rad2deg(min(q_manip_log(:,2))), rad2deg(max(q_manip_log(:,2))));

fprintf('Extension (d₃):\n');
fprintf('  Mean = %6.2f m,  Std = %5.2f m\n', ...
    mean(q_manip_log(:,3)), std(q_manip_log(:,3)));
fprintf('  Range: [%6.2f m, %6.2f m]\n', ...
    min(q_manip_log(:,3)), max(q_manip_log(:,3)));

% Check for limit saturation
rate_hits = sum(sum(abs(limit_hit_log(:,1:2)) == 1));
pos_hits = sum(sum(abs(limit_hit_log) == 2));
fprintf('\n');
if rate_hits > 0 || pos_hits > 0
    fprintf('⚠ WARNING: Joint limits reached!\n');
    fprintf('  Rate limits hit: %d times (%.1f%% of simulation)\n', ...
        rate_hits, 100*rate_hits/length(tt));
    fprintf('  Position limits hit: %d times (%.1f%% of simulation)\n', ...
        pos_hits, 100*pos_hits/length(tt));
    fprintf('  → Consider: Increasing limits or reducing gains\n');
else
    fprintf('✓ All joints within limits - motion well within operational envelope\n');
end

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('GENERATING PLOTS...\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

%% Visualization: Vessel and Platform Geometry

figure('Name', 'Vessel and Platform Configuration', 'Position', [50 50 1200 700]);

% Get vessel dimensions
Lwl = vessel.main.Lwl;  % Length at waterline
B = vessel.main.B;      % Beam
T = vessel.main.T;      % Draft

% === Top view (XY plane) ===
subplot(2,2,1); hold on; axis equal; grid on;
% Vessel outline
vessel_outline_x = [-Lwl/2, Lwl/2, Lwl/2, -Lwl/2, -Lwl/2];
vessel_outline_y = [-B/2, -B/2, B/2, B/2, -B/2];
fill(vessel_outline_x, vessel_outline_y, [0.7 0.9 1.0], ...
    'EdgeColor', 'b', 'LineWidth', 2);
% CO at origin
plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
% Platform base
plot(r_0b(1), r_0b(2), 'rs', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
% Bow direction arrow
quiver(0, 0, 15, 0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 0.8);
text(8, 3, 'BOW', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x (m) - Forward'); 
ylabel('y (m) - Starboard');
title('Top View');
legend('Vessel', 'CO', 'Platform Base', 'Bow', 'Location', 'best');
xlim([-Lwl/2-10, Lwl/2+10]); 
ylim([-B/2-5, B/2+5]);

% === Side view (XZ plane) ===
subplot(2,2,3); hold on; axis equal; grid on;
% Vessel side profile
vessel_side_x = [-Lwl/2, Lwl/2, Lwl/2, -Lwl/2, -Lwl/2];
vessel_side_z = [0, 0, T, T, 0];
fill(vessel_side_x, vessel_side_z, [0.7 0.9 1.0], ...
    'EdgeColor', 'b', 'LineWidth', 2);
% Waterline
plot([-Lwl/2-10, Lwl/2+10], [0, 0], 'c--', 'LineWidth', 2);
% CO
plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
% Platform base
plot(r_0b(1), r_0b(3), 'rs', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
% Platform with nominal extension (UPWARD, so subtract in z-down convention)
platform_z_nominal = r_0b(3) - d3_nominal;
plot([r_0b(1), r_0b(1)], [r_0b(3), platform_z_nominal], 'r-', 'LineWidth', 4);
plot(r_0b(1), platform_z_nominal, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
xlabel('x (m) - Forward'); 
ylabel('z (m) - Down');
title('Side View');
legend('Vessel', 'Waterline', 'CO', 'Platform Base', 'Extension', 'End-Effector', ...
    'Location', 'northeast');
xlim([-Lwl/2-10, Lwl/2+10]); 
ylim([-3, T+5]);
set(gca, 'YDir', 'reverse');  % z down is positive

% === 3D view ===
subplot(2,2,[2,4]); hold on; axis equal; grid on;

% Vessel hull (simplified box)
vessel_vertices = [
    -Lwl/2, -B/2, 0;   Lwl/2, -B/2, 0;   Lwl/2, B/2, 0;   -Lwl/2, B/2, 0;
    -Lwl/2, -B/2, T;   Lwl/2, -B/2, T;   Lwl/2, B/2, T;   -Lwl/2, B/2, T;
];
vessel_faces = [1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8; 1,2,3,4; 5,6,7,8];
patch('Vertices', vessel_vertices, 'Faces', vessel_faces, ...
    'FaceColor', [0.7, 0.9, 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);

% CO at origin
plot3(0, 0, 0, 'ko', 'MarkerSize', 14, 'MarkerFaceColor', 'k');

% Platform base
plot3(r_0b(1), r_0b(2), r_0b(3), 'rs', ...
    'MarkerSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 2);

% Platform extension (cylinder) - UPWARD extension
[X_cyl, Y_cyl, Z_cyl] = cylinder(0.3, 20);
Z_cyl = Z_cyl * (-d3_nominal) + r_0b(3);  % Negative for upward extension
X_cyl = X_cyl + r_0b(1);
Y_cyl = Y_cyl + r_0b(2);
surf(X_cyl, Y_cyl, Z_cyl, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'none');

% End-effector
plot3(r_0b(1), r_0b(2), platform_z_nominal, 'go', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'g', 'LineWidth', 2);

% Body frame coordinate axes at CO
quiver3(0, 0, 0, 20, 0, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.3);  % x forward
quiver3(0, 0, 0, 0, 12, 0, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.3);  % y starboard
quiver3(0, 0, 0, 0, 0, 12, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.3);  % z down
text(22, 0, 0, 'x (fwd)', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 11);
text(0, 14, 0, 'y (stbd)', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 11);
text(0, 0, 14, 'z (down)', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 11);

xlabel('x (m) - Forward'); 
ylabel('y (m) - Starboard'); 
zlabel('z (m) - Down');
title('3D Configuration (SNAME Convention)');
view(130, 25);
set(gca, 'ZDir', 'reverse');
lighting gouraud;
camlight('headlight');

%% Plot 1: Wave Spectrum and Elevation

figure('Name', 'Wave Conditions', 'Position', [100 100 1000 600]);

% Wave spectrum
subplot(2,1,1); hold on; grid on;
if spreadingFlag
    midIdx  = max(1, floor(length(mu)/2));
    qtrIdx  = max(1, floor(length(mu)/4));
    endIdx  = length(mu);
    
    plot(Omega, S_M(:, midIdx), 'LineWidth', 2, 'DisplayName', sprintf('μ = %.0f°', rad2deg(mu(midIdx))));
    plot(Omega, S_M(:, qtrIdx), 'LineWidth', 2, 'DisplayName', sprintf('μ = %.0f°', rad2deg(mu(qtrIdx))));
    plot(Omega, S_M(:, endIdx), 'LineWidth', 2, 'DisplayName', sprintf('μ = %.0f°', rad2deg(mu(endIdx))));
    
    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('ω₀ = %.3g rad/s', w0));
    legend('Location','best');
else
    plot(Omega, S_M(:,1), 'LineWidth', 2);
    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 1.5);
    legend('S(Ω)', sprintf('ω₀ = %.3g rad/s', w0), 'Location','best');
end
xlabel('Ω (rad/s)'); 
ylabel('S [m² s]');
title(sprintf('%s Wave Spectrum (Hs = %.1f m, Tp = %.1f s)', spectrumType, Hs, T0));

% Wave elevation time series
subplot(2,1,2);
plot(tt, zeta, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
xlabel('Time (s)'); 
ylabel('ζ (m)'); 
grid on;
title(sprintf('Wave Elevation (β = %.0f°, σ = %.3f m)', rad2deg(beta), std(zeta)));
ylim(1.2*[min(zeta) max(zeta)]);

%% Plot 2: Vessel Motion (6-DOF)

figure('Name', 'Vessel Motion', 'Position', [150 50 1000 800]);
DOF_txt = {'x (m)', 'y (m)', 'z (m)', 'φ (deg)', 'θ (deg)', 'ψ (deg)'};
T_scale = [1 1 1 180/pi 180/pi 180/pi];

for k = 1:6
    subplot(7,1,k);
    plot(tt, T_scale(k)*y_eta_log(:,k), 'LineWidth', 1.3, 'Color', [0 0.4470 0.7410]);
    grid on; 
    ylabel(DOF_txt{k});
    if k == 1
        title('Measured Vessel Position & Orientation (LF + WF)');
    end
end
subplot(7,1,7);
plot(tt, T_scale(6)*psi_lp_log(:), 'LineWidth', 1.3, 'Color', [0.8500 0.3250 0.0980]);
xlabel('Time (s)'); 
ylabel('ψ_{LP} (deg)'); 
grid on;
legend('Filtered heading', 'Location', 'best');

%% Plot 3: Motion Compensation Joint States

figure('Name', 'Motion Compensation Joint States', 'Position', [200 100 1100 700]);

% Roll compensation
subplot(3,1,1); hold on; grid on;
plot(tt, rad2deg(q_manip_log(:,1)), 'b', 'LineWidth', 2.5, 'DisplayName', 'Gimbal Roll (q₁)');
plot(tt, rad2deg(-vessel_rp_log(:,1)), 'r--', 'LineWidth', 1.8, 'DisplayName', '-Vessel Roll');
plot(tt, rad2deg(q1_lim(1))*ones(size(tt)), 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(tt, rad2deg(q1_lim(2))*ones(size(tt)), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Limits');
ylabel('Angle (deg)'); 
legend('Location', 'best');
title('Roll Compensation');
ylim(1.3*[min(rad2deg(q1_lim(1)), min(rad2deg(q_manip_log(:,1)))), ...
           max(rad2deg(q1_lim(2)), max(rad2deg(q_manip_log(:,1))))]);

% Pitch compensation
subplot(3,1,2); hold on; grid on;
plot(tt, rad2deg(q_manip_log(:,2)), 'b', 'LineWidth', 2.5, 'DisplayName', 'Gimbal Pitch (q₂)');
plot(tt, rad2deg(-vessel_rp_log(:,2)), 'r--', 'LineWidth', 1.8, 'DisplayName', '-Vessel Pitch');
plot(tt, rad2deg(q2_lim(1))*ones(size(tt)), 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(tt, rad2deg(q2_lim(2))*ones(size(tt)), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Limits');
ylabel('Angle (deg)');
legend('Location', 'best');
title('Pitch Compensation');
ylim(1.3*[min(rad2deg(q2_lim(1)), min(rad2deg(q_manip_log(:,2)))), ...
           max(rad2deg(q2_lim(2)), max(rad2deg(q_manip_log(:,2))))]);

% Heave compensation with tracking
subplot(3,1,3); hold on; grid on;
plot(tt, q_manip_log(:,3), 'b', 'LineWidth', 2.5, 'DisplayName', 'Extension (d₃)');

% Expected extension from feedforward
z_vessel_motion = y_eta_log(:,3);
z_offset_actual = z_vessel_motion - y_eta_log(1,3);
d3_expected = d3_nominal - K_heave_ff * z_offset_actual;
plot(tt, d3_expected, 'g--', 'LineWidth', 1.8, 'DisplayName', 'Expected (FF)');

plot(tt, d3_nominal*ones(size(tt)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Nominal');
plot(tt, d3_lim(1)*ones(size(tt)), 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(tt, d3_lim(2)*ones(size(tt)), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Limits');

ylabel('Extension (m)');
xlabel('Time (s)');
legend('Location', 'best');
title('Heave Compensation (Tracks Vessel Heave)');
ylim(1.2*[min(d3_lim(1), min(q_manip_log(:,3))), ...
           max(d3_lim(2), max(q_manip_log(:,3)))]);

%% Plot 4: Heave Compensation Analysis (Detailed)

figure('Name', 'Heave Compensation Analysis', 'Position', [250 150 1100 700]);

% Vessel heave motion
subplot(3,1,1); hold on; grid on;
plot(tt, y_eta_log(:,3), 'r', 'LineWidth', 2.5, 'DisplayName', 'Vessel Heave');
plot(tt, y_eta_log(1,3)*ones(size(tt)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Initial');
ylabel('z (m)');
legend('Location', 'best');
title('Vessel Vertical Motion (Down Positive)');

% Extension tracking vessel heave
subplot(3,1,2); hold on; grid on;
z_offset_actual = y_eta_log(:,3) - y_eta_log(1,3);
plot(tt, z_offset_actual, 'r', 'LineWidth', 2, 'DisplayName', 'Vessel Heave Δz');
plot(tt, q_manip_log(:,3) - d3_nominal, 'b', 'LineWidth', 2, 'DisplayName', 'Extension Δd₃');
plot(tt, zeros(size(tt)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero');
ylabel('Deviation (m)');
legend('Location', 'best');
title('Compensation Tracking: Extension Should Follow Vessel Heave');

% End-effector vertical position (should be constant)
subplot(3,1,3); hold on; grid on;
plot(tt, p_e_n_log(:,3), 'b', 'LineWidth', 2.5, 'DisplayName', 'Actual');
plot(tt, p_e_n_log(1,3)*ones(size(tt)), 'r--', 'LineWidth', 1.8, 'DisplayName', 'Target');
ylabel('Down (m)');
xlabel('Time (s)');
legend('Location', 'best');
title(sprintf('End-Effector Vertical Position in NED (σ = %.3f m)', std(p_e_n_log(:,3))));
grid on;

%% Plot 5: End-Effector Position

figure('Name', 'End-Effector Position in NED', 'Position', [300 200 1100 700]);

subplot(3,1,1); hold on; grid on;
plot(tt, p_e_n_log(:,1), 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
if strcmp(mc_mode, 'track')
    plot(tt, p_target_n(1)*ones(size(tt)), 'r--', 'LineWidth', 1.8);
    legend('Actual', 'Target', 'Location', 'best');
else
    plot(tt, p_e_n_log(1,1)*ones(size(tt)), 'r--', 'LineWidth', 1.8);
    legend('Actual', 'Initial', 'Location', 'best');
end
ylabel('North (m)');
title(sprintf('End-Effector Position (Mode: %s)', mc_mode));

subplot(3,1,2); hold on; grid on;
plot(tt, p_e_n_log(:,2), 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
if strcmp(mc_mode, 'track')
    plot(tt, p_target_n(2)*ones(size(tt)), 'r--', 'LineWidth', 1.8);
    legend('Actual', 'Target', 'Location', 'best');
else
    plot(tt, p_e_n_log(1,2)*ones(size(tt)), 'r--', 'LineWidth', 1.8);
    legend('Actual', 'Initial', 'Location', 'best');
end
ylabel('East (m)');

subplot(3,1,3); hold on; grid on;
plot(tt, p_e_n_log(:,3), 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
if strcmp(mc_mode, 'track')
    plot(tt, p_target_n(3)*ones(size(tt)), 'r--', 'LineWidth', 1.8);
    legend('Actual', 'Target', 'Location', 'best');
else
    plot(tt, p_e_n_log(1,3)*ones(size(tt)), 'r--', 'LineWidth', 1.8);
    legend('Actual', 'Initial', 'Location', 'best');
end
ylabel('Down (m)'); 
xlabel('Time (s)');

%% Plot 6: Motion Compensation Performance

figure('Name', 'Motion Compensation Performance', 'Position', [350 250 1100 700]);

% Position error components
subplot(3,1,1); hold on; grid on;
plot(tt, mc_error_log(:,1), 'LineWidth', 2, 'DisplayName', 'North');
plot(tt, mc_error_log(:,2), 'LineWidth', 2, 'DisplayName', 'East');
plot(tt, mc_error_log(:,3), 'LineWidth', 2, 'DisplayName', 'Down');
ylabel('Error (m)');
legend('Location', 'best');
title(sprintf('Position Error Components (Mode: %s)', mc_mode));
ylim([-1 1]*max(1.0, max(abs(mc_error_log(:)))*1.2));

% Total error magnitude
subplot(3,1,2); hold on; grid on;
plot(tt, error_norm, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980]);
ylabel('|Error| (m)'); 
title(sprintf('Total Position Error (RMS = %.3f m, Max = %.3f m)', ...
    rms(error_norm), max(error_norm)));
ylim([0 max(1.0, max(error_norm)*1.2)]);
grid on;

% Joint velocities
subplot(3,1,3); hold on; grid on;
plot(tt, rad2deg(qdot_manip_log(:,1)), 'LineWidth', 1.5, 'DisplayName', 'q₁̇ (deg/s)');
plot(tt, rad2deg(qdot_manip_log(:,2)), 'LineWidth', 1.5, 'DisplayName', 'q₂̇ (deg/s)');
plot(tt, qdot_manip_log(:,3), 'LineWidth', 1.5, 'DisplayName', 'd₃̇ (m/s)');
ylabel('Joint Velocities');
legend('Location', 'best');
xlabel('Time (s)');
title('Actual Joint Velocities');

%% Plot 7: 3D Trajectory

figure('Name', '3D End-Effector Trajectory', 'Position', [400 300 900 700]);
hold on; grid on;

% Trajectory
plot3(p_e_n_log(:,1), p_e_n_log(:,2), p_e_n_log(:,3), 'b-', 'LineWidth', 2.5);

% Start and end points
plot3(p_e_n_log(1,1), p_e_n_log(1,2), p_e_n_log(1,3), 'go', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot3(p_e_n_log(end,1), p_e_n_log(end,2), p_e_n_log(end,3), 'rs', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);

if strcmp(mc_mode, 'track')
    plot3(p_target_n(1), p_target_n(2), p_target_n(3), 'r*', ...
        'MarkerSize', 30, 'LineWidth', 3);
    legend('Trajectory', 'Start', 'End', 'Target', 'Location', 'best');
else
    legend('Trajectory', 'Start', 'End', 'Location', 'best');
end

xlabel('North (m)'); 
ylabel('East (m)'); 
zlabel('Down (m)');
title('3D End-Effector Motion in NED Frame');
axis equal; 
view(45, 30);
set(gca, 'ZDir', 'reverse');  % Down is positive

%% Plot 8: Vessel XY Position

figure('Name', 'Vessel Position', 'Position', [450 350 800 600]);
plot(eta_LF(:,1), eta_LF(:,2), 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
hold on;
plot(eta_LF(1,1), eta_LF(1,2), 'go', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot(eta_LF(end,1), eta_LF(end,2), 'rs', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(0, 0, 'kx', 'MarkerSize', 24, 'LineWidth', 4);
grid on; 
axis equal;
xlabel('X Position (m)'); 
ylabel('Y Position (m)');
legend('Trajectory', 'Start', 'End', 'Setpoint', 'Location', 'best');
title(sprintf('Vessel XY Position (β = %.0f°, Hs = %.1f m)', rad2deg(beta), Hs));

%% Summary

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('SIMULATION COMPLETE\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('All figures generated successfully!\n');
fprintf('Total computation time: %.2f seconds\n', elapsed_time);
fprintf('\n');

if rms(error_norm) < 0.5
    fprintf('✓ EXCELLENT: End-effector position very stable (RMS < 0.5 m)\n');
elseif rms(error_norm) < 1.0
    fprintf('✓ GOOD: End-effector position stable (RMS < 1.0 m)\n');
elseif rms(error_norm) < 2.0
    fprintf('⚠ ACCEPTABLE: End-effector stability could be improved (RMS < 2.0 m)\n');
else
    fprintf('✗ POOR: End-effector unstable (RMS > 2.0 m) - consider tuning gains\n');
end

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPPORTING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_e_n, v_e_n, R_e_n] = motionCompensatedPlatform(eta_vessel, nu_vessel, ...
    q_manip, qdot_manip, r_0b)
%MOTIONCOMPENSATEDPLATFORM Forward kinematics for RRP manipulator
%   Computes end-effector position, velocity, and orientation in NED frame
%   
%   Platform configuration:
%   - Two rotational gimbals (roll q1, pitch q2) at same point
%   - One prismatic joint (extension d3) normal to deck
%   - Mounted at stern of vessel
%
%   Inputs:
%       eta_vessel - Vessel position/orientation [6x1] in NED
%       nu_vessel  - Vessel velocities [6x1] in body frame
%       q_manip    - Joint positions [q1, q2, d3]
%       qdot_manip - Joint velocities
%       r_0b       - Platform base location in body frame
%
%   Outputs:
%       p_e_n - End-effector position in NED frame
%       v_e_n - End-effector velocity in NED frame
%       R_e_n - End-effector orientation (NED to end-effector)

    % Extract vessel motion
    x_n = eta_vessel(1:3);      % Vessel position in NED
    phi = eta_vessel(4);        % Roll
    theta = eta_vessel(5);      % Pitch
    psi = eta_vessel(6);        % Yaw
    
    omega_b = nu_vessel(4:6);   % Angular velocity in body frame
    v_b = nu_vessel(1:3);       % Linear velocity in body frame
    
    % Rotation matrix from body to NED (Z-Y-X Euler angles)
    R_nb = Rzyx(phi, theta, psi);
    
    % Extract manipulator configuration
    q1 = q_manip(1);    % Roll gimbal angle [rad]
    q2 = q_manip(2);    % Pitch gimbal angle [rad]
    d3 = q_manip(3);    % Prismatic extension [m]
    
    % === Forward kinematics: R(x,q1) -> R(y,q2) -> P(z,d3) ===
    
    % Rotation about x-axis (roll compensation)
    R_x_q1 = [1   0         0;
              0   cos(q1)  -sin(q1);
              0   sin(q1)   cos(q1)];
    
    % Rotation about y-axis (pitch compensation)
    R_y_q2 = [cos(q2)   0   sin(q2);
              0         1   0;
             -sin(q2)   0   cos(q2)];
    
    % Combined rotation from platform base to end-effector
    R_e0 = R_x_q1 * R_y_q2;
    
    % End-effector position relative to platform base
    % Extension along NEGATIVE z-axis (upward, opposite to SNAME down)
    p_e0 = R_e0 * [0; 0; -d3];
    
    % Position in vessel body frame
    p_eb = r_0b + p_e0;
    
    % Position in NED frame
    p_e_n = x_n + R_nb * p_eb;
    
    % End-effector orientation in NED frame
    R_e_n = R_nb * R_e0;
    
    % === Velocity kinematics ===
    
    % Geometric Jacobian
    J_manip = manipulatorJacobian(q1, q2, d3);
    
    % Manipulator contribution to velocity
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
%MANIPULATORJACOBIAN Geometric Jacobian for RRP manipulator
%   Returns 3x3 matrix relating joint velocities to end-effector velocity
%   v = J * [q1dot; q2dot; d3dot]
%
%   Updated for upward extension: p_e0 = R_e0 * [0; 0; -d3]

    c1 = cos(q1); s1 = sin(q1);
    c2 = cos(q2); s2 = sin(q2);

    % ∂p/∂q1: Velocity due to rotation about x-axis
    dp_dq1 = [0;
              -d3*c1*c2;
              d3*s1*c2];

    % ∂p/∂q2: Velocity due to rotation about y-axis
    dp_dq2 = [-d3*c2;
              d3*s1*s2;
              d3*c1*s2];

    % ∂p/∂d3: Velocity due to prismatic extension (upward)
    dp_dd3 = [-s2;
              -s1*c2;
              -c1*c2];

    J = [dp_dq1, dp_dq2, dp_dd3];
end

function R = Rzyx(phi, theta, psi)
%RZYX Rotation matrix from body to NED using Z-Y-X Euler angles
%   Standard marine convention (yaw-pitch-roll sequence)

    cphi = cos(phi); sphi = sin(phi);
    cth = cos(theta); sth = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);
    
    R = [cth*cpsi, -cphi*spsi + sphi*sth*cpsi,  sphi*spsi + cphi*sth*cpsi;
         cth*spsi,  cphi*cpsi + sphi*sth*spsi, -sphi*cpsi + cphi*sth*spsi;
         -sth,      sphi*cth,                   cphi*cth];
end

function qdot = inverseVelocityKinematics(eta_vessel, nu_vessel, q_manip, ...
    v_desired_n, r_0b, K_rp_ff, K_heave_ff)
%INVERSEVELOCITYKINEMATICS Compute joint velocities for desired end-effector velocity
%   Uses damped least-squares inverse kinematics with motion compensation
%
%   Inputs:
%       eta_vessel  - Vessel position/orientation
%       nu_vessel   - Vessel velocities
%       q_manip     - Current joint positions
%       v_desired_n - Desired end-effector velocity in NED
%       r_0b        - Platform base location in body frame
%       K_rp_ff     - Roll/pitch feedforward gain
%       K_heave_ff  - Heave feedforward gain
%
%   Output:
%       qdot - Commanded joint velocities [3x1]

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
    
    % Manipulator Jacobian in NED frame
    J_manip_0 = manipulatorJacobian(q1, q2, d3);
    J_manip_b = R_e0 * J_manip_0;
    J_manip_n = R_nb * J_manip_b;
    
    % Vessel velocity contribution in NED
    v_vessel_n = R_nb * (v_b + cross(omega_b, p_eb));
    
    % Motion compensation feedforward
    v_comp_n = -[K_rp_ff * v_vessel_n(1:2); K_heave_ff * v_vessel_n(3)];
    
    % Required relative velocity
    v_rel_n = v_desired_n + v_comp_n - v_vessel_n;
    
    % Damped least-squares inverse (prevents singularities)
    lambda = 0.01;  % Damping factor
    J_dls = J_manip_n' / (J_manip_n * J_manip_n' + lambda^2 * eye(3));
    
    qdot = J_dls * v_rel_n;
end