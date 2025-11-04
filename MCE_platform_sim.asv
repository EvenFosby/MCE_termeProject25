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
forceRaoFlag = false;
useIntegralAction = false;

% Simulation parameters
h = 0.1;
T_final = 300;
T_initTransient = 0;

t = 0:h:T_final+T_initTransient-1;
N = numel(t);

% Control objectiv
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

Kd = 2.*M*Zeta_pid*Omega_n; % - D; 

Ki = 0.10*Kp*Omega_n;

% PID controller
tau_pid = @(eta, nu, eta_int) S.*(eulerang(eta(4), eta(5), eta(6))'*(-Kp*eta ...
    - Kd*eulerang(eta(4), eta(5), eta(6))*nu - (useIntegralAction*Ki*eta_int)));

% Heading lowpass filter
T_psi = 12;
alpha = h/(T_psi + h);
psi_lp = 0;

%% Motion compensated platform configuration

% Motion compensation mode
mc_mode = 'stabilize';

% Target position {n}
p_target_n = [50; 20; -8];

% platform base location  {b}
r_0b = [30; 0; -5];

% Joints limits
q1_lim = deg2rad([-15, 15]); % Roll limits
q2_lim = deg2rad([-10, 10]); % Pitch limits
d3_lim = [3, 6];  % Cylinder extension limits

% Nominal cylinder extension
d3_nominal = 4.5; 

% Initial manipulator configuration
q_manip = [0; 0; d3_nominal];
qdot_manip = [0; 0; 0];



