%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functino computes the motion of a supply vessel with a closed-loop 
% DP controller. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
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
ship_dynamics = @(eta, nu, tau) [J6(eta) * nu; 
                                  -Minv*(D*nu + G*eta) + Minv*tau]; 

% DP controller


%% Sea state and wave spectrum 
% Sea state
Hs      = 10;               % Significant wave height [m]
gamma   = 3.3; 
beta    = deg2rad(145);     % Wave direction relative to bow [rad]

Tz = 10;            % Zero-crossing period [s]
T0 = Tz / 0.710;    % Wave spectrum modal (peak) period [s] (Fossen 2021, Eq. 10.61)
w0 = 2*pi / T0;     % Wave spectrum modal (peak) frequency [rad/s]

spectrumParam = [Hs, w0, gamma];

maxFreq = 3.0;                  % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 60;          % Number of wave frequency intervals (>50)
numDirections = 24;             % Number of wave directions (>15)

spreadingFlag = false;

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

%% Ship motion simulation
dt = 0.1;

T_final = 200;
T_initTransient = 20; 
t = (0:dt:T_final+T_initTransient-1);
N = numel(t);

% Simualtion state log
x = zeros(12,1);
x_log = zeros(N, 12);
Eta_log = zeros(N, 6);
Tau_log = zeros(N,6);

zeta_w = zeros(N_wave, 1);
TAUw_log = zeros(N_wave, 6);


