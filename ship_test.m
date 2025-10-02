clear; clc; close all;

%% ---------------- Ship & hydrostatics -----------------------------------
load supply;   % provides 'vessel'

MRB = vessel.MRB;
MA  = vessel.A(:,:,1);        % added mass slice (pick what you want)
M   = MRB + MA;               % total inertia

rho = 1025; g = 9.81;
Awp  = vessel.main.Lwl * vessel.main.B * 0.8;
GM_T = vessel.main.GM_T;
GM_L = vessel.main.GM_L;
m    = vessel.main.m;

% Linear restoring
G = diag([ M(1,1)*(0.05)^2, ...
           M(2,2)*(0.05)^2, ...
           rho*g*Awp, ...
           m*g*GM_T, ...
           m*g*GM_L, ...
           M(6,6)*(0.05)^2 ]);

% Linear damping
zeta = [1; 1; 0.20; 0.03; 0.05; 1];
Dlin = zeros(6);
for i = 1:6
    if G(i,i) > 0
        Dlin(i,i) = 2*zeta(i)*sqrt(M(i,i)*G(i,i));
    else
        Dlin(i,i) = 0;
    end
end

%% ---------------- Sea state & spectrum ----------------------------------
Hs    = 10; 
gamma = 3.3;
beta  = deg2rad(145);       % wave dir relative to bow (check convention)

Tz = 10;
T0 = Tz / 0.710;
w0 = 2*pi / T0;
spectrumParam = [Hs, w0, gamma];

numFreqIntervals = 60;
numDirections    = 24;
spreadingFlag    = false;

omegaMax = vessel.forceRAO.w(end);
[S_M, omega, Amp, ~, ~, mu] = waveDirectionalSpectrum('JONSWAP', ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% ---------------- Multi-rate simulation ---------------------------------
dt_ship = 0.01;    % fine step for vessel dynamics
dt_wave = 0.10;    % coarse step for wave forcing
ratio   = dt_wave / dt_ship;
assert(abs(ratio - round(ratio)) < 1e-12, 'dt_wave must be an integer multiple of dt_ship');
ratio = round(ratio);

T_final = 200;
t_ship  = (0:dt_ship:T_final).';
t_wave  = (0:dt_wave:T_final).';

N_ship = numel(t_ship);
N_wave = numel(t_wave);

% State x = [eta(6); nu(6)]
x = zeros(12,1);
X_log   = zeros(N_ship,12);
ETA_log = zeros(N_ship,6);
TAU_log = zeros(N_ship,6);
zeta_w  = zeros(N_wave,1);        % wave elevation (coarse log)
TAUw    = zeros(N_wave,6);        % wave forces (coarse log)

% DP setpoint and simple PD (tune as needed)
eta_ref = [0;0;0;0;0;deg2rad(0)];
Kp = diag([2e5 2e5 0 5e7 5e7 5e8]);
Kd = diag([2e4 2e4 0 5e6 5e6 5e7]);

% Dynamics function for RK4: xdot = f(x; tau)
dyn = @(x_local, M_, D_, G_, tau_) [ ...
    x_local(7:12) ; ...
    M_ \ ( -D_*x_local(7:12) - G_*x_local(1:6) + tau_ ) ];

ship_idx = 1;

for k = 1:N_wave
    % --- Evaluate wave load at coarse time using current state ---
    eta = x(1:6);  nu = x(7:12);
    psi = eta(6);
    u = nu(1); v = nu(2);
    U = hypot(u, v);      % if RAO expects forward speed use u instead

    [tau_wave, waveElev] = waveForceRAO(t_wave(k), S_M, Amp, omega, mu, ...
                                        vessel, U, psi, beta, numFreqIntervals);
    TAUw(k,:) = tau_wave(:).';
    zeta_w(k) = waveElev;

    % --- Do 'ratio' fine RK4 steps with ZOH wave load ---
    for s = 1:ratio
        % Log current state on the fine grid
        X_log(ship_idx,:)   = x.';
        ETA_log(ship_idx,:) = x(1:6).';

        % DP control (update every fine step)
        eta = x(1:6);  nu = x(7:12);
        e_eta = eta - eta_ref;
        e_eta(6) = atan2(sin(e_eta(6)), cos(e_eta(6)));  % wrap yaw error
        tau_ctrl = -Kp*e_eta - Kd*nu;

        tau = tau_wave + tau_ctrl;
        TAU_log(ship_idx,:) = tau(:).';

        % Advance one fine step with your RK4
        if ship_idx < N_ship
            x = rk4(dyn, dt_ship, x, M, Dlin, G, tau);
        end
        ship_idx = ship_idx + 1;

        % Stop if we've reached the exact end on the fine grid
        if ship_idx > N_ship
            break;
        end
    end

    if ship_idx > N_ship
        break;
    end
end

y = ETA_log;  % positions/angles on the fine grid

%% ---------------- Natural frequencies (heave/roll/pitch) ----------------
wn = zeros(1,6);  Tn = nan(1,6);
for i = 1:6
    if G(i,i) > 0
        wn(i) = sqrt(G(i,i)/M(i,i));
        Tn(i) = 2*pi/wn(i);
    end
end
fprintf('Natural frequencies (rad/s): heave=%.3f, roll=%.3f, pitch=%.3f\n', wn(3), wn(4), wn(5));
fprintf('Natural periods (s):         heave=%.1f, roll=%.1f, pitch=%.1f\n', Tn(3), Tn(4), Tn(5));

%% ---------------- Plots --------------------------------------------------
dof_labels = {'x [m] (surge)','y [m] (sway)','z [m] (heave)', ...
              '\phi [rad] (roll)','\theta [rad] (pitch)','\psi [rad] (yaw)'};

figure('Color','w');
for i = 1:6
    subplot(3,2,i)
    plot(t_ship, y(:,i), 'LineWidth', 1); grid on
    xlabel('Time [s]'); ylabel(dof_labels{i});
    title(['Response of ' dof_labels{i}])
end
sgtitle('6-DOF Ship Motion (RK4 @ 0.01 s, Waves @ 0.10 s)');

figure('Color','w');
U_est = hypot([diff(X_log(:,1))./dt_ship; 0], [diff(X_log(:,2))./dt_ship; 0]);
plot(t_ship, U_est, 'LineWidth', 1); grid on
xlabel('Time [s]'); ylabel('Horizontal speed U [m/s]');
title('Estimated horizontal speed');

figure('Color','w');
stairs(t_wave, TAUw(:,1), 'LineWidth', 1); grid on
xlabel('Time [s]'); ylabel('\tau_{wave,surge} [N]');
title('Wave surge force (ZOH @ 0.1 s)');
