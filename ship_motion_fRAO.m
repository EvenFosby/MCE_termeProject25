clear; clc; close all;
rng(1);

%% Ship & hydrostatics
load supply;

% Geometry and mass
% r_g = [0 0 0];
m = vessel.main.m; 

MRB     = vessel.MRB;
MA      = vessel.A(:,:,1);
M       = MRB + MA;

% Hydrostatics
rho  = 1025; g = 9.81;
Awp  = vessel.main.Lwl * vessel.main.B * 0.8; % Water plane displacement.
GM_T = vessel.main.GM_T; 
GM_L = vessel.main.GM_L;

% Linear restoring
G = diag([M(1,1)*0.05^2, ...
          M(2,2)*0.05^2, ...
          rho*g*Awp, ...
          m*g*GM_T, ...
          m*g*GM_L, ...
          M(6,6)*0.05^2]);

% Linear damping
zeta = [ 1; 1; 0.20; 0.03; 0.05; 1 ];         % target damping ratios
D    = zeros(6);

for i = 1:6
        % c_i = 2*zeta_i*sqrt(M_i*K_i)  (critical damping scaling)
        D(i,i) = 2*zeta(i)*sqrt(M(i,i)*G(i,i));
end

J = eulerang()
% State-space model
A = [ zeros(6) eye(6);
      -M\G   -M\D ];

B = [ zeros(6);
      M\eye(6) ];

C = [ eye(6) zeros(6) ];
D = zeros(6);

ship = @(x,tau_wave) A*x + B*tau_wave; 

%% Sea state and wave spectrum 
% Sea state
Hs      = 10;               % Significant wave height [m]
gamma   = 3.3; 
beta    = deg2rad(145);     % Wave direction relative to bow [rad]

Tz = 10;            % Zero-crossing period [s]
T0 = Tz / 0.710;    % Wave spectrum modal (peak) period [s] (Fossen 2021, Eq. 10.61)
w0 = 2*pi / T0;     % Wave spectrum modal (peak) frequency [rad/s]

spectrumParam = [Hs, w0, gamma];

numFreqIntervals = 60;          % Number of wave frequency intervals (>50)
numDirections = 24;             % Number of wave directions (>15)

spreadingFlag = false;

omegaMax = vessel.forceRAO.w(end);
[S_M, omega, Amp, ~, ~, mu] = waveDirectionalSpectrum('JONSWAP', ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% Check difference between implementations
% Wave simulation
h_wave = 0.05; % Wave simulation time step [s]
T_final = 200; % Simulation lenght [s]
T_initTransient = 20; % Remove initial transient [s]

t_wave = 0:h_wave:T_final+T_initTransient-1;
simData = zeros(length(t_wave), 7);

for i = 1:length(t_wave)
    U = 3;
    psi = deg2rad(sin(0.1 * t_wave(i)));

    [tau_wave, waveElevation] = waveForceRAO(t_wave(i), S_M, Amp, omega, mu, ...
        vessel, U, psi, beta, numFreqIntervals);

    simData(i,:) = [tau_wave' waveElevation];
end

% Plot RAO wave forces and wave elevation with const U
figure(10); clf;

DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)', ...
           'Roll (Nm)', 'Pitch (Nm)', 'Yaw (Nm)'};

t = t_wave(:);           % ensure column vector
tau_wave_all = simData(:,1:6);
zeta_w = simData(:,7);

% Plot 6-DOF wave forces
for dof = 1:6
    subplot(7,1,dof); hold on; grid on;
    plot(t, tau_wave_all(:,dof), 'LineWidth', 1.8);
    ylabel(DOF_txt{dof});
    if dof < 6
        set(gca,'XTickLabel',[]);
    end
end

% Plot wave elevation
figure(11); hold on; grid on;
plot(t, zeta_w, 'b', 'LineWidth', 1.8);
ylabel('Elevation (m)');
xlabel('Time (s)');


%% Ship motion simulation
dt_ship = 0.05; % Simulation time step, ship
dt_wave = 1; % Simulation time step, wave
ratio   = dt_wave / dt_ship;
assert(abs(ratio - round(ratio)) < 1e-12, 'dt_wave must be an integer multiple of dt_ship');
ratio = round(ratio);

T_final = 200; % Simulation lenght
T_initTransient = 20; % Remove initial transient [s]
t_ship = (0:dt_ship:T_final+T_initTransient-1).';
t_wave = (0:dt_wave:T_final+T_initTransient-1).';

N_ship = numel(t_ship);
N_wave = numel(t_wave);

% Simulation state log
x = zeros(12,1); % Ship state x=[eta(6), nu(6)]'
X_log = zeros(N_ship, 12);
ETA_log = zeros(N_ship, 6);
TAU_log = zeros(N_ship, 6);

zeta_w = zeros(N_wave, 1);
TAUw_log = zeros(N_wave, 6);

ship_idx = 1;

psi = deg2rad(0);


% Main loop
for i = 1:N_wave
    eta = x(1:6);
    nu = x(7:12);   
    psi = eta(6);   % Heading angle
    U = sqrt(nu(1)^2 + nu(2)^2);

    [tau_wave, waveElev] = waveMotionRAO(t_wave(i), S_M, Amp, omega, mu, ...
                                        vessel, 0, psi, beta, numFreqIntervals);
    TAUw_log(i,:) = tau_wave(:).';
    zeta_w(i) = waveElev;

    for k = 1:ratio
        % Log current state on the fine grid
        X_log(ship_idx,:)   = x.';
        ETA_log(ship_idx,:) = x(1:6).';
        TAU_log(ship_idx,:) = tau_wave(:).';

        % Advance one fine step with your RK4
        if ship_idx < N_ship
            x = rk4(ship, dt_ship, x, tau_wave);
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


%% -------- Report natural frequencies (heave/roll/pitch) -----------------
wn = zeros(1,6);  Tn = nan(1,6);
for i = 1:6
    if G(i,i) > 0
        wn(i) = sqrt(G(i,i)/M(i,i));   % rad/s
        Tn(i) = 2*pi/wn(i);             % s
    end
end
fprintf('Natural frequencies (rad/s): heave=%.3f, roll=%.3f, pitch=%.3f\n', wn(3), wn(4), wn(5));
fprintf('Natural periods (s):         heave=%.1f, roll=%.1f, pitch=%.1f\n', Tn(3), Tn(4), Tn(5));


%% ====== Plotting for Ship Motion Simulator ======
% Assumes these exist from your sim:
% S_M (numFreq x numDir), omega (numFreq x 1), w0, mu (1 x numDir),
% spreadingFlag (logical), Hs, beta (rad), zeta_w (N_wave x 1), t_wave,
% TAUw_log (N_wave x 6), ETA_log (N_ship x 6), t_ship, T_initTransient

% --- Helper: trim initial transient ---
idx_wave = t_wave >= T_initTransient;
idx_ship = t_ship >= T_initTransient;

t_wave_ss = t_wave(idx_wave);
t_ship_ss = t_ship(idx_ship);

zeta_w_ss = zeta_w(idx_wave);
TAUw_ss    = TAUw_log(idx_wave, :);
ETA_ss     = ETA_log(idx_ship, :);

%% 1) Wave spectrum + wave elevation (like your example)
figure(1); clf;

% Wave spectrum
subplot(2,1,1); hold on; grid on;
if spreadingFlag
    % pick a few representative directions
    id_mid  = max(1, floor(length(mu)/2));
    id_qtr  = max(1, floor(length(mu)/4));
    id_end  = length(mu);

    plot(omega, S_M(:, id_mid), 'LineWidth', 2);
    plot(omega, S_M(:, id_qtr), 'LineWidth', 2);
    plot(omega, S_M(:, id_end), 'LineWidth', 2);

    % Peak frequency marker
    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 2);

    legend( ...
        sprintf('\\mu = %.0f deg', rad2deg(mu(id_mid))), ...
        sprintf('\\mu = %.0f deg', rad2deg(mu(id_qtr))), ...
        sprintf('\\mu = %.0f deg', rad2deg(mu(id_end))), ...
        sprintf('w_0 = %.3g rad/s', w0), ...
        'Location', 'best');
else
    plot(omega, S_M(:,1), 'LineWidth', 2);
    ylo = min(S_M(:)); yhi = max(S_M(:));
    plot([w0 w0], [ylo yhi], 'k--', 'LineWidth', 2);
    legend('S(\Omega)', sprintf('w_0 = %.3g rad/s', w0), 'Location', 'best');
end
xlabel('\Omega (rad/s)');
ylabel('m^2 s');
title(sprintf('%s spectrum', 'JONSWAP')); % adjust if you vary type

% Wave elevation (steady-state)
subplot(2,1,2); hold on; grid on;
plot(t_wave_ss, zeta_w_ss, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('m');
title(sprintf('Wave Elevation for \\beta_{wave} = %.0f%c, H_s = %.2f m', ...
      rad2deg(beta), 176, Hs));

%% 2) 6-DOF 1st-order wave forces (RAO loads) on wave time grid
figure(2); clf;
DOF_txt = {'Surge (N)', 'Sway (N)', 'Heave (N)', 'Roll (Nm)', 'Pitch (Nm)', 'Yaw (Nm)'};
% Optional scaling (keep as ones unless you really want to rescale moments)
T_scale = [1 1 1 1 1 1];

for dof = 1:6
    subplot(6,1,dof); hold on; grid on;
    plot(t_wave_ss, T_scale(dof)*TAUw_ss(:,dof), 'LineWidth', 2);
    ylabel(DOF_txt{dof});
    if dof == 6
        xlabel('Time (s)');
    else
        set(gca,'XTickLabel',[]);
    end
end
sgtitle('1st-Order Wave Forces (steady-state)');

%% 3) (Optional) Ship motions \eta (positions/angles) on the ship time grid
%    Useful to visualize response; angles in rad by default.
figure(3); clf;
eta_lbl = {'\eta_1 Surge (m)','\eta_2 Sway (m)','\eta_3 Heave (m)', ...
           '\eta_4 Roll (rad)','\eta_5 Pitch (rad)','\eta_6 Yaw (rad)'};

for i = 1:6
    subplot(6,1,i); hold on; grid on;
    plot(t_ship_ss, ETA_ss(:,i), 'LineWidth', 1.8);
    ylabel(eta_lbl{i});
    if i == 6
        xlabel('Time (s)');
    else
        set(gca,'XTickLabel',[]);
    end
end
sgtitle('Ship Motions \eta (steady-state)');

% If you prefer roll/pitch/yaw in degrees, uncomment:
% ETA_deg = ETA_ss; ETA_deg(:,4:6) = rad2deg(ETA_deg(:,4:6));
% ...and plot ETA_deg instead.








% %% ==== Post-processing & plots ====
% 
% % Labels for the 6 DOFs (in your state ordering: [x y z φ θ ψ] = [surge sway heave roll pitch yaw])
% dof_labels = {'Surge [m]','Sway [m]','Heave [m]','Roll [rad]','Pitch [rad]','Yaw [rad]'};
% 
% % Remove initial transient
% idx = t_ship >= T_initTransient;
% ts  = t_ship(idx);
% Y   = ETA_log(idx,:);    % Motions (6 DOFs)
% 
% % 1) Motion responses (all 6 DOFs)
% figure('Color','w'); 
% tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
% for i = 1:6
%     nexttile
%     plot(ts, Y(:,i), 'LineWidth', 1.2); grid on
%     xlabel('Time [s]')
%     ylabel(dof_labels{i})
%     title(['Response of ' dof_labels{i}])
% end
% 
% %% 2) Wave spectrum S(omega)
% % S_M may be either Nw x Ndir (directional) or Nw x 1 (already 1D).
% % We'll integrate over directions if needed to get a 1D spectrum S(omega).
% S_omega = [];
% if ismatrix(S_M) && size(S_M,2) > 1
%     % Integrate over directions using trapezoidal rule (mu in radians)
%     % Assumes mu is monotonically increasing.
%     S_omega = trapz(mu, S_M, 2);  % -> Nw x 1
% else
%     % Already 1D spectrum
%     S_omega = S_M(:);
% end
% 
% % Convert to column vectors for plotting safety
% omega_vec = omega(:);
% S_omega   = S_omega(:);
% 
% % 3) Wave elevation time series (remove transient)
% idxw = t_wave >= T_initTransient;
% tw   = t_wave(idxw);
% zeta = zeta_w(idxw);
% 
% % Estimated significant wave height from the simulated elevation
% Hs_est = 4*sqrt(var(zeta));  % Rayleigh approx: Hs = 4*sqrt(m0) ~ 4*std(η)
% 
% % Optional: theoretical Hs line for comparison (you set Hs earlier)
% Hs_theory = Hs;
% 
% % Plot spectrum and elevation in one figure
% figure('Color','w'); 
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% 
% % (a) Spectrum S(omega)
% nexttile
% plot(omega_vec, S_omega, 'LineWidth', 1.5); grid on
% xlabel('\omega [rad/s]'); ylabel('S(\omega) [m^2 s]')
% title('Wave Spectrum')
% if exist('w0','var')
%     hold on
%     yL = ylim;
%     plot([w0 w0], yL, '--', 'LineWidth', 1); % mark peak frequency
%     legend('S(\omega)','\omega_0 (peak)','Location','best')
%     hold off
% end
% 
% % (b) Wave elevation
% nexttile
% plot(tw, zeta, 'LineWidth', 1.2); grid on
% xlabel('Time [s]'); ylabel('\eta [m]')
% title(sprintf('Wave Elevation (H_s^{est} = %.2f m, H_s = %.2f m)', Hs_est, Hs_theory))
