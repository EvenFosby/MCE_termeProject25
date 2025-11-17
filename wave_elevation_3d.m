%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Wave Elevation Visualization
% This script generates 3D surface plots of sea state realizations for both
% long-crested and short-crested waves using parameters from ship_motion_fRAO_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

rng(1);  % For reproducibility

%% Wave parameters from ship_motion_fRAO_v2.m
Hs = 2.5;           % Significant wave height [m]
Tz = 6;             % Zero-crossing period [s]
T0 = Tz / 0.710;    % Peak period [s]
w0 = 2*pi / T0;     % Peak frequency [rad/s]
gamma = 3.3;        % JONSWAP peak enhancement factor
beta = deg2rad(140); % Mean wave direction [rad]

% Visualization parameters
n_zscale = 5;       % Z-axis scaling factor (n times Hs for better visibility)

% Simulation parameters
g = 9.81;           % Gravity [m/s^2]
h = Inf;            % Water depth (infinite = deep water)
numFreqIntervals = 200;  % Number of wave components (200 for better visualization)
numDirections = 24;      % Number of directional bins

% Frequency range
maxFreq = 3.0;      % Maximum frequency [rad/s]
w_min = 0.2;        % Minimum frequency [rad/s]
w = linspace(w_min, maxFreq, numFreqIntervals);
dw = w(2) - w(1);

% Directional range (±90° around mean direction)
theta = linspace(beta - pi/2, beta + pi/2, numDirections);
dth = theta(2) - theta(1);

%% JONSWAP spectrum calculation
% Using the JONSWAP formulation
S = zeros(size(w));
for i = 1:length(w)
    if w(i) <= w0
        sigma = 0.07;
    else
        sigma = 0.09;
    end

    alpha = 5/16 * Hs^2 * w0^4;  % Approximate scaling
    S(i) = alpha * w(i)^(-5) * exp(-5/4 * (w(i)/w0)^(-4)) * ...
           gamma^(exp(-(w(i) - w0)^2 / (2*sigma^2*w0^2)));
end

%% Directional spreading function
% Using cos^(2s) spreading function
s_theta = 5;  % Spreading parameter (higher = more narrow)
D = cos((theta - beta)/2).^(2*s_theta);
D(abs(theta - beta) > pi/2) = 0;
D = D / (sum(D) * dth);  % Normalize

%% Spatial grid
Lx = 500;  % Domain length in x [m]
Ly = 500;  % Domain length in y [m]
Nx = 400;  % Number of points in x
Ny = 400;  % Number of points in y

x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

%% Component wavenumbers and phases
k = w2k(w, h, g);  % Solve dispersion relation

% Random phases for each frequency-direction combination
phi = 2*pi * rand(numFreqIntervals, numDirections);

%% LONG-CRESTED WAVES
% Component amplitudes (no directional spreading)
A_long = sqrt(2 * S(:) * dw);  % Only frequency dependence

% Sea-surface elevation at snapshot time
t0 = 50;  % Time snapshot [s]
eta_long = zeros(Ny, Nx);

fprintf('Computing long-crested wave elevation...\n');
for i = 1:numFreqIntervals
    % All components propagate in mean direction beta
    kx = k(i) * cos(beta);
    ky = k(i) * sin(beta);
    eta_long = eta_long + A_long(i) * cos(kx*X + ky*Y - w(i)*t0 + phi(i,1));
end

%% SHORT-CRESTED WAVES
% Component amplitudes from 2D spectrum: A = sqrt(2 * S(w) * D(th) * dw * dth)
A_short = sqrt(2 * (S(:) * D(:).') * dw * dth);  % numFreqIntervals x numDirections

% Sea-surface elevation at snapshot time
eta_short = zeros(Ny, Nx);

fprintf('Computing short-crested wave elevation...\n');
for i = 1:numFreqIntervals
    for j = 1:numDirections
        kx = k(i) * cos(theta(j));
        ky = k(i) * sin(theta(j));
        eta_short = eta_short + A_short(i,j) * cos(kx*X + ky*Y - w(i)*t0 + phi(i,j));
    end
end

%% Verification
fprintf('\n=== Wave Statistics ===\n');
fprintf('Long-crested waves:\n');
fprintf('  std(eta) = %.3f m (expected ~ Hs/4 = %.3f m)\n', std(eta_long(:)), Hs/4);
fprintf('  max(eta) = %.3f m, min(eta) = %.3f m\n', max(eta_long(:)), min(eta_long(:)));

fprintf('\nShort-crested waves:\n');
fprintf('  std(eta) = %.3f m (expected ~ Hs/4 = %.3f m)\n', std(eta_short(:)), Hs/4);
fprintf('  max(eta) = %.3f m, min(eta) = %.3f m\n', max(eta_short(:)), min(eta_short(:)));

%% Plot JONSWAP Spectrum
figure('Color','w', 'Position', [100 100 800 500]);
plot(w, S, 'b-', 'LineWidth', 2);
hold on;
xline(w0, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\omega [rad/s]', 'FontSize', 12);
ylabel('S_\eta(\omega) [m^2 s]', 'FontSize', 12);
title(sprintf('JONSWAP Spectrum (H_s = %.1f m, T_z = %.1f s, \\gamma = %.1f)', ...
    Hs, Tz, gamma), 'FontSize', 14);
legend('S(\omega)', '\omega_0', 'Location', 'northeast', 'FontSize', 11);

%% Plot Directional Spreading
figure('Color','w', 'Position', [100 100 800 500]);
plot(rad2deg(theta), D, 'b-', 'LineWidth', 2);
hold on;
xline(rad2deg(beta), 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Direction [deg]', 'FontSize', 12);
ylabel('D(\theta) [1/rad]', 'FontSize', 12);
title(sprintf('Directional Spreading Function (s = %d, \\beta = %.0f°)', ...
    s_theta, rad2deg(beta)), 'FontSize', 14);
legend('D(\theta)', 'Mean direction', 'Location', 'northeast', 'FontSize', 11);

%% Plot Long-Crested Waves
figure('Color','w', 'Position', [100 100 1000 700]);
surf(x, y, eta_long, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
shading interp;
colormap jet;
colorbar;
view(-32, 30);
grid on;
xlabel('x [m]', 'FontSize', 12);
ylabel('y [m]', 'FontSize', 12);
zlabel('z [m]', 'FontSize', 12);
zlim([-n_zscale*Hs, n_zscale*Hs]);  % Scale z-axis for better visibility
title(sprintf('Sea state realization, %d wave components (Long-Crested)', ...
    numFreqIntervals), 'FontSize', 14);
set(gca, 'FontSize', 11);

%% Plot Short-Crested Waves
figure('Color','w', 'Position', [100 100 1000 700]);
surf(x, y, eta_short, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
shading interp;
colormap jet;
colorbar;
view(-32, 30);
grid on;
xlabel('x [m]', 'FontSize', 12);
ylabel('y [m]', 'FontSize', 12);
zlabel('z [m]', 'FontSize', 12);
zlim([-n_zscale*Hs, n_zscale*Hs]);  % Scale z-axis for better visibility
title(sprintf('Sea state realization, %d wave components (Short-Crested)', ...
    numFreqIntervals*numDirections), 'FontSize', 14);
set(gca, 'FontSize', 11);

%% Comparison Plot (Side by Side)
figure('Color','w', 'Position', [50 50 1600 600]);

% Long-crested
subplot(1,2,1);
surf(x, y, eta_long, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
shading interp;
colormap jet;
colorbar;
view(-32, 30);
grid on;
xlabel('x [m]', 'FontSize', 11);
ylabel('y [m]', 'FontSize', 11);
zlabel('z [m]', 'FontSize', 11);
zlim([-n_zscale*Hs, n_zscale*Hs]);  % Scale z-axis for better visibility
title(sprintf('Long-Crested (%d components)', numFreqIntervals), 'FontSize', 12);
set(gca, 'FontSize', 10);

% Short-crested
subplot(1,2,2);
surf(x, y, eta_short, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
shading interp;
colormap jet;
colorbar;
view(-32, 30);
grid on;
xlabel('x [m]', 'FontSize', 11);
ylabel('y [m]', 'FontSize', 11);
zlabel('z [m]', 'FontSize', 11);
zlim([-n_zscale*Hs, n_zscale*Hs]);  % Scale z-axis for better visibility
title(sprintf('Short-Crested (%d components)', numFreqIntervals*numDirections), 'FontSize', 12);
set(gca, 'FontSize', 10);

sgtitle(sprintf('Wave Elevation Comparison (H_s = %.1f m, T_z = %.1f s, \\beta = %.0f°)', ...
    Hs, Tz, rad2deg(beta)), 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nPlots generated successfully!\n');

%% ========================================================================
%% USING MSS TOOLBOX waveDirectionalSpectrum FUNCTION
%% ========================================================================

fprintf('\n=== Using MSS Toolbox waveDirectionalSpectrum ===\n');

% Parameters for MSS toolbox function
spectrumType = 'JONSWAP';
spectrumParam = [Hs, w0, gamma];
omegaMax = maxFreq;

% Long-crested waves (spreadingFlag = false)
spreadingFlag_long = false;
[S_M_long, Omega_long, Amp_long, ~, ~, mu_long] = waveDirectionalSpectrum(spectrumType, ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag_long, numDirections);

% Short-crested waves (spreadingFlag = true)
spreadingFlag_short = true;
[S_M_short, Omega_short, Amp_short, ~, ~, mu_short] = waveDirectionalSpectrum(spectrumType, ...
    spectrumParam, numFreqIntervals, omegaMax, spreadingFlag_short, numDirections);

fprintf('MSS toolbox spectra computed.\n');
fprintf('  Long-crested:  S_M size = [%d x %d]\n', size(S_M_long, 1), size(S_M_long, 2));
fprintf('  Short-crested: S_M size = [%d x %d]\n', size(S_M_short, 1), size(S_M_short, 2));
fprintf('  Omega: %d frequencies\n', length(Omega_long));
fprintf('  mu: %d directions\n', length(mu_short));

% Random phases for MSS wave components
phi_mss = 2*pi * rand(length(Omega_long), length(mu_short));

% Wavenumbers for MSS frequencies
k_mss = w2k(Omega_long, h, g);

%% Long-crested waves using MSS toolbox
t0_mss = 50;  % Time snapshot [s]
eta_long_mss = zeros(Ny, Nx);

fprintf('Computing long-crested wave elevation (MSS toolbox)...\n');
for i = 1:length(Omega_long)
    kx = k_mss(i) * cos(beta);
    ky = k_mss(i) * sin(beta);
    eta_long_mss = eta_long_mss + Amp_long(i,1) * cos(kx*X + ky*Y - Omega_long(i)*t0_mss + phi_mss(i,1));
end

%% Short-crested waves using MSS toolbox
eta_short_mss = zeros(Ny, Nx);

fprintf('Computing short-crested wave elevation (MSS toolbox)...\n');
for i = 1:length(Omega_long)
    for j = 1:length(mu_short)
        kx = k_mss(i) * cos(mu_short(j));
        ky = k_mss(i) * sin(mu_short(j));
        eta_short_mss = eta_short_mss + Amp_short(i,j) * cos(kx*X + ky*Y - Omega_long(i)*t0_mss + phi_mss(i,j));
    end
end

%% Verification for MSS toolbox waves
fprintf('\n=== MSS Toolbox Wave Statistics ===\n');
fprintf('Long-crested waves (MSS):\n');
fprintf('  std(eta) = %.3f m (expected ~ Hs/4 = %.3f m)\n', std(eta_long_mss(:)), Hs/4);
fprintf('  max(eta) = %.3f m, min(eta) = %.3f m\n', max(eta_long_mss(:)), min(eta_long_mss(:)));

fprintf('\nShort-crested waves (MSS):\n');
fprintf('  std(eta) = %.3f m (expected ~ Hs/4 = %.3f m)\n', std(eta_short_mss(:)), Hs/4);
fprintf('  max(eta) = %.3f m, min(eta) = %.3f m\n', max(eta_short_mss(:)), min(eta_short_mss(:)));

%% Plot MSS Toolbox Spectra
figure('Color','w', 'Position', [100 100 1400 500]);

% Long-crested spectrum
subplot(1,2,1);
plot(Omega_long, S_M_long(:,1), 'b-', 'LineWidth', 2);
hold on;
xline(w0, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\omega [rad/s]', 'FontSize', 12);
ylabel('S_\eta(\omega) [m^2 s]', 'FontSize', 12);
title('MSS Toolbox: Long-Crested Spectrum', 'FontSize', 13);
legend('S(\omega)', '\omega_0', 'Location', 'northeast', 'FontSize', 11);

% Short-crested spectrum (multiple directions)
subplot(1,2,2);
hold on;
midIdx  = max(1, floor(length(mu_short)/2));
qtrIdx  = max(1, floor(length(mu_short)/4));
endIdx  = length(mu_short);

plot(Omega_short, S_M_short(:, midIdx), 'LineWidth', 2);
plot(Omega_short, S_M_short(:, qtrIdx), 'LineWidth', 2);
plot(Omega_short, S_M_short(:, endIdx), 'LineWidth', 2);

ylo = min(S_M_short(:)); yhi = max(S_M_short(:));
xline(w0, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('\omega [rad/s]', 'FontSize', 12);
ylabel('S_\eta(\omega,\mu) [m^2 s/rad]', 'FontSize', 12);
title('MSS Toolbox: Short-Crested Spectrum', 'FontSize', 13);
legend(sprintf('\\mu = %.0f°', rad2deg(mu_short(midIdx))), ...
       sprintf('\\mu = %.0f°', rad2deg(mu_short(qtrIdx))), ...
       sprintf('\\mu = %.0f°', rad2deg(mu_short(endIdx))), ...
       '\omega_0', 'Location', 'northeast', 'FontSize', 11);

sgtitle(sprintf('MSS Toolbox Wave Spectra (H_s = %.1f m, T_z = %.1f s)', ...
    Hs, Tz), 'FontSize', 14, 'FontWeight', 'bold');

%% Plot MSS Toolbox Wave Elevations
figure('Color','w', 'Position', [50 50 1600 600]);

% Long-crested (MSS)
subplot(1,2,1);
surf(x, y, eta_long_mss, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
shading interp;
colormap jet;
colorbar;
view(-32, 30);
grid on;
xlabel('x [m]', 'FontSize', 11);
ylabel('y [m]', 'FontSize', 11);
zlabel('z [m]', 'FontSize', 11);
zlim([-n_zscale*Hs, n_zscale*Hs]);
title(sprintf('Long-Crested (MSS: %d components)', length(Omega_long)), 'FontSize', 12);
set(gca, 'FontSize', 10);

% Short-crested (MSS)
subplot(1,2,2);
surf(x, y, eta_short_mss, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.15);
shading interp;
colormap jet;
colorbar;
view(-32, 30);
grid on;
xlabel('x [m]', 'FontSize', 11);
ylabel('y [m]', 'FontSize', 11);
zlabel('z [m]', 'FontSize', 11);
zlim([-n_zscale*Hs, n_zscale*Hs]);
title(sprintf('Short-Crested (MSS: %d components)', length(Omega_long)*length(mu_short)), 'FontSize', 12);
set(gca, 'FontSize', 10);

sgtitle(sprintf('MSS Toolbox Wave Elevations (H_s = %.1f m, T_z = %.1f s, \\beta = %.0f°)', ...
    Hs, Tz, rad2deg(beta)), 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nMSS Toolbox plots generated successfully!\n');

%% Helper function: Dispersion relation solver
function k = w2k(w, h, g)
% Solve dispersion relation: w^2 = g * k * tanh(k * h)
% For deep water (h = Inf): k = w^2 / g
    k = w.^2 / g;  % Deep-water initial guess

    if isfinite(h)
        % Newton-Raphson iteration for finite depth
        for iter = 1:30
            f = g * k .* tanh(k * h) - w.^2;
            df = g * (tanh(k * h) + k * h .* sech(k * h).^2);
            k = k - f ./ df;
        end
    end
end
