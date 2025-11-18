%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motion Compensated Platform Simulation with MRU Sensor Noise
% Loads recorded vessel motion data and adds realistic MRU measurement noise
% MRU specifications based on Kongsberg Maritime state-of-the-art sensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% Load recorded vessel motion data
fprintf('Loading recorded vessel motion data...\n');
load('vessel_motion_data.mat');

% Extract data from structure
t = motion_data.time;
eta_true = motion_data.eta;      % True positions [x, y, z, phi, theta, psi]
nu_true = motion_data.nu;        % True velocities [u, v, w, p, q, r]
N = length(t);
h = t(2) - t(1);                   % Sample time

fprintf('  Data loaded: %d samples, dt = %.4f s, duration = %.1f s\n', N, h, t(end));
fprintf('  %s\n', motion_data.description);

%% MRU Noise Specifications (Kongsberg Maritime state-of-the-art)
% Based on Kongsberg Seapath 380+ and MRU 5+ specifications

% Position noise standard deviations (RMS)
sigma_eta = [
    0.10;      % Surge (x) - 10 cm (position not directly measured, from integration)
    0.10;      % Sway (y) - 10 cm (position not directly measured, from integration)
    0.05;      % Heave (z) - 5 cm (primary measurement, very accurate)
    deg2rad(0.03);  % Roll (phi) - 0.03 degrees RMS
    deg2rad(0.03);  % Pitch (theta) - 0.03 degrees RMS
    deg2rad(0.10);  % Yaw (psi) - 0.10 degrees RMS (with gyrocompass)
];

% Velocity noise standard deviations (RMS)
sigma_nu = [
    0.05;      % Surge velocity (u) - 5 cm/s
    0.05;      % Sway velocity (v) - 5 cm/s
    0.03;      % Heave velocity (w) - 3 cm/s (heave rate is well measured)
    deg2rad(0.05);  % Roll rate (p) - 0.05 deg/s
    deg2rad(0.05);  % Pitch rate (q) - 0.05 deg/s
    deg2rad(0.08);  % Yaw rate (r) - 0.08 deg/s
];

% Optional: Small sensor bias (can drift slowly over time)
% Set to zero for no bias, or small values for realistic drift
bias_eta = zeros(6, 1);
bias_nu = zeros(6, 1);

% Optional: Bias drift rates (random walk)
drift_rate_eta = sigma_eta * 0.001;  % Very slow drift
drift_rate_nu = sigma_nu * 0.001;    % Very slow drift

fprintf('\nMRU Noise Model (Kongsberg Maritime state-of-the-art):\n');
fprintf('Position noise (RMS):\n');
fprintf('  Surge/Sway: %.1f cm, Heave: %.1f cm\n', sigma_eta(1)*100, sigma_eta(3)*100);
fprintf('  Roll/Pitch: %.3f deg, Yaw: %.3f deg\n', rad2deg(sigma_eta(4)), rad2deg(sigma_eta(6)));
fprintf('Velocity noise (RMS):\n');
fprintf('  Linear velocities: %.1f cm/s\n', sigma_nu(1)*100);
fprintf('  Angular velocities: %.3f deg/s\n', rad2deg(sigma_nu(4)));

%% Generate MRU measurements with noise
fprintf('\nGenerating MRU measurements with sensor noise...\n');

% Preallocate measurement arrays
eta_mru = zeros(N, 6);
nu_mru = zeros(N, 6);

% Preallocate bias arrays
bias_eta_log = zeros(N, 6);
bias_nu_log = zeros(N, 6);

% Set random seed for reproducibility
rng(42);

% Generate measurements for each time step
for k = 1:N
    % Update bias with random walk (slow drift)
    if k > 1
        bias_eta = bias_eta + sqrt(h) * drift_rate_eta .* randn(6, 1);
        bias_nu = bias_nu + sqrt(h) * drift_rate_nu .* randn(6, 1);
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

%% Visualization
fprintf('\nGenerating plots...\n');

% Create time vector relative to start
tt = t - t(1);

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
    sgtitle('Heave Motion - Critical for Motion Compensation');
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
    sgtitle('Position Error Distributions (should be Gaussian)');
end

fprintf('Done!\n');

%% Save MRU data for use in motion compensation simulations
mru_data.time = t;
mru_data.eta_true = eta_true;
mru_data.nu_true = nu_true;
mru_data.eta_mru = eta_mru;
mru_data.nu_mru = nu_mru;
mru_data.eta_error = eta_error;
mru_data.nu_error = nu_error;
mru_data.sigma_eta = sigma_eta;
mru_data.sigma_nu = sigma_nu;
mru_data.description = 'MRU measurements with Kongsberg Maritime state-of-the-art noise model';

save('vessel_mru_data.mat', 'mru_data');
fprintf('\nMRU data saved to vessel_motion_data.mat\n');
