clear; close all; clc;
rng(1);
load supply;

% Simulation parameters
h = 0.1;
T_final = 200;
T_initTransient = 50;
t = 0:h:T_final+T_initTransient-1;
N = numel(t);

% Sea state
Hs = 2.5;
Tz = 10;
beta = deg2rad(140);
T0 = Tz / 0.710;
w0 = 2*pi/T0;
gamma = 3.3;
spectrumParameters = [Hs, w0, gamma];

maxFreq = 3.0;              % Maximum frequency in RAO computations (rad/s) 
numFreqIntervals = 100;     % Number of wave frequency intervals (>50)
numDirections = 24;         % Number of wave directions (>15)

spectrumType = 'JONSWAP';
spreadingFlag = 1;

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
    spectrumParameters, numFreqIntervals, omegaMax, spreadingFlag, numDirections);

%% MAIN LOOP
t = 0:h:T_final+T_initTransient-1;  % Time vector
simdata = zeros(length(t),7);        % Pre-allocate table

for i = 1:length(t)
    U = 5;                           % Time-varying ship speed (m/s)
    psi = deg2rad(sin(0.1 * t(i)));  % Time-varying heading angle (rad)
    
    % 6-DOF generalized wave forces
    [tau_wave1, waveElevation] = waveForceRAO(t(i), ...
        S_M, Amp, Omega, mu, vessel, U, psi, beta, numFreqIntervals);
    
    simdata(i,:) = [tau_wave1' waveElevation];
end

% Remove initial transient
idx_start = find(t >= T_initTransient, 1);
t_plot = t(idx_start:end) - T_initTransient;
simdata_plot = simdata(idx_start:end, :);

%% POST-PROCESSING AND PLOTTING

% DOF labels
dof_labels = {'Surge (X)', 'Sway (Y)', 'Heave (Z)', 'Roll (K)', 'Pitch (M)', 'Yaw (N)'};
dof_units = {'N', 'N', 'N', 'Nm', 'Nm', 'Nm'};

%% 1. WAVE SPECTRUM PLOT
figure('Name', 'Wave Spectrum Analysis', 'Position', [100 100 1200 800]);

% Plot 1: Wave spectrum
subplot(2,2,1)
if size(S_M, 2) == 1
    % Omnidirectional spectrum
    plot(Omega, S_M, 'b-', 'LineWidth', 2);
    xlabel('Frequency \omega (rad/s)', 'FontSize', 11);
    ylabel('Spectral Density S(\omega) (m^2s/rad)', 'FontSize', 11);
    title('Wave Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
else
    % Directional spectrum - plot sum over all directions
    S_total = sum(S_M, 2);
    plot(Omega, S_total, 'b-', 'LineWidth', 2);
    xlabel('Frequency \omega (rad/s)', 'FontSize', 11);
    ylabel('Spectral Density S(\omega) (m^2s/rad)', 'FontSize', 11);
    title('Wave Spectrum (Integrated over all directions)', 'FontSize', 12, 'FontWeight', 'bold');
end
grid on;
xlim([0 max(Omega)]);

% Plot 2: Significant wave height comparison
subplot(2,2,2)
% Calculate Hs from simulated wave elevation
eta = simdata_plot(:, 7);
Hs_simulated = 4 * std(eta);

bar_data = [Hs, Hs_simulated];
bar(bar_data);
set(gca, 'XTickLabel', {'Input H_s', 'Simulated H_s'});
ylabel('Significant Wave Height (m)', 'FontSize', 11);
title('H_s Comparison', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([0 max(bar_data)*1.2]);

% Add text labels on bars
for i = 1:length(bar_data)
    text(i, bar_data(i) + 0.1, sprintf('%.3f m', bar_data(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

% Calculate error
error_percent = abs(Hs_simulated - Hs) / Hs * 100;
text(1.5, max(bar_data)*1.1, sprintf('Error: %.2f%%', error_percent), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);

% Plot 3: Wave elevation time series
subplot(2,2,3)
plot(t_plot, eta, 'b-', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Wave Elevation \eta (m)', 'FontSize', 11);
title('Wave Elevation Time History', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0 min(100, T_final)]);

% Plot 4: Wave elevation statistics
subplot(2,2,4)
histogram(eta, 50, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none');
hold on;
% Overlay theoretical Gaussian distribution (manual calculation)
eta_mean = mean(eta);
eta_std = std(eta);
eta_range = linspace(min(eta), max(eta), 200);
% Manual Gaussian PDF: f(x) = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((x-mu)/sigma)^2)
pdf_theoretical = (1/(eta_std*sqrt(2*pi))) * exp(-0.5*((eta_range-eta_mean)/eta_std).^2);
plot(eta_range, pdf_theoretical, 'r-', 'LineWidth', 2);
xlabel('Wave Elevation \eta (m)', 'FontSize', 11);
ylabel('Probability Density', 'FontSize', 11);
title('Wave Elevation Distribution', 'FontSize', 12, 'FontWeight', 'bold');
legend('Simulated', 'Theoretical Gaussian', 'Location', 'best');
grid on;

%% 2. WAVE FORCES TIME SERIES
figure('Name', 'Wave Forces on Vessel', 'Position', [150 50 1400 900]);

for dof = 1:6
    subplot(3, 2, dof)
    plot(t_plot, simdata_plot(:, dof), 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', 10);
    ylabel([dof_labels{dof} ' (' dof_units{dof} ')'], 'FontSize', 10);
    title(['Wave Force: ' dof_labels{dof}], 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    xlim([0 min(100, T_final)]);
    
    % Add statistics text
    force_rms = rms(simdata_plot(:, dof));
    force_max = max(abs(simdata_plot(:, dof)));
    text(0.98, 0.95, sprintf('RMS: %.2f\nMax: %.2f', force_rms, force_max), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'BackgroundColor', 'w');
end

%% 3. FORCE SPECTRUM ANALYSIS
figure('Name', 'Force Spectrum Analysis', 'Position', [200 100 1400 900]);

for dof = 1:6
    subplot(3, 2, dof)
    
    % Calculate power spectral density using Welch's method
    [psd, freq] = pwelch(simdata_plot(:, dof), hamming(512), 256, 512, 1/h);
    
    % Convert frequency from Hz to rad/s
    omega_force = 2*pi*freq;
    
    plot(omega_force, psd, 'LineWidth', 1.5);
    xlabel('Frequency \omega (rad/s)', 'FontSize', 10);
    ylabel(['PSD (' dof_units{dof} ')^2/(rad/s)'], 'FontSize', 10);
    title(['Force Spectrum: ' dof_labels{dof}], 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    xlim([0 2]);
    
    % Find peak frequency
    [peak_val, peak_idx] = max(psd);
    peak_omega = omega_force(peak_idx);
    hold on;
    plot(peak_omega, peak_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    text(peak_omega, peak_val, sprintf('  Peak: %.3f rad/s', peak_omega), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom');
end

%% 4. COMPREHENSIVE SUMMARY FIGURE
figure('Name', 'Summary Statistics', 'Position', [250 150 1200 600]);

% Plot 1: RMS values of all forces
subplot(1,3,1)
force_rms = rms(simdata_plot(:, 1:6), 1);
bar(force_rms);
set(gca, 'XTickLabel', {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'});
ylabel('RMS Force/Moment', 'FontSize', 11);
title('RMS Wave Forces', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xtickangle(45);

% Add value labels
for i = 1:6
    text(i, force_rms(i), sprintf('%.1f', force_rms(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9);
end

% Plot 2: Maximum values of all forces
subplot(1,3,2)
force_max = max(abs(simdata_plot(:, 1:6)), [], 1);
bar(force_max);
set(gca, 'XTickLabel', {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'});
ylabel('Max Force/Moment', 'FontSize', 11);
title('Maximum Wave Forces', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xtickangle(45);

% Add value labels
for i = 1:6
    text(i, force_max(i), sprintf('%.1f', force_max(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9);
end

% Plot 3: Summary statistics table
subplot(1,3,3)
axis off;
summary_text = sprintf('SIMULATION SUMMARY\n\n');
summary_text = [summary_text sprintf('Sea State Parameters:\n')];
summary_text = [summary_text sprintf('  H_s (input): %.2f m\n', Hs)];
summary_text = [summary_text sprintf('  H_s (simulated): %.2f m\n', Hs_simulated)];
summary_text = [summary_text sprintf('  Error: %.2f%%\n', error_percent)];
summary_text = [summary_text sprintf('  T_z: %.2f s\n', Tz)];
summary_text = [summary_text sprintf('  Wave direction: %.1f°\n', rad2deg(beta))];
summary_text = [summary_text sprintf('\nSimulation Parameters:\n')];
summary_text = [summary_text sprintf('  Duration: %.1f s\n', T_final)];
summary_text = [summary_text sprintf('  Time step: %.2f s\n', h)];
summary_text = [summary_text sprintf('  Ship speed: %.1f m/s\n', U)];
summary_text = [summary_text sprintf('\nWave Elevation:\n')];
summary_text = [summary_text sprintf('  Mean: %.3f m\n', mean(eta))];
summary_text = [summary_text sprintf('  Std Dev: %.3f m\n', std(eta))];
summary_text = [summary_text sprintf('  Max: %.3f m\n', max(eta))];
summary_text = [summary_text sprintf('  Min: %.3f m\n', min(eta))];

text(0.1, 0.95, summary_text, 'FontSize', 10, 'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%% Print summary to console
fprintf('\n=================================================\n');
fprintf('WAVE FORCE SIMULATION SUMMARY\n');
fprintf('=================================================\n\n');
fprintf('Sea State:\n');
fprintf('  Input H_s:      %.3f m\n', Hs);
fprintf('  Simulated H_s:  %.3f m\n', Hs_simulated);
fprintf('  Error:          %.2f%%\n', error_percent);
fprintf('  T_z:            %.2f s\n', Tz);
fprintf('  Wave direction: %.1f°\n\n', rad2deg(beta));

fprintf('Wave Forces (RMS values):\n');
for dof = 1:6
    fprintf('  %s: %.2f %s\n', dof_labels{dof}, force_rms(dof), dof_units{dof});
end

fprintf('\nWave Forces (Maximum values):\n');
for dof = 1:6
    fprintf('  %s: %.2f %s\n', dof_labels{dof}, force_max(dof), dof_units{dof});
end
fprintf('\n=================================================\n');