clear; clc; close all;

load('vessel_motion_data.mat');

% Simulation parameters
t = motion_data.time;
dt = t(2) - t(1);
N = length(t);

% Extract motion data
% eta = [x, y, z, phi, theta, psi]' (position and Euler angles in NED)
% nu = [u, v, w, p, q, r]' (linear and angular velocities in body frame)
eta = motion_data.eta;
nu = motion_data.nu;

% Check dimensions and transpose if necessary (expect 6 x N)
if size(eta, 1) > size(eta, 2)
    eta = eta';
    fprintf('Transposed eta to 6 x %d\n', size(eta, 2));
end
if size(nu, 1) > size(nu, 2)
    nu = nu';
    fprintf('Transposed nu to 6 x %d\n', size(nu, 2));
end

%% Robot base position in ship frame (fixed mounting)
% Robot base is mounted on the deck:
%   - 10 m forward from ship origin (x_s = 10)
%   - 2 m to starboard (y_s = 2)
%   - 5 m above ship origin (z_s = -5, negative because z_s points down)
p_b_s = [-30; 1; 0];

%% Compute frame transformations for entire time series

% Preallocate arrays
R_sn = zeros(3, 3, N);  % Ship to NED rotation matrices
R_bn = zeros(3, 3, N);  % Robot base to NED rotation matrices
p_s_n = zeros(3, N);    % Ship position in NED
p_b_n = zeros(3, N);    % Robot base position in NED

% Rotation matrix from ship to robot base (constant)
R_bs = diag([-1, 1, -1]);  % 180° rotation about y-axis

fprintf('Computing frame transformations for %d time steps...\n', N);

for i = 1:N
    % Extract ship position and orientation
    x = eta(1, i);
    y = eta(2, i);
    z = eta(3, i);
    phi = eta(4, i);      % Roll
    theta = eta(5, i);    % Pitch
    psi = eta(6, i);      % Yaw (heading)

    % Ship position in NED
    p_s_n(:, i) = [x; y; z];

    % Rotation matrix from NED to ship frame (ZYX Euler angles)
    R_sn(:,:,i) = Rzyx(phi, theta, psi);

    % Rotation matrix from NED to robot base frame
    R_bn(:,:,i) = R_bs * R_sn(:,:,i);

    % Robot base position in NED
    R_ns = R_sn(:,:,i).';  % Transpose = inverse for rotation matrix
    p_b_n(:, i) = p_s_n(:, i) + R_ns * p_b_s;
end

fprintf('Frame transformations complete.\n\n');

%% Extract Euler angles for plotting
phi_deg = rad2deg(eta(4, :));
theta_deg = rad2deg(eta(5, :));
psi_deg = rad2deg(eta(6, :));

%% Plotting

%% Figure 1: Ship and Robot Base Trajectories in NED Frame (3D)
figure('Name', 'Trajectories in NED Frame', 'Position', [100 100 1200 800]);

subplot(2, 2, 1);
plot3(p_s_n(1,:), p_s_n(2,:), p_s_n(3,:), 'b-', 'LineWidth', 2);
hold on;
plot3(p_b_n(1,:), p_b_n(2,:), p_b_n(3,:), 'r-', 'LineWidth', 2);
grid on;
xlabel('North [m]');
ylabel('East [m]');
zlabel('Down [m]');
title('3D Trajectories in NED Frame');
legend('Ship Origin', 'Robot Base', 'Location', 'best');
axis equal;
set(gca, 'ZDir', 'reverse');

% Add start and end markers
plot3(p_s_n(1,1), p_s_n(2,1), p_s_n(3,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot3(p_s_n(1,end), p_s_n(2,end), p_s_n(3,end), 'gs', 'MarkerSize', 10, 'LineWidth', 2);
plot3(p_b_n(1,1), p_b_n(2,1), p_b_n(3,1), 'mo', 'MarkerSize', 10, 'LineWidth', 2);

subplot(2, 2, 2);
plot(p_s_n(1,:), p_s_n(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(p_b_n(1,:), p_b_n(2,:), 'r-', 'LineWidth', 2);
plot(p_s_n(1,1), p_s_n(2,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(p_s_n(1,end), p_s_n(2,end), 'gs', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('North [m]');
ylabel('East [m]');
title('Top View (North-East Plane)');
legend('Ship Origin', 'Robot Base', 'Start', 'End', 'Location', 'best');
axis equal;

subplot(2, 2, 3);
plot(t, p_s_n(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(t, p_b_n(3,:), 'r-', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Down [m]');
title('Vertical Position (Heave)');
legend('Ship Origin', 'Robot Base', 'Location', 'best');

subplot(2, 2, 4);
plot(t, vecnorm(p_s_n(1:2,:) - p_s_n(1:2,1)), 'b-', 'LineWidth', 2);
hold on;
plot(t, vecnorm(p_b_n(1:2,:) - p_b_n(1:2,1)), 'r-', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Horizontal Distance [m]');
title('Horizontal Displacement from Start');
legend('Ship Origin', 'Robot Base', 'Location', 'best');

%% Figure 2: Ship Orientation (Euler Angles)
figure('Name', 'Ship Orientation', 'Position', [150 150 1200 600]);

subplot(3, 1, 1);
plot(t, phi_deg, 'LineWidth', 2);
grid on;
ylabel('Roll \phi [deg]');
title('Ship Orientation (Euler Angles)');

subplot(3, 1, 2);
plot(t, theta_deg, 'LineWidth', 2);
grid on;
ylabel('Pitch \theta [deg]');

subplot(3, 1, 3);
plot(t, psi_deg, 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Yaw \psi [deg]');

%% Figure 3: Coordinate Frames Visualization at Selected Time Instants
figure('Name', 'Coordinate Frames', 'Position', [200 200 1400 800]);

% Select time instants to visualize (e.g., start, 25%, 50%, 75%, end)
time_indices = round(linspace(1, N, 5));
frame_scale = 5;  % Length of frame axes [m]

for idx = 1:length(time_indices)
    subplot(2, 3, idx);
    i = time_indices(idx);

    % Plot trajectories up to current time
    plot3(p_s_n(1,1:i), p_s_n(2,1:i), p_s_n(3,1:i), 'b-', 'LineWidth', 1);
    hold on;
    plot3(p_b_n(1,1:i), p_b_n(2,1:i), p_b_n(3,1:i), 'r-', 'LineWidth', 1);

    % Ship frame axes
    ship_origin = p_s_n(:, i);
    R_ns = R_sn(:,:,i).';
    x_s_axis = R_ns * [frame_scale; 0; 0];
    y_s_axis = R_ns * [0; frame_scale; 0];
    z_s_axis = R_ns * [0; 0; frame_scale];

    quiver3(ship_origin(1), ship_origin(2), ship_origin(3), ...
            x_s_axis(1), x_s_axis(2), x_s_axis(3), 0, 'b-', 'LineWidth', 2);
    quiver3(ship_origin(1), ship_origin(2), ship_origin(3), ...
            y_s_axis(1), y_s_axis(2), y_s_axis(3), 0, 'g-', 'LineWidth', 2);
    quiver3(ship_origin(1), ship_origin(2), ship_origin(3), ...
            z_s_axis(1), z_s_axis(2), z_s_axis(3), 0, 'r-', 'LineWidth', 2);

    % Robot base frame axes
    robot_origin = p_b_n(:, i);
    R_nb = R_bn(:,:,i).';
    x_b_axis = R_nb * [frame_scale; 0; 0];
    y_b_axis = R_nb * [0; frame_scale; 0];
    z_b_axis = R_nb * [0; 0; frame_scale];

    quiver3(robot_origin(1), robot_origin(2), robot_origin(3), ...
            x_b_axis(1), x_b_axis(2), x_b_axis(3), 0, 'b--', 'LineWidth', 1.5);
    quiver3(robot_origin(1), robot_origin(2), robot_origin(3), ...
            y_b_axis(1), y_b_axis(2), y_b_axis(3), 0, 'g--', 'LineWidth', 1.5);
    quiver3(robot_origin(1), robot_origin(2), robot_origin(3), ...
            z_b_axis(1), z_b_axis(2), z_b_axis(3), 0, 'r--', 'LineWidth', 1.5);

    grid on;
    xlabel('North [m]');
    ylabel('East [m]');
    zlabel('Down [m]');
    title(sprintf('t = %.2f s', t(i)));
    axis equal;
    %set(gca, 'ZDir', 'reverse');
    view(45, 20);

    if idx == 1
        legend('Ship Traj', 'Robot Traj', ...
               'x_s (fwd)', 'y_s (stbd)', 'z_s (down)', ...
               'x_b (back)', 'y_b (stbd)', 'y_b (up)', ...
               'Location', 'best', 'FontSize', 7);
    end
end

%% Figure 4: Relative Position of Robot Base in Ship Frame
figure('Name', 'Robot Base in Ship Frame', 'Position', [250 250 1200 400]);

% Transform robot base position back to ship frame at each time
% (should be constant since it's fixed to the ship)
p_b_s_check = zeros(3, N);
for i = 1:N
    p_b_s_check(:, i) = R_sn(:,:,i) * (p_b_n(:, i) - p_s_n(:, i));
end

% Debug: Check orthogonality of R_sn at a few time points
fprintf('\n=== Debugging R_sn orthogonality ===\n');
test_indices = [1, round(N/2), N];
for idx = test_indices
    R_test = R_sn(:,:,idx);
    identity_check = R_test * R_test';
    fprintf('Time index %d: ||R*R'' - I|| = %.2e\n', idx, norm(identity_check - eye(3), 'fro'));
end

% Debug: Check specific transformation at one time point
i_test = round(N/2);
fprintf('\n=== Debugging transformation at t = %.2f s ===\n', t(i_test));
fprintf('p_b_s (input):           [%.4f, %.4f, %.4f]\n', p_b_s);
fprintf('p_b_s_check (recovered): [%.4f, %.4f, %.4f]\n', p_b_s_check(:, i_test));
fprintf('Difference:              [%.2e, %.2e, %.2e]\n', p_b_s_check(:, i_test) - p_b_s);

subplot(1, 3, 1);
plot(t, p_b_s_check(1,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('x_s [m]');
title('Robot Base Position in Ship Frame (should be constant)');

subplot(1, 3, 2);
plot(t, p_b_s_check(2,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('y_s [m]');
title('Robot Base Position in Ship Frame');

subplot(1, 3, 3);
plot(t, p_b_s_check(3,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('z_s [m]');
title('Robot Base Position in Ship Frame');

%% Display statistics
fprintf('=== Simulation Statistics ===\n');
fprintf('Total time: %.2f s\n', t(end));
fprintf('Time step: %.4f s\n', dt);
fprintf('Number of samples: %d\n', N);
fprintf('\n--- Ship Motion Range ---\n');
fprintf('Position (NED):\n');
fprintf('  North: [%.2f, %.2f] m\n', min(p_s_n(1,:)), max(p_s_n(1,:)));
fprintf('  East:  [%.2f, %.2f] m\n', min(p_s_n(2,:)), max(p_s_n(2,:)));
fprintf('  Down:  [%.2f, %.2f] m\n', min(p_s_n(3,:)), max(p_s_n(3,:)));
fprintf('Orientation:\n');
fprintf('  Roll:  [%.2f, %.2f] deg\n', min(phi_deg), max(phi_deg));
fprintf('  Pitch: [%.2f, %.2f] deg\n', min(theta_deg), max(theta_deg));
fprintf('  Yaw:   [%.2f, %.2f] deg\n', min(psi_deg), max(psi_deg));
fprintf('\n--- Robot Base Motion Range (in NED) ---\n');
fprintf('Position (NED):\n');
fprintf('  North: [%.2f, %.2f] m\n', min(p_b_n(1,:)), max(p_b_n(1,:)));
fprintf('  East:  [%.2f, %.2f] m\n', min(p_b_n(2,:)), max(p_b_n(2,:)));
fprintf('  Down:  [%.2f, %.2f] m\n', min(p_b_n(3,:)), max(p_b_n(3,:)));
fprintf('\n--- Robot Base Position in Ship Frame (verification) ---\n');
fprintf('Expected: [%.2f, %.2f, %.2f] m\n', p_b_s);
fprintf('Computed mean: [%.4f, %.4f, %.4f] m\n', mean(p_b_s_check, 2));
fprintf('Standard deviation: [%.2e, %.2e, %.2e] m\n', std(p_b_s_check, 0, 2));

%% Helper Functions

function R = Rzyx(phi, theta, psi)
    % RZYX Create rotation matrix from ZYX Euler angles (roll-pitch-yaw)
    % Rotation from NED frame to body frame
    % Input: phi (roll), theta (pitch), psi (yaw) in radians
    % Output: R = Rz(psi) * Ry(theta) * Rx(phi)

    cpsi = cos(psi);   spsi = sin(psi);
    cth = cos(theta);  sth = sin(theta);
    cphi = cos(phi);   sphi = sin(phi);

    R = [cpsi*cth, -spsi*cphi + cpsi*sth*sphi,  spsi*sphi + cpsi*sth*cphi;
         spsi*cth,  cpsi*cphi + spsi*sth*sphi, -cpsi*sphi + spsi*sth*cphi;
         -sth,      cth*sphi,                   cth*cphi];
end
