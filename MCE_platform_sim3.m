clear; clc; close all;

load('vessel_mru_data.mat');

% Simulation parameters
t = mru_data.time;
dt = t(2) - t(1);
N = length(t);

% Extract motion data
eta = mru_data.eta_mru'; % Position and Euler angles in NED
nu = mru_data.nu_mru'; % Linear and angular velocities in body frame

%% Robot base position in ship frame
p_b_s = [-30; 0; -3];

%% Compute frame transformations using homogeneous transforms
% Preallocate arrays for homogeneous transformation matrices
T_n_s = zeros(4, 4, N);  % Ship pose in NED
T_n_b = zeros(4, 4, N);  % Robot base in NED

% Preallocate position and rotation arrays for plotting
p_s_n = zeros(3, N);    % Ship position in NED
p_b_n = zeros(3, N);    % Robot base position in NED
R_sn = zeros(3, 3, N);  % Ship rotation matrices (NED to ship)
R_bn = zeros(3, 3, N);  % Robot base rotation matrices (NED to base)

%% Manipulator Kinematics (3-DOF: RRP configuration)
% Joint parameters
% q1: Revolute joint 1 (rotation about x-axis)
% q2: Revolute joint 2 (rotation about y-axis)
% d3: Prismatic joint 3 (translation along z-axis)
% c3: Constant offset on prismatic joint (end-effector middel position)

c3 = 6; % Constant offset [m]

% Motion compensation: joint trajectories will be computed to stabilize end-effector
% Preallocate joint trajectory arrays
q1_traj = zeros(N,1);  % Joint 1: rotation about x-axis (compensates roll)
q2_traj = zeros(N,1);  % Joint 2: rotation about y-axis (compensates pitch)
d3_traj = zeros(N,1);  % Joint 3: prismatic along z-axis (compensates heave)

% Preallocate end-effector data
T_n_e = zeros(4, 4, N);  % End-effector pose in NED frame
p_e_n = zeros(3, N);      % End-effector position in NED
R_en = zeros(3, 3, N);    % End-effector rotation (NED to end-effector)

% Define fixed transformation: robot base pose in ship frame
base_orientation = deg2rad(180); % 180 degrees rotation about z-axis
R_bs_x = diag([1, -1, -1]); % 180° rotation about x-axis
R_bs_z = [cos(base_orientation), -sin(base_orientation), 0;
          sin(base_orientation),  cos(base_orientation), 0;
          0,                     0,                    1];
R_bs = R_bs_z * R_bs_x; % Combined rotation from ship to robot base

% Target end-effector position in NED (desired stable position)
phi_0 = eta(4, 1);
theta_0 = eta(5, 1);
psi_0 = eta(6, 1);
R_sn_0 = Rzyx(phi_0, theta_0, psi_0);
p_s_n_0 = eta(1:3, 1);

% Initial robot base position in NED
R_ns_0 = R_sn_0.';
p_b_n_0 = p_s_n_0 + R_ns_0 * p_b_s;

% Initial end-effector position (q1=0, q2=0, d3=c3)
R_bn_0 = R_bs * R_sn_0;
R_nb_0 = R_bn_0.';
p_e_b_0 = [0; 0; c3]; 
p_e_n_target = p_b_n_0 + R_nb_0 * p_e_b_0;

% Robot base pose expressed in ship frame
T_s_b = [R_bs,    p_b_s;
         0 0 0,   1    ];

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

    % Ship pose in NED frame
    R_ns = R_sn(:,:,i).';
    T_n_s(:,:,i) = [R_ns,     p_s_n(:,i);
                    0 0 0,    1          ];

    % Robot base pose in NED frame
    T_n_b(:,:,i) = T_n_s(:,:,i) * T_s_b;

    % Extract rotation (NED to robot base) and position
    R_nb = T_n_b(1:3, 1:3, i);  % Rotation from base to NED
    R_bn(:,:,i) = R_nb.';        % Rotation from NED to base
    p_b_n(:, i) = T_n_b(1:3, 4, i);  % Robot base origin in NED

    %% MOTION COMPENSATION
    % Counteract ship orientation to keep end-effector aligned with NED
    q1 = (phi - phi_0)*cos(base_orientation) - (theta - theta_0)*sin(base_orientation);
    q2 = (phi - phi_0)*sin(base_orientation) + (theta - theta_0)*cos(base_orientation);

    % Target position relative to robot base in NED frame
    p_target_rel = p_e_n_target - p_b_n(:, i);

    % Transform to robot base frame to find required z-extension
    p_target_base = R_bn(:,:,i) * p_target_rel;
        
    % Compute rotation matrix from base to end-effector (first two joints)
    R_12 = [1,  0,        0;
            0,  cos(q1), -sin(q1);
            0,  sin(q1),  cos(q1)] * ...
           [cos(q2),  0,  sin(q2);
            0,        1,  0;
           -sin(q2),  0,  cos(q2)];

    % Rotate target position into frame 2 to isolate prismatic joint
    p_target_rot = R_12.' * p_target_base;
    d3 = p_target_rot(3) - c3;

    % Store computed joint values
    q1_traj(i) = q1;
    q2_traj(i) = q2;
    d3_traj(i) = d3;

    % Frame 1 relative to robot base (rotation about x-axis)
    R_1b = [1,  0,        0;
            0,  cos(q1), -sin(q1);
            0,  sin(q1),  cos(q1)];
    p_1b = [0; 0; 0];
    T_b1 = [R_1b,    p_1b;
            0 0 0,   1    ];

    % Frame 2 relative to frame 1 (rotation about y-axis)
    R_21 = [cos(q2),  0,  sin(q2);
            0,        1,  0;
           -sin(q2),  0,  cos(q2)];
    p_21 = [0; 0; 0];
    T_12 = [R_21,    p_21;
            0 0 0,   1    ];

    % Frame 3 (end-effector) relative to frame 2 (prismatic joint)
    R_32 = eye(3);
    p_32 = [0; 0; c3 + d3];
    T_23 = [R_32,    p_32;
            0 0 0,   1    ];

    % Forward kinematics: end-effector relative to robot base
    T_be = T_b1 * T_12 * T_23;

    % End-effector pose in NED frame
    T_n_e(:,:,i) = T_n_b(:,:,i) * T_be;

    % Extract position and rotation
    p_e_n(:, i) = T_n_e(1:3, 4, i);
    R_ne = T_n_e(1:3, 1:3, i);
    R_en(:,:,i) = R_ne.';
end

% Compute end-effector position error (deviation from target)
p_e_error = p_e_n - repmat(p_e_n_target, 1, N);
p_e_error_norm = vecnorm(p_e_error);

%% Extract Euler angles for plotting
phi_deg = rad2deg(eta(4, :));
theta_deg = rad2deg(eta(5, :));
psi_deg = rad2deg(eta(6, :));

%% Plotting - Create tabbed figure interface

% Create main figure window with tabs
main_fig = figure('Name', 'MCE Platform Simulation with Manipulator', ...
                  'Position', [50 50 1400 900], ...
                  'NumberTitle', 'off');
tab_group = uitabgroup(main_fig);

%% Tab 1: Trajectories in NED Frame
tab1 = uitab(tab_group, 'Title', 'Trajectories');

axes('Parent', tab1);
subplot(2, 2, 1);
plot3(p_s_n(1,:), p_s_n(2,:), p_s_n(3,:), 'b-', 'LineWidth', 2);
hold on;
plot3(p_b_n(1,:), p_b_n(2,:), p_b_n(3,:), 'r-', 'LineWidth', 2);
plot3(p_e_n(1,:), p_e_n(2,:), p_e_n(3,:), 'm-', 'LineWidth', 2);
grid on;
xlabel('North [m]');
ylabel('East [m]');
zlabel('Down [m]');
title('3D Trajectories in NED Frame');
legend('Ship Origin', 'Robot Base', 'End-Effector', 'Location', 'best');
axis equal;
set(gca, 'ZDir', 'reverse', 'YDir', 'reverse');

% Add start and end markers
plot3(p_s_n(1,1), p_s_n(2,1), p_s_n(3,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot3(p_s_n(1,end), p_s_n(2,end), p_s_n(3,end), 'gs', 'MarkerSize', 10, 'LineWidth', 2);
plot3(p_b_n(1,1), p_b_n(2,1), p_b_n(3,1), 'co', 'MarkerSize', 10, 'LineWidth', 2);
plot3(p_e_n(1,1), p_e_n(2,1), p_e_n(3,1), 'mo', 'MarkerSize', 10, 'LineWidth', 2);

subplot(2, 2, 2);
plot(p_s_n(1,:), p_s_n(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(p_b_n(1,:), p_b_n(2,:), 'r-', 'LineWidth', 2);
plot(p_e_n(1,:), p_e_n(2,:), 'm-', 'LineWidth', 2);
plot(p_s_n(1,1), p_s_n(2,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(p_s_n(1,end), p_s_n(2,end), 'gs', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('North [m]');
ylabel('East [m]');
title('Top View (North-East Plane)');
legend('Ship Origin', 'Robot Base', 'End-Effector', 'Start', 'End', 'Location', 'best');
axis equal;

subplot(2, 2, 3);
plot(t, p_s_n(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(t, p_b_n(3,:), 'r-', 'LineWidth', 2);
plot(t, p_e_n(3,:), 'm-', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Down [m]');
title('Vertical Position (Heave)');
legend('Ship Origin', 'Robot Base', 'End-Effector', 'Location', 'best');

subplot(2, 2, 4);
plot(t, vecnorm(p_s_n(1:2,:) - p_s_n(1:2,1)), 'b-', 'LineWidth', 2);
hold on;
plot(t, vecnorm(p_b_n(1:2,:) - p_b_n(1:2,1)), 'r-', 'LineWidth', 2);
plot(t, vecnorm(p_e_n(1:2,:) - p_e_n(1:2,1)), 'm-', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Horizontal Distance [m]');
title('Horizontal Displacement from Start');
legend('Ship Origin', 'Robot Base', 'End-Effector', 'Location', 'best');

%% Tab 2: Ship Orientation
tab2 = uitab(tab_group, 'Title', 'Ship Orientation');
axes('Parent', tab2);

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

%% Tab 3: Manipulator Joint Angles (Motion Compensation)
tab3 = uitab(tab_group, 'Title', 'Joint Angles');
axes('Parent', tab3);

subplot(3, 1, 1);
plot(t, rad2deg(q1_traj), 'LineWidth', 2);
hold on;
plot(t, phi_deg, '--', 'LineWidth', 1.5);
grid on;
ylabel('Angle [deg]');
title('Manipulator Joint Angles (Motion Compensation)');
legend('q_1 (compensation)', '\phi (ship roll)', 'Location', 'best');

subplot(3, 1, 2);
plot(t, rad2deg(q2_traj), 'LineWidth', 2);
hold on;
plot(t, theta_deg, '--', 'LineWidth', 1.5);
grid on;
ylabel('Angle [deg]');
legend('q_2 (compensation)', '\theta (ship pitch)', 'Location', 'best');

subplot(3, 1, 3);
plot(t, d3_traj, 'LineWidth', 2);
hold on;
plot(t, p_s_n(3,:) - p_s_n(3,1), '--', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Extension [m]');
legend('d_3 (compensation)', 'Ship heave', 'Location', 'best');

%% Tab 4: Motion Compensation Performance
tab4 = uitab(tab_group, 'Title', 'Compensation Performance');
axes('Parent', tab4);

subplot(2, 2, 1);
plot(t, p_e_error(1,:), 'LineWidth', 2);
grid on;
ylabel('North Error [m]');
title('End-Effector Position Error (Deviation from Target)');

subplot(2, 2, 2);
plot(t, p_e_error(2,:), 'LineWidth', 2);
grid on;
ylabel('East Error [m]');

subplot(2, 2, 3);
plot(t, p_e_error(3,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Down Error [m]');

subplot(2, 2, 4);
plot(t, p_e_error_norm * 1000, 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Total Position Error [mm]');
title(sprintf('RMS Error: %.2f mm', rms(p_e_error_norm)*1000));

%% Tabs 5-9: Coordinate Frames at Different Time Instants

% Select time instants to visualize
time_indices = round(linspace(1, N, 5));
frame_scale = 8;  % Length of frame axes [m]
time_labels = {'Start', '25%', '50%', '75%', 'End'};

for idx = 1:length(time_indices)
    % Create a tab for this time instant (tabs 5-9)
    tab_frame = uitab(tab_group, 'Title', sprintf('Frames: %s', time_labels{idx}));
    axes('Parent', tab_frame);

    i = time_indices(idx);

    % Plot trajectories up to current time
    plot3(p_s_n(1,1:i), p_s_n(2,1:i), p_s_n(3,1:i), 'b-', 'LineWidth', 1.5);
    hold on;
    plot3(p_b_n(1,1:i), p_b_n(2,1:i), p_b_n(3,1:i), 'r-', 'LineWidth', 1.5);
    plot3(p_e_n(1,1:i), p_e_n(2,1:i), p_e_n(3,1:i), 'm-', 'LineWidth', 1.5);

    % Ship frame axes (extract from T_n_s)
    ship_origin = T_n_s(1:3, 4, i);
    R_ns = T_n_s(1:3, 1:3, i);
    x_s_axis = R_ns * [frame_scale; 0; 0];
    y_s_axis = R_ns * [0; frame_scale; 0];
    z_s_axis = R_ns * [0; 0; frame_scale];

    quiver3(ship_origin(1), ship_origin(2), ship_origin(3), ...
            x_s_axis(1), x_s_axis(2), x_s_axis(3), 0, 'b-', 'LineWidth', 3);
    quiver3(ship_origin(1), ship_origin(2), ship_origin(3), ...
            y_s_axis(1), y_s_axis(2), y_s_axis(3), 0, 'g-', 'LineWidth', 3);
    quiver3(ship_origin(1), ship_origin(2), ship_origin(3), ...
            z_s_axis(1), z_s_axis(2), z_s_axis(3), 0, 'r-', 'LineWidth', 3);

    % Robot base frame axes (extract from T_n_b)
    robot_origin = T_n_b(1:3, 4, i);
    R_nb = T_n_b(1:3, 1:3, i);
    x_b_axis = R_nb * [frame_scale; 0; 0];
    y_b_axis = R_nb * [0; frame_scale; 0];
    z_b_axis = R_nb * [0; 0; frame_scale];

    quiver3(robot_origin(1), robot_origin(2), robot_origin(3), ...
            x_b_axis(1), x_b_axis(2), x_b_axis(3), 0, 'b--', 'LineWidth', 2.5);
    quiver3(robot_origin(1), robot_origin(2), robot_origin(3), ...
            y_b_axis(1), y_b_axis(2), y_b_axis(3), 0, 'g--', 'LineWidth', 2.5);
    quiver3(robot_origin(1), robot_origin(2), robot_origin(3), ...
            z_b_axis(1), z_b_axis(2), z_b_axis(3), 0, 'r--', 'LineWidth', 2.5);

    % End-effector frame axes (extract from T_n_e)
    ee_origin = T_n_e(1:3, 4, i);
    R_ne = T_n_e(1:3, 1:3, i);
    x_e_axis = R_ne * [frame_scale*0.6; 0; 0];
    y_e_axis = R_ne * [0; frame_scale*0.6; 0];
    z_e_axis = R_ne * [0; 0; frame_scale*0.6];

    quiver3(ee_origin(1), ee_origin(2), ee_origin(3), ...
            x_e_axis(1), x_e_axis(2), x_e_axis(3), 0, 'b:', 'LineWidth', 2);
    quiver3(ee_origin(1), ee_origin(2), ee_origin(3), ...
            y_e_axis(1), y_e_axis(2), y_e_axis(3), 0, 'g:', 'LineWidth', 2);
    quiver3(ee_origin(1), ee_origin(2), ee_origin(3), ...
            z_e_axis(1), z_e_axis(2), z_e_axis(3), 0, 'r:', 'LineWidth', 2);

    % Draw manipulator links
    plot3([robot_origin(1), ee_origin(1)], ...
          [robot_origin(2), ee_origin(2)], ...
          [robot_origin(3), ee_origin(3)], 'k-', 'LineWidth', 3);

    % Add origin markers
    plot3(ship_origin(1), ship_origin(2), ship_origin(3), 'bo', ...
          'MarkerSize', 12, 'MarkerFaceColor', 'b');
    plot3(robot_origin(1), robot_origin(2), robot_origin(3), 'ro', ...
          'MarkerSize', 12, 'MarkerFaceColor', 'r');
    plot3(ee_origin(1), ee_origin(2), ee_origin(3), 'mo', ...
          'MarkerSize', 12, 'MarkerFaceColor', 'm');

    grid on;
    xlabel('North [m]', 'FontSize', 12);
    ylabel('East [m]', 'FontSize', 12);
    zlabel('Down [m]', 'FontSize', 12);
    title(sprintf('Coordinate Frames at t = %.2f s (q1=%.1f°, q2=%.1f°, d3=%.2fm)', ...
          t(i), rad2deg(q1_traj(i)), rad2deg(q2_traj(i)), d3_traj(i)), ...
          'FontSize', 14, 'FontWeight', 'bold');
    axis equal;
    set(gca, 'ZDir', 'reverse', 'YDir', 'reverse', 'FontSize', 11);
    view(45, 20);

    legend('Ship Traj', 'Robot Traj', 'EE Traj', ...
           'x_s (fwd)', 'y_s (stbd)', 'z_s (down)', ...
           'x_b (back)', 'y_b (stbd)', 'z_b (up)', ...
           'x_e', 'y_e', 'z_e', 'Manipulator', ...
           'Ship Origin', 'Robot Origin', 'End-Effector', ...
           'Location', 'eastoutside', 'FontSize', 10);
end

%% Tab 10: Verification - Robot Base in Ship Frame
tab10 = uitab(tab_group, 'Title', 'Verification');
axes('Parent', tab10);

% Transform robot base position back to ship frame using homogeneous transforms
% T_n_s transforms from ship to NED, so we need its inverse
p_b_s_check = zeros(3, N);
for i = 1:N
    % Inverse of T_n_s to get T_s_n (NED to ship transform)
    % For T = [R, p; 0 1], T^-1 = [R', -R'*p; 0 1]
    R_ns = T_n_s(1:3, 1:3, i);
    p_s_n_i = T_n_s(1:3, 4, i);

    T_s_n = [R_ns.',  -R_ns.' * p_s_n_i;
             0 0 0,   1                 ];

    % Transform robot origin from NED to ship frame
    p_b_n_homo = [p_b_n(:, i); 1];
    p_b_s_homo = T_s_n * p_b_n_homo;
    p_b_s_check(:, i) = p_b_s_homo(1:3);
end

subplot(1, 3, 1);
plot(t, p_b_s_check(1,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('x_s [m]');
title('Robot Base x-position in Ship Frame (should be constant)');
ylim([p_b_s(1)-1, p_b_s(1)+1]);

subplot(1, 3, 2);
plot(t, p_b_s_check(2,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('y_s [m]');
title('Robot Base y-position in Ship Frame (should be constant)');
ylim([p_b_s(2)-1, p_b_s(2)+1]);

subplot(1, 3, 3);
plot(t, p_b_s_check(3,:), 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('z_s [m]');
title('Robot Base z-position in Ship Frame (should be constant)');
ylim([p_b_s(3)-1, p_b_s(3)+1]);

%% Display statistics
fprintf('\n=== Simulation Statistics ===\n');
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
fprintf('\n--- End-Effector Motion Range (in NED) ---\n');
fprintf('Position (NED):\n');
fprintf('  North: [%.2f, %.2f] m\n', min(p_e_n(1,:)), max(p_e_n(1,:)));
fprintf('  East:  [%.2f, %.2f] m\n', min(p_e_n(2,:)), max(p_e_n(2,:)));
fprintf('  Down:  [%.2f, %.2f] m\n', min(p_e_n(3,:)), max(p_e_n(3,:)));
fprintf('\n--- Manipulator Joint Ranges (Motion Compensation) ---\n');
fprintf('Joint angles:\n');
fprintf('  q1: [%.2f, %.2f] deg (vs ship roll: [%.2f, %.2f] deg)\n', ...
        min(rad2deg(q1_traj)), max(rad2deg(q1_traj)), min(phi_deg), max(phi_deg));
fprintf('  q2: [%.2f, %.2f] deg (vs ship pitch: [%.2f, %.2f] deg)\n', ...
        min(rad2deg(q2_traj)), max(rad2deg(q2_traj)), min(theta_deg), max(theta_deg));
fprintf('  d3: [%.2f, %.2f] m (vs ship heave: [%.2f, %.2f] m)\n', ...
        min(d3_traj), max(d3_traj), ...
        min(p_s_n(3,:) - p_s_n(3,1)), max(p_s_n(3,:) - p_s_n(3,1)));
fprintf('\n--- Motion Compensation Performance ---\n');
fprintf('Target end-effector position: [%.2f, %.2f, %.2f] m\n', p_e_n_target);
fprintf('Actual end-effector position range:\n');
fprintf('  North: [%.2f, %.2f] m (error: %.3f m)\n', ...
        min(p_e_n(1,:)), max(p_e_n(1,:)), max(abs(p_e_error(1,:))));
fprintf('  East:  [%.2f, %.2f] m (error: %.3f m)\n', ...
        min(p_e_n(2,:)), max(p_e_n(2,:)), max(abs(p_e_error(2,:))));
fprintf('  Down:  [%.2f, %.2f] m (error: %.3f m)\n', ...
        min(p_e_n(3,:)), max(p_e_n(3,:)), max(abs(p_e_error(3,:))));
fprintf('Position error statistics:\n');
fprintf('  RMS error: %.2f mm\n', rms(p_e_error_norm) * 1000);
fprintf('  Max error: %.2f mm\n', max(p_e_error_norm) * 1000);
fprintf('  Mean error: %.2f mm\n', mean(p_e_error_norm) * 1000);
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