%% EXAMPLE: Frame Transformation from NED to Robot Base
% This script demonstrates how to use the coordinate frame transformation
% functions for a robot mounted on a ship deck.

clear; clc; close all;

%% Example 1: Basic transformation at a single time instant

% Ship orientation (example: heading 30°, roll 5°, pitch 2°)
psi = deg2rad(30);   % Yaw (heading)
phi = deg2rad(5);    % Roll
theta = deg2rad(2);  % Pitch

% Create rotation matrix from NED to ship frame using ZYX Euler angles
R_sn = Rzyx(phi, theta, psi);

% Ship position in NED frame [North, East, Down] in meters
p_s_n = [100; 50; 5];

% Robot base position in ship frame [x_s, y_s, z_s] in meters
% Example: robot base is 10m forward, 2m to starboard, 3m above ship origin
p_b_s = [10; 2; -3];  % z_s is negative because ship z points down

%% Compute transformations

% Get rotation and position of robot base in NED frame
[R_bn, p_b_n] = ned_to_robot_base(R_sn, p_s_n, p_b_s);

fprintf('Ship position in NED: [%.2f, %.2f, %.2f] m\n', p_s_n);
fprintf('Robot base position in NED: [%.2f, %.2f, %.2f] m\n', p_b_n);
fprintf('\n');

%% Example 2: Transform a point from NED to robot base frame

% Point of interest in NED frame (e.g., target position)
p_target_n = [120; 55; 10];

% Transform to robot base coordinates
p_target_b = transform_point_to_robot_frame(p_target_n, R_bn, p_b_n);

fprintf('Target in NED frame: [%.2f, %.2f, %.2f] m\n', p_target_n);
fprintf('Target in robot base frame: [%.2f, %.2f, %.2f] m\n', p_target_b);
fprintf('\n');

%% Example 3: Using homogeneous transformations

% Get 4x4 homogeneous transformation matrix
T_bn = ned_to_robot_base_homogeneous(R_sn, p_s_n, p_b_s);

% Transform point using homogeneous coordinates
p_target_n_hom = [p_target_n; 1];
p_target_b_hom = T_bn * p_target_n_hom;
p_target_b_check = p_target_b_hom(1:3);

fprintf('Target in robot base (homogeneous method): [%.2f, %.2f, %.2f] m\n', ...
        p_target_b_check);
fprintf('Difference from direct method: %.2e m\n', norm(p_target_b - p_target_b_check));
fprintf('\n');

%% Example 4: Time series transformation (for ship motion simulation)

% Simulate 100 time steps
N = 100;
t = linspace(0, 10, N);

% Preallocate arrays
R_bn_series = zeros(3, 3, N);
p_b_n_series = zeros(3, N);

% Example: ship following a sinusoidal motion
for i = 1:N
    % Time-varying ship orientation
    psi_t = deg2rad(30) + 0.1*sin(2*pi*0.1*t(i));    % Yaw oscillation
    phi_t = deg2rad(5)*sin(2*pi*0.2*t(i));           % Roll oscillation
    theta_t = deg2rad(2)*sin(2*pi*0.15*t(i));        % Pitch oscillation

    R_sn_t = Rzyx(phi_t, theta_t, psi_t);

    % Time-varying ship position (moving north)
    p_s_n_t = [100 + 2*t(i); 50; 5 + 0.5*sin(2*pi*0.3*t(i))];

    % Robot base fixed in ship frame
    p_b_s_fixed = [10; 2; -3];

    % Compute robot base frame transformation
    [R_bn_series(:,:,i), p_b_n_series(:,i)] = ...
        ned_to_robot_base(R_sn_t, p_s_n_t, p_b_s_fixed);
end

% Plot robot base trajectory in NED frame
figure('Name', 'Robot Base Trajectory in NED Frame');
plot3(p_b_n_series(1,:), p_b_n_series(2,:), p_b_n_series(3,:), 'b-', 'LineWidth', 2);
grid on;
xlabel('North [m]');
ylabel('East [m]');
zlabel('Down [m]');
title('Robot Base Trajectory in NED Frame');
axis equal;
set(gca, 'ZDir', 'reverse');  % Down is positive

%% Helper function: ZYX Euler angle rotation matrix
function R = Rzyx(phi, theta, psi)
    % RZYX Create rotation matrix from ZYX Euler angles (yaw-pitch-roll)
    % R rotates from inertial frame to body frame
    % phi: roll, theta: pitch, psi: yaw

    cpsi = cos(psi);   spsi = sin(psi);
    cth = cos(theta);  sth = sin(theta);
    cphi = cos(phi);   sphi = sin(phi);

    R = [cpsi*cth, -spsi*cphi + cpsi*sth*sphi,  spsi*sphi + cpsi*sth*cphi;
         spsi*cth,  cpsi*cphi + spsi*sth*sphi, -cpsi*sphi + spsi*sth*cphi;
         -sth,      cth*sphi,                   cth*cphi];
end
