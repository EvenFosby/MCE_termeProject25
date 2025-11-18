clear; clc; close all;

load('vessel_motion_data.mat'); % Load vessel motion data

% Simulation parameters
t = motion_data.time; % Time vector
dt = t(2) - t(1); % Time step
N = length(t); % Number of time steps

% Extract vessel motion data
eta = motion_data.eta'; % Vessel position and orientation
nu = motion_data.nu'; % Vessel velocities

% Manipulator base position in body frame
p_0_b = [-30; 0; -3]; % Example position of manipulator base in body frame

% Preallocate arrays
R_bn = zeros(3, 3, N);  % Ship to NED rotation matrices
R_0n = zeros(3, 3, N);  % Robot base to NED rotation matrices
p_b_n = zeros(3, N);    % Ship position in NED
p_0_n = zeros(3, N);    % Robot base position in NED

% Constant rotations
    R_0b = Rzyx(0, pi, 0);

for i = 1:N
    % Extract vessel position
    x = eta(1, i);
    y = eta(2, i);
    z = eta(3, i);
    % Extract Euler angles
    phi = eta(4, i);   % Roll
    theta = eta(5, i); % Pitch
    psi = eta(6, i);   % Yaw

    % Ship body frame to NED frame
    p_b_n(:, :, i) = [x; y; z];
    R_bn(:, :, i) = Rzyx(phi, theta, psi);

    T_bn = [R_bn(:, :, i), p_b_n(:, :, i);
            eye(1,3),                   1];

    % Robot base frame to NED frame
    R_0n(:, :, i) = R_bn(:, :, i) * R_0b;

    % Robot base position in NED frame
    p_0_n(:, :, i) = p_b_n(:, :, i) + R_bn(:, :, i) * p_0_b;

    % NED to manipulator base frame

end

q1 = 0; % rad
q2 = 0; % rad
d3 = 0; c3 = 6; % m

R_1b = [1 0 0; 
        0 cos(q1) -sin(q1);
        0 sin(q1)  cos(q1)];

p_1_b = [0; 0; 0];

T_1b = [R_1b,    p_1_b;
         0 0 0,   1    ];

R_21 = [cos(q2) 0 sin(q2);
        0       1 0;
       -sin(q2) 0 cos(q2)];
p_21 = [0; 0; 0];

T_21 = [R_21,    p_21;
         0 0 0,   1    ];

R_32 = [1 0 0;
        0 1 0;
        0 0 1];
p_32 = [0; 0; c3 + d3];

T_32 = [R_32,    p_32;
         0 0 0,   1    ];


function R = Rzyx(phi, theta, psi)
    % RZYX Create rotation matrix from ZYX Euler angles (roll-pitch-yaw)
    % Input: phi (roll), theta (pitch), psi (yaw) in radians
    % Output: R = Rz(psi) * Ry(theta) * Rx(phi)

    cpsi = cos(psi);   spsi = sin(psi);
    cth = cos(theta);  sth = sin(theta);
    cphi = cos(phi);   sphi = sin(phi);

    R = [cpsi*cth, -spsi*cphi + cpsi*sth*sphi,  spsi*sphi + cpsi*sth*cphi;
         spsi*cth,  cpsi*cphi + spsi*sth*sphi, -cpsi*sphi + spsi*sth*cphi;
         -sth,      cth*sphi,                   cth*cphi];
end

