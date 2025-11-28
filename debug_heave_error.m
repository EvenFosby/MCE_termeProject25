% Debug script to analyze the heave error source
clear; clc;

load('vessel_motion_data.mat');

% Simulation parameters
t = motion_data.time;
N = length(t);
eta = motion_data.eta';

%% Robot base position in ship frame
p_b_s = [-30; 0; -3];
c3 = 6; % Constant offset [m]

% Base orientation
base_orientation = deg2rad(180);
R_bs_x = diag([1, -1, -1]);
R_bs_z = [cos(base_orientation), -sin(base_orientation), 0;
          sin(base_orientation),  cos(base_orientation), 0;
          0,                     0,                    1];
R_bs = R_bs_z * R_bs_x;

% Initial conditions
phi_0 = eta(4, 1);
theta_0 = eta(5, 1);
psi_0 = eta(6, 1);
R_sn_0 = Rzyx(phi_0, theta_0, psi_0);
p_s_n_0 = eta(1:3, 1);

% Initial robot base position in NED
R_ns_0 = R_sn_0.';
p_b_n_0 = p_s_n_0 + R_ns_0 * p_b_s;

% Initial end-effector position (with all joints at zero)
R_bn_0 = R_bs * R_sn_0;
R_nb_0 = R_bn_0.';
p_e_b_0 = [0; 0; c3];
p_e_n_target = p_b_n_0 + R_nb_0 * p_e_b_0;

fprintf('=== Initial Conditions at t=0 ===\n');
fprintf('Target end-effector (NED): [%.4f, %.4f, %.4f] m\n', p_e_n_target);
fprintf('p_e_b_0 (base frame): [%.4f, %.4f, %.4f] m\n', p_e_b_0);
fprintf('\nInitial ship orientation: phi=%.4f, theta=%.4f, psi=%.4f rad\n', phi_0, theta_0, psi_0);

% Robot base pose expressed in ship frame
T_s_b = [R_bs,    p_b_s;
         0 0 0,   1    ];

% Analyze specific time points
test_indices = [1, round(N/4), round(N/2), round(3*N/4), N];
test_labels = {'t=0', 't=25%', 't=50%', 't=75%', 't=end'};

fprintf('\n=== Detailed Analysis at Key Time Points ===\n');

for idx = 1:length(test_indices)
    i = test_indices(idx);

    % Extract ship orientation
    phi = eta(4, i);
    theta = eta(5, i);
    psi = eta(6, i);

    % Ship position in NED
    p_s_n = eta(1:3, i);

    % Ship rotation and transformation
    R_sn = Rzyx(phi, theta, psi);
    R_ns = R_sn.';
    T_n_s = [R_ns,     p_s_n;
             0 0 0,    1     ];

    % Robot base pose in NED frame
    T_n_b = T_n_s * T_s_b;
    R_nb = T_n_b(1:3, 1:3);
    R_bn = R_nb.';
    p_b_n = T_n_b(1:3, 4);

    % Compute compensation angles
    q1 = (phi - phi_0)*cos(base_orientation) - (theta - theta_0)*sin(base_orientation);
    q2 = (phi - phi_0)*sin(base_orientation) + (theta - theta_0)*cos(base_orientation);

    % Target position relative to robot base in NED frame
    p_target_rel = p_e_n_target - p_b_n;

    % Transform to robot base frame
    p_target_base = R_bn * p_target_rel;

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

    % Forward kinematics verification
    R_1b = [1,  0,        0;
            0,  cos(q1), -sin(q1);
            0,  sin(q1),  cos(q1)];
    R_21 = [cos(q2),  0,  sin(q2);
            0,        1,  0;
           -sin(q2),  0,  cos(q2)];
    R_32 = eye(3);
    p_32 = [0; 0; c3 + d3];

    T_b1 = [R_1b,    [0;0;0];
            0 0 0,   1      ];
    T_12 = [R_21,    [0;0;0];
            0 0 0,   1      ];
    T_23 = [R_32,    p_32;
            0 0 0,   1    ];

    % Forward kinematics
    T_be = T_b1 * T_12 * T_23;
    T_n_e = T_n_b * T_be;

    p_e_n = T_n_e(1:3, 4);
    p_e_b_actual = T_be(1:3, 4);

    % Error analysis
    error_ned = p_e_n - p_e_n_target;
    error_base = R_bn * error_ned;

    fprintf('\n--- %s (t=%.2fs) ---\n', test_labels{idx}, t(i));
    fprintf('Ship orientation: phi=%.4f, theta=%.4f (deg: %.2f, %.2f)\n', phi, theta, rad2deg(phi), rad2deg(theta));
    fprintf('Joint angles: q1=%.4f, q2=%.4f, d3=%.4f (deg: %.2f, %.2f)\n', q1, q2, d3, rad2deg(q1), rad2deg(q2));
    fprintf('p_target_base: [%.4f, %.4f, %.4f]\n', p_target_base);
    fprintf('p_target_rot (frame2): [%.4f, %.4f, %.4f]\n', p_target_rot);
    fprintf('p_e_b_actual: [%.4f, %.4f, %.4f]\n', p_e_b_actual);
    fprintf('Target Down (NED): %.4f m\n', p_e_n_target(3));
    fprintf('Actual Down (NED): %.4f m\n', p_e_n(3));
    fprintf('Down error: %.6f m = %.2f mm\n', error_ned(3), error_ned(3)*1000);
    fprintf('Error in base frame: [%.6f, %.6f, %.6f] m\n', error_base);
end

% Helper function
function R = Rzyx(phi, theta, psi)
    cpsi = cos(psi);   spsi = sin(psi);
    cth = cos(theta);  sth = sin(theta);
    cphi = cos(phi);   sphi = sin(phi);

    R = [cpsi*cth, -spsi*cphi + cpsi*sth*sphi,  spsi*sphi + cpsi*sth*cphi;
         spsi*cth,  cpsi*cphi + spsi*sth*sphi, -cpsi*sphi + spsi*sth*cphi;
         -sth,      cth*sphi,                   cth*cphi];
end
