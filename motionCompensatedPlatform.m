function [p_e_n, v_e_n, R_e_n] = motionCompensatedPlatform(eta_vessel, nu_vessel, q_manip, qdot_manip, vessel_param)
% MOTIONCOMPENSATEDPLATFORM Computes the end-effector pose and velocity of 
% an RRP manipulator mounted on a DP vessel in wave motion. 
%
% Inputs: 
%   eta_vessel - [6x1] Vessel pose in NED: [x, y, z, phi, theta, psi]
%   nu_vessel    - [6x1] Vessel velocities in body frame: [u, v, w, p, q, r]
%   q_manip      - [3x1] Manipulator joint coordinates: [q1, q2, d3]
%                  q1: azimuth rotation (yaw), q2: elevation (pitch), d3: extension
%   qdot_manip   - [3x1] Manipulator joint velocities
%   r_0b         - [3x1] Position of platform base {0} relative to vessel {b}
%   vessel_params - struct with vessel mass properties (optional, for dynamics)
%
% Outputs:
%   p_e_n - [3x1] End-effector position in NED frame
%   v_e_n - [3x1] End-effector velocity in NED frame
%   R_e_n - [3x3] End-effector rotation matrix (NED to end-effector)

% Extract vessel motion
x_n     = eta_vessel(1:3);  % Position {NED}
phi     = eta_vessel(4);    % Roll
theta   = eta_vessel(5);    % Pitch
psi     = eta_vessel(6);    % Yaw

v_b         = nu_vessel(1:3);   % Linear velocity {b}
omega_b     = nu_vessel(4:6);   % Angular velocity {b}

% Rotation from {b} to {NED}
R_nb = Rzyx(phi, theta, psi);

% Manipulator config
q1 = q_manip(1);    % Rotation angle about x0 (Roll)
q2 = q_manip(2);    % Rotation angle about y1 (Pitch)
q3 = q_manip(3);    % Extention along z2 (Heave)

q1dot = qdot_manip(1);
q2dot = qdot_manip(2);
q3dot = qdot_manip(3);

% Manipulator constants
d3 = 1; % cylinder fixed lenght

%% Forward Kinematics 
%  Platform to end-effector

% Rotation about x0
Rx_q1 = [1   0        0;
         0   cos(q1) -sin(q1);
         0   sin(q1)  cos(q1)];

% Rotation about y1
Ry_q2 = [cos(q2)    0   sin(q2);
         0          1   0;
        -sin(q2)    0   cos(q2)];

% Rotation form base {0} to end-effector {ee}
R_e0 = Rx_q1 * Ry_q2;

% End-effector position relativ to base
p_e0 = R_e0 * [0 0 d3+q3]';

% Position in vessel body frame
p_eb = r_0b + p_e0;

% Rotation from NED to EE
R_en = R_nb * R_e0;

% Position in NED frame
p_en = x_n + R_nb * R_e0;

%% Velocity Kinematics
% Manipulator Jacobian
J = JacobianRRP(q1,q2,q3);

% End-effector velocity
V_e0_0  = J * qdot_manip;               % {0} to {e} in {0} frame
V_e0_b  = R_e0 * V_e0_0;                % {0} to {e} in {b} frame

V_vessel_contrib = v_b + cross(omega_b, p_eb);
V_eb    = V_vessel_contrib + V_e0_b;    % Velocity in body

V_en    = R_nb *V_eb;                   % Velocity in NED
end


function J = JacobianRRP(q1, q2, q3)
% ∂p/∂q1
dp_dq1 = [0;
         -q3*cos(q1)*cos(q2);
         -q3*sin(q1)*cos(q2)];
% ∂p/∂q2
dp_dq2 = [q3*cos(q2); 
          q3*sin(q1)*sin(q2);
         -q3*cos(q1)*sin(q2)];
% ∂p/∂d3
dp_dq3 = [sin(q2);
         -sin(q1)*cos(q2);
          cos(q1)*cos(q2)];

J = [dp_dq1 dp_dq2 dp_dq3];
end

% function J_omega = angVelJacobian(q1, q2)
% R_x_q1 = [1   0         0;
%           0   cos(q1)  -sin(q1);
%           0   sin(q1)   cos(q1)];
% 
% x0 = [1; 0; 0];
% y1 = R_x_q1 * [0; 1; 0];
% 
% J_omega = [x0, y1];
% end