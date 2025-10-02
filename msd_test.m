clear; clc; close all;

%% Constants
% Ship Parameters
m = 6000e+3; %ship weight = 6000 ton
L = 90.6; B= 22.0; T = 7.5; 
r_g = [0 0 0];
R44 = 0.35*B;
R55 = 0.25*L;
R66 = 0.25*L;

I_g = m*diag([R44^2 R55^2 R66^2]);
S = Smtrx(r_g);

MRB = [m*eye(3) -m*S;
       m*S      I_g - m*S^2];

D = diag([0 0 0 0 0 0]);
K = diag([0 0 0 0 0 0]);

% external forces
tau0 = [0 0 0 0 0 0]';
tau1 = [0 0 1.0e+3 5.0e+2 4.5e+4 0]'; 
omega = 0;  
phi   = zeros(6,1);


% state space 
A = [zeros(6) eye(6);
    -inv(MRB)*K -inv(MRB)*D];

B = [zeros(6) eye(6)/MRB]';

C = [eye(6) zeros(6)];

D = zeros(6);

sys = ss(A, B, C, D);

% time series
t = (0:0.05:200)';          % column vector
tau = cos(omega*t + phi.'); % N x 6
tau = tau .* tau0.';
tau(1,:) = tau1;
y = lsim(sys, tau, t, zeros(12,1));

% Plot all 6 DOFs in a single page
dof_labels = {'x [m] (surge)','y [m] (sway)','z [m] (heave)', ...
              '\phi [rad] (roll)','\theta [rad] (pitch)','\psi [rad] (yaw)'};

figure;
for i = 1:6
    subplot(3,2,i)                 % 3 rows × 2 cols of subplots
    plot(t, y(:,i), 'LineWidth',1)
    grid on
    xlabel('Time [s]')
    ylabel(dof_labels{i})
    title(['Response of ' dof_labels{i}])
end
sgtitle('6-DOF Ship Motion Responses')
