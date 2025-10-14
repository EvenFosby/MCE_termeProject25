% Motion-Stabilized Platform Simulation
clc; clear; close all;

%% Parameters
L3 = 1.0;        % Fixed length of cylinder (meters)
d3 = 0.5;        % Extendable length (meters)
theta1 = deg2rad(20);  % Roll angle (degrees → radians)
theta2 = deg2rad(15);  % Pitch angle (degrees → radians)

%% DH Parameters: [theta, d, a, alpha]
DH = [theta1, 0, 0, -pi/2;
      theta2, 0, 0,  pi/2;
      0,      L3 + d3, 0, 0];

%% Function to compute DH transformation
dh_transform = @(theta, d, a, alpha) [
    cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
    sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
    0,           sin(alpha),             cos(alpha),            d;
    0,           0,                      0,                     1];

%% Forward Kinematics
T01 = dh_transform(DH(1,1), DH(1,2), DH(1,3), DH(1,4));
T12 = dh_transform(DH(2,1), DH(2,2), DH(2,3), DH(2,4));
T23 = dh_transform(DH(3,1), DH(3,2), DH(3,3), DH(3,4));

T02 = T01 * T12;
T03 = T02 * T23;

%% Extract positions
origin = [0; 0; 0];
p1 = T01(1:3,4);
p2 = T02(1:3,4);
p3 = T03(1:3,4);  % End effector

%% Plotting
figure;
plot3([origin(1) p1(1) p2(1) p3(1)], ...
      [origin(2) p1(2) p2(2) p3(2)], ...
      [origin(3) p1(3) p2(3) p3(3)], 'o-', 'LineWidth', 2);
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Motion-Stabilized Platform Forward Kinematics');
axis equal;
view(135, 30);

% Show coordinate frames
hold on;
frames = {eye(4), T01, T02, T03};
colors = {'r','g','b'};
for i = 1:length(frames)
    T = frames{i};
    origin = T(1:3,4);
    for j = 1:3
        dir = T(1:3,j);
        quiver3(origin(1), origin(2), origin(3), ...
                dir(1), dir(2), dir(3), 0.2, colors{j}, 'LineWidth', 1.5);
    end
end
legend('Links', 'X-axis', 'Y-axis', 'Z-axis');
