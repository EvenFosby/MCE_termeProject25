
% ---------- manipulator_plot.m ----------
% Plot a 3-DOF manipulator defined by:
% T = [ R_x(a)*R_y(b),  R_x(a)*R_y(b)*[0;0;c+d];
%       0 0 0,          1 ]
% where a=q1, b=q2 (revolute), c=q3 (prismatic), d=constant offset.

% Example pose (radians, meters)
a = deg2rad(50);   % q1
b = deg2rad(20);  % q2
c = 0.25;   % q3 (prismatic extension)
d = 3.5;   % constant offset

manipulator_plot(a,b,c,d);

function manipulator_plot(a, b, c, d)
    if nargin < 4, d = 0.2; end   % default constant offset (meters)
    if nargin < 3, error('Usage: manipulator_plot(a,b,c,dOptional)'); end

    % Rotation matrices
    Rx = @(th)[1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
    Ry = @(th)[cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];

    R = Rx(a) * Ry(b);
    p = R * [0;0;(c + d)];     % end-effector position

    % Plot setup
    figure('Color','w'); clf; hold on; grid on;
    axis equal
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(35,25);

    % Workspace bounds (auto from length)
    L = max(0.1, abs(c)+abs(d)) + 0.2;
    axis([-L L -L L 0 L*1.5]);

    % Draw base frame at origin
    drawFrame(eye(3), [0;0;0], 0.12*max(0.3,L), 'Base');

    % Draw oriented frame at the wrist/end-effector
    drawFrame(R, p, 0.12*max(0.3,L), 'EE');

    % Link is just a prismatic slide along current z after two rotations
    plot3([0 p(1)], [0 p(2)], [0 p(3)], 'LineWidth', 3);

    title(sprintf('3-DOF Manipulator  |  a=q1=%.2f rad, b=q2=%.2f rad, c=q3=%.2f', a,b,c));
    hold off;
end

% ----- Helper: draw a coordinate frame (triad) -----
function drawFrame(R, p, len, labelStr)
    % R: 3x3 rotation, p: 3x1 position, len: axis length
    o = p(:).';
    x = p + R(:,1)*len;
    y = p + R(:,2)*len;
    z = p + R(:,3)*len;

    quiver3(o(1),o(2),o(3), x(1)-o(1), x(2)-o(2), x(3)-o(3), 0, 'LineWidth', 2);
    quiver3(o(1),o(2),o(3), y(1)-o(1), y(2)-o(2), y(3)-o(3), 0, 'LineWidth', 2);
    quiver3(o(1),o(2),o(3), z(1)-o(1), z(2)-o(2), z(3)-o(3), 0, 'LineWidth', 2);

    text(o(1), o(2), o(3), ['  ' labelStr], 'FontWeight','bold');
end
% ---------- end of file ----------


