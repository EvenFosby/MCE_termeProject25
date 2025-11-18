function [R_bn, p_b_n] = ned_to_robot_base(R_sn, p_s_n, p_b_s)
% NED_TO_ROBOT_BASE Compute transformation from NED frame to robot base frame
%
% Inputs:
%   R_sn  : 3x3 rotation matrix from NED (n) to ship body frame (s)
%           Transforms vectors: v_s = R_sn * v_n
%   p_s_n : 3x1 ship origin expressed in NED frame [m]
%   p_b_s : 3x1 robot base origin expressed in ship frame [m]
%
% Outputs:
%   R_bn  : 3x3 rotation matrix from NED to robot base frame (b)
%           Transforms vectors: v_b = R_bn * v_n
%   p_b_n : 3x1 robot base origin expressed in NED frame [m]
%
% Frame conventions:
%   - NED frame (n): North-East-Down
%   - Ship frame (s): x_s forward, y_s starboard, z_s down (SNAME)
%   - Robot base (b): x_b backward (-x_s), y_b starboard (y_s), z_b up (-z_s)
%                     (180° rotation about ship y-axis)

    % Rotation from ship frame to robot base frame
    % R_bs rotates 180° about y-axis: x_b = -x_s, y_b = y_s, z_b = -z_s
    R_bs = diag([-1, 1, -1]);

    % Rotation from NED to robot base frame
    R_bn = R_bs * R_sn;

    % Rotation from ship frame to NED frame (inverse rotation)
    R_ns = R_sn.';

    % Robot base origin in NED frame
    % Transform robot base position from ship frame to NED and add ship origin
    p_b_n = p_s_n + R_ns * p_b_s;

end
