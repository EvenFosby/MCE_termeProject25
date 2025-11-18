function T_bn = ned_to_robot_base_homogeneous(R_sn, p_s_n, p_b_s)
% NED_TO_ROBOT_BASE_HOMOGENEOUS Compute 4x4 homogeneous transform NED to robot base
%
% Inputs:
%   R_sn  : 3x3 rotation matrix from NED (n) to ship body frame (s)
%   p_s_n : 3x1 ship origin expressed in NED frame [m]
%   p_b_s : 3x1 robot base origin expressed in ship frame [m]
%
% Output:
%   T_bn  : 4x4 homogeneous transformation matrix from NED to robot base frame
%           Transforms homogeneous points: p_b_hom = T_bn * p_n_hom
%
% Frame conventions:
%   - NED frame (n): North-East-Down
%   - Ship frame (s): x_s forward, y_s starboard, z_s down (SNAME)
%   - Robot base (b): x_b backward, y_b starboard, z_b up (180° rotation about y)

    % Homogeneous transform from NED to ship frame
    T_sn = [R_sn,    p_s_n;
            0 0 0,   1    ];

    % Rotation from ship to robot base (180° about y-axis)
    R_bs = diag([-1, 1, -1]);

    % Homogeneous transform from ship to robot base
    T_bs = [R_bs,    p_b_s;
            0 0 0,   1    ];

    % Homogeneous transform from NED to robot base
    T_bn = T_bs * T_sn;

end
