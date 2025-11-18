function p_b = transform_point_to_robot_frame(p_n, R_bn, p_b_n)
% TRANSFORM_POINT_TO_ROBOT_FRAME Transform a point from NED to robot base frame
%
% Inputs:
%   p_n   : 3x1 or 3xN point(s) expressed in NED frame [m]
%   R_bn  : 3x3 rotation matrix from NED to robot base frame
%   p_b_n : 3x1 robot base origin expressed in NED frame [m]
%
% Output:
%   p_b   : 3x1 or 3xN point(s) expressed in robot base frame [m]
%
% Formula:
%   p_b = R_bn * (p_n - p_b_n)

    % Handle single point or multiple points
    if size(p_n, 2) == 1
        % Single point
        p_b = R_bn * (p_n - p_b_n);
    else
        % Multiple points
        N = size(p_n, 2);
        p_b = zeros(3, N);
        for i = 1:N
            p_b(:,i) = R_bn * (p_n(:,i) - p_b_n);
        end
    end

end
