function [theta, phi] = jisuanangle(pos1, pos2)
%JISUANANGLE 计算pos1相对pos2的方位角(theta)和俯仰角(phi)
%   输入：
%       pos1 - 待测点三维坐标，格式为[X1,Y1,Z1]（行向量）
%       pos2 - 参考点三维坐标，格式为[X2,Y2,Z2]（行向量）
%   输出：
%       theta - 方位角（度），范围[-180,180]，以X轴为0度，逆时针为正
%       phi   - 俯仰角（度），范围[-90,90]，向上为正，向下为负

% 计算相对坐标差（从pos2指向pos1的矢量）
p = pos1 - pos2;
dx = p(1);  % X方向差
dy = p(2);  % Y方向差
dz = p(3);  % Z方向差

% 计算水平距离（xy平面），增加微小值避免除以0
r_xy = sqrt(dx^2 + dy^2 + 1e-10);  

% 计算方位角（转换为角度）
theta = atan2(dy, dx) * 180 / pi;

% 计算俯仰角（转换为角度）
phi = atan2(dz, r_xy) * 180 / pi;

end