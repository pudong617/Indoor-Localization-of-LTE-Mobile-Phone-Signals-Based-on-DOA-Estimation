% 结合三个阵列获得的MUSIC谱，联合锁定待测目标三维位置
clc; clear; close all;

% 1. 环境与坐标定义
roomsize = [10, 10, 3]; % 房屋尺寸（长宽高）
BS_pos = [1, 1, 3;      % 基站 1
          1, 9, 3;      % 基站 2
          9, 9, 3];     % 基站 3
num_BS = size(BS_pos, 1);
N_RB = 100;
N_subframes = 5;
f_c = 2.4e9;
% target_pos = [5, 6, 1]; % 待测物体位置
% grid_res = 0.1;
% [est_pos, ~] = FusionLocalization3D(target_pos, BS_pos, roomsize, 'ISSM', grid_res, -10, false);
%[txGrid, f_list, N_subc, M_symbols] = LTE_generate(N_RB, N_subframes, f_c);