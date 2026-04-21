function Rxx_smoothed = spatial_smoothing_2d(X_signal, N_side, sub_N)
% SPATIAL_SMOOTHING_2D 二维前向空间平滑算法 (2D-FSS)Forward-Backward Spatial Smoothing, FBSS
% 输入:
%   X_signal - 原始接收信号矩阵 (尺寸: N_side^2 x M快拍数)
%   N_side   - 原始正方形阵列的边长天线数 (例如 4)
%   sub_N    - 划分的子阵列边长天线数 (例如 3)
% 输出:
%   Rxx_smoothed - 平滑后的协方差矩阵 (尺寸: sub_N^2 x sub_N^2)

    [N_total, M_snapshots] = size(X_signal);
    if N_total ~= N_side^2
        error('输入信号的维度与物理阵列大小不匹配！');
    end
    if sub_N >= N_side
        error('子阵边长必须小于原阵列边长！');
    end

    % 计算X和Y方向上可以划分的子阵数量
    L_x = N_side - sub_N + 1;
    L_y = N_side - sub_N + 1;
    L_total = L_x * L_y; % 总子阵数

    % 子阵的协方差矩阵初始化
    sub_N_total = sub_N^2;
    Rxx_smoothed = zeros(sub_N_total, sub_N_total);

    % 遍历每一个子阵
    for lx = 0 : L_x - 1
        for ly = 0 : L_y - 1
            % 提取当前子阵包含的天线在原始信号中的索引
            sub_array_idx = zeros(1, sub_N_total);
            count = 1;
            
            % 根据你生成位置的规则：外层x，内层y
            for kx = 0 : sub_N - 1
                for ky = 0 : sub_N - 1
                    % 计算全局 x 和 y 坐标索引,即此时的这根天线在天线阵列中的位置
                    global_x = lx + kx;
                    global_y = ly + ky;
                    
                    % 映射回一维索引 (公式对应你生成 pos 时的逻辑)
                    % idx = x_idx * N_side + y_idx + 1
                    linear_idx = global_x * N_side + global_y + 1; 
                    
                    sub_array_idx(count) = linear_idx; %将子阵列的所有天线位置在POS中的索引找出
                    count = count + 1;
                end
            end
            
            % 提取该子阵的接收信号
            X_sub = X_signal(sub_array_idx, :);
            
            % 计算该子阵的协方差矩阵并累加
            R_sub = (X_sub * X_sub') / M_snapshots;
            Rxx_smoothed = Rxx_smoothed + R_sub;
        end
    end

    % 求平均
    Rxx_smoothed = Rxx_smoothed / L_total;
    % --- 增加后向平滑 (Backward Smoothing) ---
    % 构造交换矩阵 J (副对角线为1，其余为0)
    J = fliplr(eye(sub_N_total)); 
    % 计算前向-后向平均协方差矩阵
    Rxx_smoothed = 0.5 * (Rxx_smoothed + J * conj(Rxx_smoothed) * J);
end