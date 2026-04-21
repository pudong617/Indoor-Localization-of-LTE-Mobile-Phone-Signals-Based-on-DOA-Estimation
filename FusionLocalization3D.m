function [est_pos, P_fusion_3D] = FusionLocalization3D(target_pos, BS_pos, roomsize, ISSM, grid_res, SNR, plot_flag)
% FUSIONLOCALIZATION3D 多基站联合 3D MUSIC 空间概率热力图融合定位
%
% 输入参数:
%   target_pos - 目标的真实三维坐标 [X, Y, Z]
%   BS_pos     - 所有基站的三维坐标矩阵 (N_BS x 3)
%   roomsize   - 房间尺寸 [长, 宽, 高]
%   grid_res   - 三维体素网格分辨率 (如 0.1 米)
%   SNR        - 仿真信噪比 (dB)
%   plot_flag  - 是否绘图 (布尔值: true 或 false)。蒙特卡洛时务必设为 false!
%
% 输出参数:
%   est_pos     - 融合算法估计出的目标三维坐标 [X, Y, Z]
%   P_fusion_3D - 完整的 3D 融合概率矩阵，用于后期分析

    % 1. 初始化变量
    num_BS = size(BS_pos, 1);
    theta_range = [-180, 180]; 
    P_linear_list = cell(1, num_BS);
    
    % 2. 模拟各个基站获取 3D MUSIC 空间谱
    for i = 1:num_BS
        % 计算真实的相对角度
        [true_theta, true_phi] = jisuanangle(target_pos, BS_pos(i,:));
        
        % 静默调用 DOA 估计函数 
        switch ISSM
            case 'ISSM'
                [P_linear, ~] = DOAtianxian(BS_pos(i,:), true_theta, true_phi, SNR, 6, theta_range,false);
            case 'noISSM'
                [P_linear, ~] = DOAtianxian_no(BS_pos(i,:), true_theta, true_phi, SNR, 6, theta_range,false);
            otherwise
                [P_linear, ~] = DOAtianxian(BS_pos(i,:), true_theta, true_phi, SNR, 6, theta_range,false);
        end
        P_linear_list{i} = P_linear;
    end

    % 3. 生成三维体素网格
    x_vec = 0 : grid_res : roomsize(1);
    y_vec = 0 : grid_res : roomsize(2);
    z_vec = 0 : grid_res : roomsize(3);
    [X_vol, Y_vol, Z_vol] = meshgrid(x_vec, y_vec, z_vec);
    
    % 初始化全局连乘矩阵
    P_fusion_3D = ones(size(X_vol));
    
    % 匹配 DOAtianxian 内部的插值网格向量
    theta_search_vec = linspace(theta_range(1), theta_range(2), 361);
    phi_search_vec   = linspace(-90, 0, 91); % 严格匹配你修改后的 [-90, 0] 范围
    
    % 4. 空间概率热力图融合
    for i = 1:num_BS
        % 计算空间网格点相对于当前基站的角度
        dX = X_vol - BS_pos(i, 1);
        dY = Y_vol - BS_pos(i, 2);
        dZ = Z_vol - BS_pos(i, 3);
        R_xy = sqrt(dX.^2 + dY.^2 + 1e-10);
        
        Theta_vol = atan2(dY, dX) * 180 / pi;
        Phi_vol   = atan2(dZ, R_xy) * 180 / pi;
        
        % 插值映射
        P_mapped = interp2(phi_search_vec, theta_search_vec, P_linear_list{i}, ...
                           Phi_vol, Theta_vol, 'linear', 0);
                       
        % 乘法融合
        P_fusion_3D = P_fusion_3D .* P_mapped;
    end
    
    % 5. 寻找极值点输出坐标
    [~, max_idx] = max(P_fusion_3D(:));
    [iy, ix, iz] = ind2sub(size(P_fusion_3D), max_idx);
    est_pos = [x_vec(ix), y_vec(iy), z_vec(iz)];
    
    % 6. (可选) 绘图逻辑
    if nargin < 7
        plot_flag = false; % 默认不绘图
    end
    
    if plot_flag
        fprintf('>>> 融合定位估计坐标: (%.2f, %.2f, %.2f) <<<\n', est_pos(1), est_pos(2), est_pos(3));
        fprintf('>>> 目标真实物理坐标: (%.2f, %.2f, %.2f) <<<\n\n', target_pos(1), target_pos(2), target_pos(3));
        
        figure('Color', 'w', 'Name', '3D Multi-Array Fusion');
        hold on; grid on; view(3);
        %slice(X_vol, Y_vol, Z_vol, P_fusion_3D, est_pos(1), est_pos(2), est_pos(3));
        %shading interp; colormap('jet'); colorbar;
        
        plot3(BS_pos(:,1), BS_pos(:,2), BS_pos(:,3), 'b^', 'MarkerSize', 12, 'LineWidth', 2,'DisplayName', '基站');
        plot3(target_pos(1), target_pos(2), target_pos(3), 'ko', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', '真实目标');
        plot3(est_pos(1), est_pos(2), est_pos(3), 'r*', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', '估计目标');
        
        num_BS = size(BS_pos, 1); % 获取基站数量
        for i = 1:num_BS
            % plot3 需要 x, y, z 分开的向量：[起点坐标, 终点坐标]
            plot3([BS_pos(i, 1), est_pos(1)], ... % X 坐标向量
                  [BS_pos(i, 2), est_pos(2)], ... % Y 坐标向量
                  [BS_pos(i, 3), est_pos(3)], ... % Z 坐标向量
                  'r--', 'LineWidth', 1, 'HandleVisibility', 'off'); % 绘制黑色虚线，不显示在图例中
        end
        % ------------------------------------------------
        
        % 添加图例
        legend('show', 'Location', 'northeast');

        axis equal; axis([0 roomsize(1) 0 roomsize(2) 0 roomsize(3)]);
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title('三维空间多基站概率相乘融合定位');
    end
end