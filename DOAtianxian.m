function [P_music_global,P_music_db] = DOAtianxian(pos_first, theta_true, phi_true, SNR, N_side, thetasize, plot_flag)
% DOATIANXIAN 适配 3GPP LTE 标准宽带信号与 EPA 信道的 2D MUSIC 估计算法

%% 1. 系统与阵列参数初始化
c = 3e8;                   
f_c = 2.4e9;               % 中心频率 2.4GHz
lambda_c = c / f_c;        
d = lambda_c / 2;          % 物理天线间距按中心频点设定

% 生成正方形接收阵列坐标
N_total = N_side^2;
pos = zeros(N_total, 3);
idx = 1;
for x_idx = 0:N_side-1
    for y_idx = 0:N_side-1
        pos(idx, 1) = pos_first(1) + x_idx * d;  %所有天线的x轴坐标
        pos(idx, 2) = pos_first(2) + y_idx * d;  %所有天线的y轴坐标
        pos(idx, 3) = pos_first(3);              %所有天线的z轴坐标
        idx = idx + 1;
    end
end

% 空间平滑所需的子阵列坐标生成
if N_side < 4
    error('天线阵列太小，N_side 至少需要为 4 才能解相干 7 条多径！');
else
    sub_N = N_side - 2;  %需要所有子阵列数目大于7 Lx = Nx - Mx + 1 >=3; Mx <= Nx-2
end
pos_sub = zeros(sub_N^2, 3);
idx_sub = 1;
for x_idx = 0:sub_N-1
    for y_idx = 0:sub_N-1
        pos_sub(idx_sub, 1) =  x_idx * d; %相对坐标
        pos_sub(idx_sub, 2) =  y_idx * d;
        pos_sub(idx_sub, 3) = 0;
        idx_sub = idx_sub + 1;
    end
end

%% 2. 触发基带信号生成与信道衰落 
N_RB = 100;        % 1000 个资源块 (对应约 3MHz 带宽，在仿真和性能间取良好平衡)
N_subframes = 5;  % 仿真 5 个子帧 (5ms, 提供充足快拍数)

% (a) 获取纯净的 LTE 上行资源网格
[txGrid, f_list, N_subc, ~] = LTE_generate(N_RB, N_subframes, f_c);

% (b) 过 3GPP EPA 信道，投射到接收天线阵列
X_rx = apply_3gpp_multipath(txGrid, f_list, pos, theta_true, phi_true, SNR);

%% -----------------------------------------------------------
%% 3. 抽取导频子载波进行 ISSM (非相干子空间平滑) - 【极致向量化版】
N_music_subc = 32; % 抽取 32 个子载波
step = floor(N_subc / N_music_subc);
music_indices = 1:step:N_subc;

theta_search = linspace(thetasize(1), thetasize(2), 361);  
phi_search = linspace(-90, 0, 91);       
P_music_global = zeros(length(theta_search), length(phi_search));
D_music = 14; 

% ==================== 向量化核心：预处理搜索网格 ====================
% 使用 ndgrid 铺平所有搜索角度，避免在循环内部重复计算三角函数
[Theta_grid_search, Phi_grid_search] = ndgrid(theta_search, phi_search);
Theta_vec = Theta_grid_search(:).'; % 展平为 1 x N_angles
Phi_vec = Phi_grid_search(:).';     % 展平为 1 x N_angles

th_rad = Theta_vec * pi / 180;
ph_rad = Phi_vec * pi / 180;

% 计算统一的空间波矢方向余弦矩阵 K_mat (尺寸: 3 x N_angles)
K_mat = [cos(th_rad) .* cos(ph_rad); 
         sin(th_rad) .* cos(ph_rad); 
         sin(ph_rad)];
% =================================================================

% 遍历抽取出的子载波
for idx_k = 1:length(music_indices)
    k = music_indices(idx_k);
    f_k = f_list(k);
    X_k = X_rx(:, :, k); 
    lambda_f = c / f_k; 
    
    % 二维空间平滑
    Rxx_k = spatial_smoothing_2d(X_k, N_side, sub_N);
    
    % 特征分解
    [V, Lambda] = eig(Rxx_k);            
    [~, idx_eig] = sort(diag(Lambda), 'descend');  
    V = V(:, idx_eig);                     
    Vn = V(:, D_music+1:sub_N^2); % 提取噪声子空间 (尺寸: sub_N^2 x 噪声维数)
    
    % ------------------- 极速向量化扫描 -------------------
    % 1. 一次性算出当前频点下，所有搜索角度的导向矢量矩阵
    % 尺寸运算: (sub_N^2 x 3) * (3 x N_angles) = sub_N^2 x N_angles
    A_search = exp(-1j * (2 * pi / lambda_f) * (pos_sub * K_mat));
    
    % 2. 一次性计算所有角度在噪声子空间上的投影
    % 尺寸运算: Vn' * A_search = 噪声维数 x N_angles
    Proj = Vn' * A_search; 
    
    % 3. 按列求绝对值平方和 ( a'*(Vn*Vn')*a 替代写法)
    denominator = sum(abs(Proj).^2, 1); 
    
    % 4. 求倒数得到该频点的伪谱，并塑形回二维矩阵
    P_music_k_vec = 1 ./ max(denominator, 1e-10);
    P_music_k = reshape(P_music_k_vec, length(theta_search), length(phi_search));
    % ------------------------------------------------------
    
    % 宽带谱累加
    P_music_global = P_music_global + P_music_k;
end
% 归一化伪谱 (dB) 以增强对比度
P_music_db = 10 * log10(P_music_global / max(P_music_global(:)));

%% 绘图、结果分析功能
if plot_flag
    % 设定算法参数
    threshold_dB = -15;    % 只保留距离最高峰 15dB 以内的峰 (过滤掉那些 -20dB 的背景底噪突起)
    search_radius = 1;     % 搜索半径(网格点数)。如果你的网格步长是1度，半径5就代表寻找 +/-5 度范围内的局部最高点
    
    % 调用找峰算法
    [peak_thetas, peak_phis, peak_vals_db] = find_2d_peaks(P_music_global, theta_search, phi_search, threshold_dB, search_radius);
    
    fprintf('\n共找到 %d 个有效的 MUSIC 局部峰值。\n', length(peak_thetas));
    for p = 1:length(peak_thetas)
        fprintf('峰 %d: 方位角 %.1f, 俯仰角 %.1f (相对谱高: %.2f dB)\n', ...
                p, peak_thetas(p), peak_phis(p), peak_vals_db(p));
    end
    fprintf('真实 LOS 预设角度 : (%.1f, %.1f)\n', theta_true, phi_true);
    
    
    
    %% 5. 可视化 MUSIC 空间谱
    figure('Color', 'w', 'Name', '2D MUSIC Spectrum');
    % 使用 meshgrid 生成用于绘图的坐标矩阵
    [Phi_grid, Theta_grid] = meshgrid(phi_search, theta_search);
    
    
    
    % 绘制三维曲面
    surf(Theta_grid, Phi_grid, P_music_db, 'EdgeColor', 'none');
    shading interp; % 平滑色彩
    
    % 标注估计值
    hold on;
    %plot3(theta_true,phi_true,0,'bo','MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'true pos');
    
    %% --- 提取信道冲激响应 (CIR) 并在时域画出包络 ---
    % 1. 频域信道估计 (与原代码相同)
    rx_ant1 = squeeze(X_rx(1, :, :)).'; 
    H_est = rx_ant1 ./ txGrid; 
    H_est_mean = mean(H_est, 2); % 此时长度为 N_subc，频率排列为 [-B/2, ..., +B/2]
    
    % 2. 严谨的时域超分辨补零 (Zero-padding)
    N_ifft = 8192; 
    H_padded = zeros(N_ifft, 1);
    
    % 将原本的有效频域数据放入 8192 点的物理中心
    % 这样两端（高频部分）自动补零，不破坏基带频谱形状
    start_idx = floor((N_ifft - N_subc)/2) + 1;
    H_padded(start_idx : start_idx + N_subc - 1) = H_est_mean;
    
    % 3. IFFT 处理核心逻辑
    % (a) ifftshift: 将中心 0Hz 移动到数组的第一个元素，迎合 MATLAB ifft 的死板要求
    H_ready = ifftshift(H_padded);
    
    % (b) 进行 IFFT，得到高分辨率时域信号
    h_t_raw = ifft(H_ready);
    
    % (c) fftshift: 将时域信号从 [0~T] 重新居中为 [-T/2 ~ T/2] 
    % 这样既能看到正时延（多径），也能看到可能由于负向时间偏移造成的波形
    h_t = fftshift(h_t_raw);
    
    % 4. 计算绝对真实的时间轴 (以 0 为绝对对称中心)
    delta_f = 15e3; % LTE 子载波间隔 15kHz
    dt = 1 / (N_ifft * delta_f); % IFFT 后的时域采样间隔
    
    % 生成以 0 为中心的时间轴 (ns)
    t_axis = (-N_ifft/2 : N_ifft/2 - 1) * dt * 1e9; 
    
    % 5. 可视化 CIR 包络
    figure('Color', 'w', 'Name', 'Time-Domain Channel Impulse Response (CIR)');
    plot(t_axis, abs(h_t), 'b-', 'LineWidth', 2, 'DisplayName', 'Estimated CIR Envelope');
    hold on; grid on;
    
    % % 标注出 3GPP InH 信道的理论真实时延位置作为红线参考
    % tau_3gpp_ns = [0, 10, 20, 30, 45, 55, 70];
    % for i = 1:length(tau_3gpp_ns)
    %     if i == 1
    %         xline(tau_3gpp_ns(i), 'r--', 'LineWidth', 1.5, 'DisplayName', 'True EPA Delays');
    %     else
    %         xline(tau_3gpp_ns(i), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    %     end
    % end
    
    % 设置视图范围，为了看清 0 ns 附近的对齐情况，我们往负轴稍微看一点 (-50 ns)
    xlim([-450, 450]); 
    xlabel('Time Delay (ns)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude |h(t)|', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('CIR Envelope (Bandwidth: %.1f MHz)', (N_subc * delta_f)/1e6), 'FontSize', 14);
    legend('Location', 'northeast');
end

end
