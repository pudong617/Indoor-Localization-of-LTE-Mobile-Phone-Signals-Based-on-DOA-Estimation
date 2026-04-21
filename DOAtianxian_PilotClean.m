function [P_music_global, P_music_db] = DOAtianxian_PilotClean(pos_first, theta_true, phi_true, SNR, N_side, thetasize, plot_flag)
% DOATIANXIAN_PILOTCLEAN 基于时域导频提取与多径剔除的 2D-MUSIC 算法
% 严格遵循物理层接收逻辑，利用 DMRS 导频估计 CIR 并加窗隔离多径

%% 1. 系统与阵列参数初始化
c = 3e8;                   
f_c = 2.4e9;               
lambda_c = c / f_c;        
d = lambda_c / 2;          

% 生成正方形接收阵列坐标 (使用全阵列)
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

%% 2. 触发基带信号生成与信道衰落
N_RB = 100;        
N_subframes = 5;  
[txGrid, f_list, N_subc, M_symbols] = LTE_generate(N_RB, N_subframes, f_c);
X_rx = apply_3gpp_multipath(txGrid, f_list, pos, theta_true, phi_true, SNR);

%% 3. 【核心】提取导频并计算纯净信道频域响应 H(f)
% PUSCH 的 DMRS 导频固定位于每个时隙的第 4 个符号
dmrs_symbols = [4, 11, 18, 25, 32, 39, 46, 53, 60, 67];
N_dmrs = length(dmrs_symbols);

% 提取接收端的导频符号矩阵 [N_ant, N_dmrs, N_subc]
Y_dmrs = X_rx(:, dmrs_symbols, :); 
% 提取本地已知的导频复本 [N_subc, N_dmrs]
X_pilot_local = txGrid(:, dmrs_symbols); 

H_est_raw = zeros(N_total, N_dmrs, N_subc);

% 严格执行信道估计 (消除数据调制)
for i_sym = 1:N_dmrs
    % 取出当前导频符号的所有天线接收数据: [N_ant, N_subc]
    Y_sym = squeeze(Y_dmrs(:, i_sym, :)); 
    % 取出本地复本对应的列，转置为行向量以便广播: [1, N_subc]
    X_pilot_sym = X_pilot_local(:, i_sym).'; 
    % 矩阵除法: 消除相位旋转，得到 H
    H_est_raw(:, i_sym, :) = Y_sym ./ X_pilot_sym;
end

% 在 10 个导频符号上求平均，有效压制高斯白噪声
H_mean_f = squeeze(mean(H_est_raw, 2)); % 结果尺寸: [N_total, N_subc]

%% 4. 【核心】时域 IFFT 加窗剔除多径 (Time-Domain Gating)
N_ifft = 8192; 
H_clean_f = zeros(N_total, N_subc);

% 为了绘图对比，保留第一根天线的原始和清理后的时域波形
h_t_raw_plot = zeros(1, N_ifft); 
h_t_clean_plot = zeros(1, N_ifft);

for i_ant = 1:N_total
    % 频域补零 (Zero-padding)
    H_ant = H_mean_f(i_ant, :);
    H_padded = zeros(1, N_ifft);
    start_idx = floor((N_ifft - N_subc)/2) + 1;
    H_padded(start_idx : start_idx + N_subc - 1) = H_ant;
    
    % IFFT 变至时域 (高分辨率)
    h_t = fftshift(ifft(ifftshift(H_padded)));
    
    % --- 时域加窗逻辑 ---
    % 1. 寻找直射径（室内视距下通常是最强峰）
    [~, peak_idx] = max(abs(h_t)); 
    
    % 2. 构建矩形窗：保留主峰前后的安全裕度，剔除远端多径拖尾
    window_radius = 20; % 容差窗口宽度（可根据带宽和采样率微调）
    h_t_clean = zeros(1, N_ifft);
    
    % 确保窗口索引不越界
    idx_left = max(1, peak_idx - window_radius);
    idx_right = min(N_ifft, peak_idx + window_radius);
    
    % 仅保留直射径窗口内的能量，多径被强制置零
    h_t_clean(idx_left : idx_right) = h_t(idx_left : idx_right);
    
    if i_ant == 1
        h_t_raw_plot = h_t;
        h_t_clean_plot = h_t_clean;
    end
    
    % 3. FFT 变回频域，得到没有多径干扰的纯净 H(f)
    H_padded_clean = fftshift(fft(ifftshift(h_t_clean)));
    H_clean_f(i_ant, :) = H_padded_clean(start_idx : start_idx + N_subc - 1);
end

%% 5. 全阵列 2D MUSIC 估计 (无需空间平滑)
% 抽取子载波以减少计算量
N_music_subc = 32; 
step = floor(N_subc / N_music_subc);
music_indices = 1:step:N_subc;

theta_search = linspace(thetasize(1), thetasize(2), 361);  
phi_search = linspace(-90, 0, 91);       
P_music_global = zeros(length(theta_search), length(phi_search));

% 向量化预计算搜索网格
[Theta_grid, Phi_grid] = ndgrid(theta_search, phi_search);
th_rad = Theta_grid(:).' * pi / 180;
ph_rad = Phi_grid(:).' * pi / 180;
K_mat = [cos(th_rad) .* cos(ph_rad); 
         sin(th_rad) .* cos(ph_rad); 
         sin(ph_rad)];

% 因为多径已经被时域彻底切除，信道中只剩下 1 条直射径
D_music = 1; 

% 遍历选定的子载波频点计算伪谱
for idx_k = 1:length(music_indices)
    k = music_indices(idx_k);
    f_k = f_list(k);
    lambda_f = c / f_k; 
    
    % 取出该频点下的纯净天线响应 (列向量)
    H_k = H_clean_f(:, k);
    
    % 构建协方差矩阵 (外积)
    Rxx_k = H_k * H_k';
    
    % 特征分解
    [V, Lambda] = eig(Rxx_k);            
    [~, idx_eig] = sort(diag(Lambda), 'descend');  
    V = V(:, idx_eig);                     
    
    % 提取噪声子空间 (从 D_music+1 到 N_total)
    Vn = V(:, D_music+1:N_total); 
    
    % 计算导向矢量
    A_search = exp(-1j * (2 * pi / lambda_f) * (pos * K_mat));
    
    % 空间谱投影计算
    Proj = Vn' * A_search; 
    denominator = sum(abs(Proj).^2, 1); 
    
    P_music_k_vec = 1 ./ max(denominator, 1e-10);
    P_music_k = reshape(P_music_k_vec, length(theta_search), length(phi_search));
    
    P_music_global = P_music_global + P_music_k;
end

P_music_db = 10 * log10(P_music_global / max(P_music_global(:)));

%% 6. 绘图：验证多径剔除效果
if plot_flag
    [peak_thetas, peak_phis, peak_vals_db] = find_2d_peaks(P_music_global, theta_search, phi_search, -15, 1);
    fprintf('\n共找到 %d 个有效的 MUSIC 局部峰值。\n', length(peak_thetas));
    fprintf('真实 LOS 预设角度 : (%.1f, %.1f)\n', theta_true, phi_true);
    
    % --- 绘制时域加窗对比图 ---
    figure('Color', 'w', 'Name', 'Time-Domain Multipath Gating');
    delta_f = 15e3; 
    dt = 1 / (N_ifft * delta_f); 
    t_axis = (-N_ifft/2 : N_ifft/2 - 1) * dt * 1e9; 
    
    plot(t_axis, abs(h_t_raw_plot), 'k-', 'LineWidth', 1.5, 'DisplayName', '原始 CIR (包含多径)');
    hold on; grid on;
    plot(t_axis, abs(h_t_clean_plot), 'r-', 'LineWidth', 2, 'DisplayName', '加窗后的 CIR (仅剩直射径)');
    
    xlim([-100, 300]); 
    xlabel('Time Delay (ns)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude |h(t)|', 'FontSize', 12, 'FontWeight', 'bold');
    title('时域多径切除对比', 'FontSize', 14);
    legend('Location', 'northeast');

    % --- 绘制纯净 MUSIC 谱 ---
    figure('Color', 'w', 'Name', 'Clean 2D MUSIC Spectrum');
    [Phi_grid_plot, Theta_grid_plot] = meshgrid(phi_search, theta_search);
    surf(Theta_grid_plot, Phi_grid_plot, P_music_db, 'EdgeColor', 'none');
    shading interp; colormap jet; colorbar;
    hold on;
    plot3(theta_true, phi_true, max(P_music_db(:)), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '真实坐标');
    xlabel('方位角 \theta'); ylabel('俯仰角 \phi'); zlabel('伪谱幅度 (dB)');
    title('剔除多径后的全阵列 MUSIC 空间谱');
end
end