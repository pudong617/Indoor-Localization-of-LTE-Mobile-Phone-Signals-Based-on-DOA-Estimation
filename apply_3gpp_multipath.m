function X_rx = apply_3gpp_multipath(txGrid, f_list, pos, theta_true, phi_true, SNR)
% APPLY_3GPP_MULTIPATH 应用 3GPP InH 簇-射线 (Cluster-Ray) 结构信道模型
% 参考标准: 3GPP TR 38.901 Section 7.5 / TR 36.873 Section 7.3
% 输入参数：
%   txGrid     - 纯净的 LTE 频域网格 (N_subc x M_symbols)
%   f_list     - 子载波频率列表
%   pos        - 接收端物理天线坐标矩阵 (N_ant x 3)
%   theta_true - 目标真实方位角 (直射径角度)
%   phi_true   - 目标真实俯仰角 (直射径角度)
%   SNR        - 信噪比 (dB)
% 输出：
%   X_rx       - 三维接收信号张量 (N_ant x M_symbols x N_subc)
    
    c = 3e8;
    [N_subc, M_symbols] = size(txGrid); 
    N_ant = size(pos, 1);
    
    %% 1. 3GPP 室内多径抽头 (CDL Taps / Clusters)
    % 包含 14 个簇 (Clusters) 的中心时延和平均功率
    tau_3gpp = [0, 0, 0.7, 12.24, 27.26, 28.1, 36.08, 51.92, 35.5, 80.84, 158.74, 188.48, 194.16, 250.5] * 1e-9;
    pdb_3gpp = [-0.2, -13.5, -18.8, -21, -22.8, -17.9, -20.1, -21.9, -22.9, -27.8, -23.6, -24.8, -30, -27.7]; 
    amp_linear = 10.^(pdb_3gpp / 20); % 簇的总幅度系数
    N_clusters = length(pdb_3gpp);
    
    %% 2. 空间角度映射 (Macro Angles using Laplacian)
    theta_clusters = zeros(1, N_clusters); 
    phi_clusters   = zeros(1, N_clusters);
    
    % 簇1 (LOS直射径)：对齐真实坐标角度
    theta_clusters(1) = theta_true; 
    phi_clusters(1)   = phi_true;
    
    % 簇2-14：拉普拉斯角度扩展 (假设全局 ASA=25度, ZSA=12度)
    b_theta = 25 / sqrt(2);
    b_phi   = 12 / sqrt(2);
    for p = 2:N_clusters
        u_th = rand() - 0.5;
        u_ph = rand() - 0.5;
        laplacian_th = -b_theta * sign(u_th) * log(1 - 2*abs(u_th));
        laplacian_ph = -b_phi   * sign(u_ph) * log(1 - 2*abs(u_ph));
        
        theta_clusters(p) = theta_true + laplacian_th;
        phi_clusters(p)   = phi_true   + laplacian_ph;
    end
    
    %% 3. 3GPP 簇内射线结构定义 (20 Rays per Cluster)
    % 协议规定的 20 条射线的固定偏移系数 (Table 7.3-3 / 7.5-3)
    ray_offset_bases = [0.0447, 0.1413, 0.2492, 0.3715, 0.5129, 0.6797, 0.8844, 1.1481, 1.5195, 2.1551];
    ray_offsets = zeros(1, 20);
    ray_offsets(1:2:19) = ray_offset_bases;
    ray_offsets(2:2:20) = -ray_offset_bases;
    
    % InH 场景的簇内角度扩展标准差 (Cluster ASA / ZSA) 依据 Table 7.3-6a
    c_ASA = 11; % 方位角簇内扩展 (度)
    c_ZSA = 9;  % 俯仰角簇内扩展 (度)

    %% 4. 宽带信号与信道相乘
    X_rx_clean = zeros(N_ant, M_symbols, N_subc);
    
    for k = 1:N_subc
        f_k = f_list(k);
        lambda = c / f_k; 
        X_k_clean = zeros(N_ant, M_symbols);
        S_k = txGrid(k, :); 
        
        for c_idx = 1:N_clusters
            % 区分直射径簇与反射多径簇
            if c_idx == 1
                % 簇 1: 纯 LOS 信号 (时延0，功率 -0.2dB)
                % 只有 1 条射线，代表纯粹的直射方向，无衰落
                num_rays = 1;
                ray_thetas = theta_clusters(c_idx);
                ray_phis   = phi_clusters(c_idx);
                ray_amps   = amp_linear(c_idx); 
                ray_phases = 1; % 绝对相位参考
                
            else
                % 簇 2 - 14: 散射多径簇 (包含簇2：时延0，功率 -13.5dB)
                % 每个簇包含 20 条射线，引起瑞利衰落和空间相干弥散
                num_rays = 20;
                ray_thetas = theta_clusters(c_idx) + c_ASA * ray_offsets;
                ray_phis   = phi_clusters(c_idx)   + c_ZSA * ray_offsets;
                
                % 簇功率均分给 20 条射线，引入随机相位
                ray_amps = (amp_linear(c_idx) / sqrt(20)) * ones(1, 20);
                ray_phases = exp(1j * 2 * pi * rand(1, 20)); 
            end
            
            % 将当前簇内的所有射线叠加到接收天线上
            for r = 1:num_rays
                th_rad = ray_thetas(r) * pi / 180;
                ph_rad = ray_phis(r) * pi / 180;
                
                kx = cos(th_rad) * cos(ph_rad);
                ky = cos(ph_rad) * sin(th_rad);
                kz = sin(ph_rad);
                k_vec = (2 * pi / lambda) * [kx; ky; kz];
                
                % 阵列流形
                A_p = exp(-1j * (pos * k_vec)); 
                
                % 幅度 * 随机初始相位 * 频域时延相位扭转
                H_f = ray_amps(r) * ray_phases(r) * exp(-1j * 2 * pi * f_k * tau_3gpp(c_idx));
                
                % 累加该射线的接收信号
                X_k_clean = X_k_clean + A_p * (S_k * H_f);
            end
        end
        X_rx_clean(:, :, k) = X_k_clean; 
    end
    
    %% 5. 添加高斯白噪声
    signal_power = mean(abs(X_rx_clean(:)).^2);
    noise_var = signal_power / (10^(SNR/10));
    N_noise = sqrt(noise_var / 2) * (randn(size(X_rx_clean)) + 1j * randn(size(X_rx_clean)));
    
    X_rx = X_rx_clean + N_noise;
end