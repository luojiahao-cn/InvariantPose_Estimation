%--------------------------------------------------------------------------
% 多磁铁实验：分析不同磁铁数量/布局下的定位能力
% 参考 expA_main.m 的组织方式，比较 4 种算法在 4 类工况下的误差
% 维度： [配置(磁铁数量/布局), 工况, 算法, 样本]
% 输出：results/expC_main.mat
%--------------------------------------------------------------------------

clc; clear; close all;
addpath('utils')
addpath('Functions')

rng(2025);

%% 参数设置
num_points = 10;            % 每个配置的采样数量
offset = 0.02;              % 初值偏离度（米），固定一个代表性偏差
labels = {'无干扰噪声','仅干扰','仅噪声','干扰和噪声'};
alg_labels = {'LM','ELM','Rlm','Ours'};
alg_colors = lines(numel(alg_labels)); 

% 传感器相对主体的安装阵列（与 expA 保持一致，零均值）
d_list_e = [
    [0; 0; 0], ...
    [1e-3; 0; 0], ...
    [-1e-3; 0; 0], ...
    [0; 1e-3; 0], ...
    [0; -1e-3; 0], ...
    [0; 0; 1e-3], ...
    [0; 0; -1e-3]
    ];
row_means = mean(d_list_e, 2);
d_list = d_list_e - row_means;
num_sensors = size(d_list, 2);

% LM 选项
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off');

% 工作空间与边界
workspace_center = [0; 0; 0];
workspace_radius = 0.5;
lb_p = workspace_center - workspace_radius;
ub_p = workspace_center + workspace_radius;

% 噪声与干扰
noise_level = 3e-5;               % 高斯白噪声标准差（T）
disturbance = [0; 0; 2e-3];       % 常偏干扰（T）


% 固定磁铁参数
m_pos = [
    [-0.3; 0; 0], ...
    [0.2; 0; 0], ...
    [0; 0.25; 0], ...
    [0; -0.25; 0]
    ];
m_hat = [
    [1; 0; 0], ...
    [0; 0; 1], ...
    [0; 1; 0], ...
    [0; 0; -1]
    ];
m_hat = m_hat ./ vecnorm(m_hat);
m_norm = [8e2, 8e2, 8e2, 8e2];
num_magnets = size(m_pos, 2);


num_alg = numel(alg_labels);
num_cond = numel(labels);

% 误差存储 [条件, 算法, 样本]
err_pos = zeros(num_cond, num_alg, num_points);
err_rot = zeros(num_cond, num_alg, num_points);

%% 主循环


hwb = waitbar(0, '正在运行多磁铁实验，请稍候...');
total_samples = num_points;

for pt_idx = 1:num_points
    % 进度提示
    percent_sample = pt_idx / num_points;
    elapsed = toc;
    if pt_idx > 1
        avg_time = elapsed / pt_idx;
        est_left = (num_points - pt_idx) * avg_time;
    else
        est_left = 0;
    end
    waitbar(percent_sample, hwb, ...
        sprintf('进度 %.2f%% (%d/%d)\n已用时 %.1fs，预计剩余 %.1fs', ...
            percent_sample*100, pt_idx, num_points, elapsed, est_left));

    % 随机生成真值姿态与位置
    theta_true = -pi + 2*pi*rand(3,1);
    p_true = -1e-1 + 2e-1*rand(3,1);
    R_true = MatrixExp3(VecToso3(theta_true));
    sensor_positions = p_true + R_true * d_list;

    % 生成磁场（世界坐标），叠加所有磁铁
    B_total = zeros(3, num_sensors);
    for sensor_idx = 1:num_sensors
        B_t = zeros(3,1);
        for magnet_idx = 1:num_magnets
            r = sensor_positions(:, sensor_idx) - m_pos(:, magnet_idx);
            moment_unit = m_hat(:, magnet_idx);
            moment_mag = m_norm(magnet_idx);
            [B_single, ~] = dipole_b_and_gradb(r, moment_unit, moment_mag);
            B_t = B_t + B_single;
        end
        B_total(:, sensor_idx) = B_t;
    end
    % 转到机体系（传感器测到的 b）
    b_total = R_true' * B_total;

    % 初值（固定偏离度 offset）
    init_error = -1 + 2*rand(3,1);
    p_init = p_true + offset * init_error;
    theta_init = theta_true + offset * 10 * (-pi + 2*pi*init_error);
    R_init = MatrixExp3(VecToso3(theta_init));

    % 四种工况
    for mode = 1:num_cond
        b_meas = b_total;
        switch mode
            case 1 % 无干扰噪声
                % 保持 b_total
            case 2 % 仅干扰
                b_meas = b_meas + disturbance;
            case 3 % 仅噪声
                b_meas = b_meas + randn(size(b_meas)) * noise_level;
            case 4 % 干扰和噪声
                b_meas = b_meas + disturbance + randn(size(b_meas)) * noise_level;
        end

        % 算法比较
        % 1) LM
        [p_lm, R_lm, ~] = estimate_pose_lm( ...
            b_meas, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
        err_pos(mode, 1, pt_idx) = norm(p_lm - p_true);
        err_rot(mode, 1, pt_idx) = norm(R_lm - R_true, 'fro');

        % 2) ELM
        [p_elm, R_elm, ~] = estimate_pose_elm( ...
            b_meas, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
        err_pos(mode, 2, pt_idx) = norm(p_elm - p_true);
        err_rot(mode, 2, pt_idx) = norm(R_elm - R_true, 'fro');

        % 3) Ours（注意该函数接口需要真值用于中间步骤/评估）
        [p_ours, R_ours, ~] = estimate_pose_ours( ...
            b_meas, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p );
        err_pos(mode, 4, pt_idx) = norm(p_ours - p_true);
        err_rot(mode, 4, pt_idx) = norm(R_ours - R_true, 'fro');

        % 4) Rlm（先用 Ours 的位置，再估计旋转）
        [p_Rlm, R_Rlm, ~] = estimate_R_lm( ...
            b_meas, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );
        err_pos(mode, 3, pt_idx) = norm(p_Rlm - p_true);
        err_rot(mode, 3, pt_idx) = norm(R_Rlm - R_true, 'fro');
    end
end
close(hwb);
fprintf('全部多磁铁实验完成，总耗时 %.2fs\n', toc);

% 保存结果
save('../results/expC_main.mat', 'err_pos', 'err_rot', 'num_cond', 'num_alg', ...
    'num_points', 'labels', 'alg_labels', 'm_pos', 'm_hat', 'm_norm');

% 提示：可根据需要编写单独绘图脚本对 err_pos/err_rot 在不同 Nm、不同工况下做箱线图/均值±方差带等可视化。
