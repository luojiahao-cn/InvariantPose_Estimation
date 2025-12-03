% 收敛域测绘实验（Basin of Attraction）
% 对比：LM/ELM/Rlm/Ours
% 输出：热力图（半径→成功率）、耗时箱线图

clc; clear; close all;
addpath('utils')
addpath('Functions')

rng(2025);

%% 参数设置（可修改）
num_pose_samples = 10;      % 真值位姿数量
num_init_per_pose = 50;     % 每个真值的初值数量
radius_list = linspace(1e-3, 50e-3, 10); % 初值扰动半径（米）
angle_list = linspace(1, 30, 6);         % 初值扰动角度（度）

success_thresh_pos = 3e-3;   % 成功判定阈值（位置3mm）
success_thresh_rot = deg2rad(5); % 成功判定阈值（姿态5°）

algorithms = {'LM', 'ELM', 'Rlm', 'Ours'};
alg_colors = {'b', 'g', 'm', 'r'};

results = struct();

%% 磁铁参数定义
m_pos = [
    [-0.3; 0; 0], ...
    [0.2; 0; 0]
    ];
m_hat = [
    [1; 0; 0], ...
    [0; 0; 1]
    ];
m_hat = m_hat ./ vecnorm(m_hat);
m_norm = [8e2, 8e2];

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
num_magnets = size(m_pos, 2);

% 优化器参数
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'off', ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6, ...
    'MaxIter', 1000, ...
    'MaxFunctionEvaluations', 10000);

workspace_center = [0; 0; 0];
workspace_radius = 1;
lb_p = workspace_center - workspace_radius;
ub_p = workspace_center + workspace_radius;

%% 主循环
% 预分配结果数组（拆分为独立变量）
lm      = false(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
elm     = false(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
Rlm     = false(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
ours    = false(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
lm_time   = zeros(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
elm_time  = zeros(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
Rlm_time  = zeros(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);
ours_time = zeros(length(radius_list), length(angle_list), num_pose_samples, num_init_per_pose);

total_iter = num_pose_samples * length(radius_list) * length(angle_list) * num_init_per_pose;
iter_count = 0;
tic; % 总计时

for pose_idx = 1:num_pose_samples
    pose_tic = tic;
    % 随机生成真值
    theta_true = pi*rand(3,1);
    p_true = -0.05 + 0.1*rand(3,1);
    R_true = MatrixExp3(VecToso3(theta_true));

    % 计算传感器全局位置
    sensor_positions = p_true + R_true * d_list;

    % 生成磁场测量数据
    B_total = zeros(3, num_sensors);
    b_total = zeros(3, num_sensors);
    for sensor_idx = 1:num_sensors
        B_t = zeros(3, 1);
        for magnet_idx = 1:num_magnets
            r = sensor_positions(:, sensor_idx) - m_pos(:, magnet_idx);
            moment_unit = m_hat(:, magnet_idx);
            moment_mag = m_norm(magnet_idx);
            [B_single, ~] = dipole_b_and_gradb(r, moment_unit, moment_mag);
            B_t = B_t + B_single;
        end
        B_total(:, sensor_idx) = B_t;
        b_total(:, sensor_idx) = R_true' * B_t;
    end

    % 添加噪声和地磁干扰
    % B_total = B_total + randn(size(B_total)) * 1e-6;
    % b_total = b_total + randn(size(b_total)) * 1e-6;
    % B_total = B_total + [0; 0; 1e-6];
    % b_total = b_total + [0; 0; 1e-6];

    for r_idx = 1:length(radius_list)
        radius = radius_list(r_idx);
        for a_idx = 1:length(angle_list)
            angle_deg = angle_list(a_idx);
            angle_rad = deg2rad(angle_deg);

            parfor init_idx = 1:num_init_per_pose
                % 随机扰动初值
                p_init = p_true + radius * randn(3,1);
                theta_init = theta_true + angle_rad * randn(3,1);
                R_init = MatrixExp3(VecToso3(theta_init));

                % LM
                t_lm_start = tic;
                [p_lm, R_lm_, ~] = estimate_pose_lm( ...
                    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
                t_lm = toc(t_lm_start);

                % ELM
                t_elm_start = tic;
                [p_elm, R_elm_, ~] = estimate_pose_elm( ...
                    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
                t_elm = toc(t_elm_start);

                % Rlm
                t_Rlm_start = tic;
                [p_Rlm_, R_Rlm_, ~] = estimate_R_lm( ...
                    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options );
                t_Rlm = toc(t_Rlm_start);

                % Ours
                t_ours_start = tic;
                [p_ours_, R_ours_, ~] = estimate_pose_ours( ...
                    b_total, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p );
                t_ours = toc(t_ours_start);

                % 误差
                pos_err_lm   = norm(p_lm   - p_true);
                rot_err_lm   = norm(R_lm_  - R_true, 'fro');
                pos_err_elm  = norm(p_elm  - p_true);
                rot_err_elm  = norm(R_elm_ - R_true, 'fro');
                pos_err_Rlm  = norm(p_Rlm_ - p_true);
                rot_err_Rlm  = norm(R_Rlm_ - R_true, 'fro');
                pos_err_ours = norm(p_ours_ - p_true);
                rot_err_ours = norm(R_ours_ - R_true, 'fro');

                % 成功判定
                succ_lm   = (pos_err_lm   < success_thresh_pos) && (rot_err_lm   < success_thresh_rot);
                succ_elm  = (pos_err_elm  < success_thresh_pos) && (rot_err_elm  < success_thresh_rot);
                succ_Rlm  = (pos_err_Rlm  < success_thresh_pos) && (rot_err_Rlm  < success_thresh_rot);
                succ_ours = (pos_err_ours < success_thresh_pos) && (rot_err_ours < success_thresh_rot);

                % 统计
                lm(r_idx, a_idx, pose_idx, init_idx)   = succ_lm;
                elm(r_idx, a_idx, pose_idx, init_idx)  = succ_elm;
                Rlm(r_idx, a_idx, pose_idx, init_idx)  = succ_Rlm;
                ours(r_idx, a_idx, pose_idx, init_idx) = succ_ours;

                lm_time(r_idx, a_idx, pose_idx, init_idx)   = t_lm;
                elm_time(r_idx, a_idx, pose_idx, init_idx)  = t_elm;
                Rlm_time(r_idx, a_idx, pose_idx, init_idx)  = t_Rlm;
                ours_time(r_idx, a_idx, pose_idx, init_idx) = t_ours;
            end
        end
    end
    pose_time = toc(pose_tic);
    est_left = (num_pose_samples-pose_idx)*pose_time;
    fprintf('Pose %d/%d 完成，用时 %.2fs，预计剩余 %.2fs\n', pose_idx, num_pose_samples, pose_time, est_left);
end

fprintf('全部实验完成，总耗时 %.2fs\n', toc);

% parfor外赋值回results结构体
results.lm      = lm;
results.elm     = elm;
results.Rlm     = Rlm;
results.ours    = ours;
results.lm_time   = lm_time;
results.elm_time  = elm_time;
results.Rlm_time  = Rlm_time;
results.ours_time = ours_time;

%% 统计与可视化
% 成功率统计
lm_success_rate   = squeeze(mean(mean(results.lm,4),3));   % [radius, angle]
elm_success_rate  = squeeze(mean(mean(results.elm,4),3));
Rlm_success_rate  = squeeze(mean(mean(results.Rlm,4),3));
ours_success_rate = squeeze(mean(mean(results.ours,4),3));

% 中位耗时统计
lm_median_time   = squeeze(median(median(results.lm_time,4),3));
elm_median_time  = squeeze(median(median(results.elm_time,4),3));
Rlm_median_time  = squeeze(median(median(results.Rlm_time,4),3));
ours_median_time = squeeze(median(median(results.ours_time,4),3));

figure('Name','收敛域热力图','Color','white');
subplot(2,2,1);
imagesc(radius_list*1e3, angle_list, lm_success_rate'*100);
colorbar; xlabel('初值扰动半径 (mm)'); ylabel('初值扰动角度 (deg)');
title('LM成功率 (%)'); set(gca,'YDir','normal');

subplot(2,2,2);
imagesc(radius_list*1e3, angle_list, elm_success_rate'*100);
colorbar; xlabel('初值扰动半径 (mm)'); ylabel('初值扰动角度 (deg)');
title('ELM成功率 (%)'); set(gca,'YDir','normal');

subplot(2,2,3);
imagesc(radius_list*1e3, angle_list, Rlm_success_rate'*100);
colorbar; xlabel('初值扰动半径 (mm)'); ylabel('初值扰动角度 (deg)');
title('Rlm成功率 (%)'); set(gca,'YDir','normal');

subplot(2,2,4);
imagesc(radius_list*1e3, angle_list, ours_success_rate'*100);
colorbar; xlabel('初值扰动半径 (mm)'); ylabel('初值扰动角度 (deg)');
title('Ours成功率 (%)'); set(gca,'YDir','normal');

figure('Name','耗时分布','Color','white');
boxplot([lm_median_time(:), elm_median_time(:), Rlm_median_time(:), ours_median_time(:)], ...
    'Labels', {'LM','ELM','Rlm','Ours'});
ylabel('中位耗时 (秒)');
title('不同算法中位耗时分布');

% 统计所有扰动下的误差
num_r = length(radius_list);
num_a = length(angle_list);

% 预分配误差数组
lm_pos_errs   = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);
elm_pos_errs  = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);
Rlm_pos_errs  = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);
ours_pos_errs = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);

lm_rot_errs   = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);
elm_rot_errs  = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);
Rlm_rot_errs  = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);
ours_rot_errs = zeros(num_r, num_a, num_pose_samples, num_init_per_pose);

% 重新计算误差（只统计，不判定成功）
for pose_idx = 1:num_pose_samples
    theta_true = pi*rand(3,1);
    p_true = -0.05 + 0.1*rand(3,1);
    R_true = MatrixExp3(VecToso3(theta_true));
    for r_idx = 1:num_r
        for a_idx = 1:num_a
            for init_idx = 1:num_init_per_pose
                % 初值
                p_init = p_true + radius_list(r_idx) * randn(3,1);
                theta_init = theta_true + deg2rad(angle_list(a_idx)) * randn(3,1);
                R_init = MatrixExp3(VecToso3(theta_init));

                % LM
                [p_lm, R_lm_, ~] = estimate_pose_lm( ...
                    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
                lm_pos_errs(r_idx,a_idx,pose_idx,init_idx) = norm(p_lm - p_true);
                lm_rot_errs(r_idx,a_idx,pose_idx,init_idx) = norm(R_lm_ - R_true, 'fro');

                % ELM
                [p_elm, R_elm_, ~] = estimate_pose_elm( ...
                    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
                elm_pos_errs(r_idx,a_idx,pose_idx,init_idx) = norm(p_elm - p_true);
                elm_rot_errs(r_idx,a_idx,pose_idx,init_idx) = norm(R_elm_ - R_true, 'fro');

                % Rlm
                [p_Rlm_, R_Rlm_, ~] = estimate_R_lm( ...
                    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options );
                Rlm_pos_errs(r_idx,a_idx,pose_idx,init_idx) = norm(p_Rlm_ - p_true);
                Rlm_rot_errs(r_idx,a_idx,pose_idx,init_idx) = norm(R_Rlm_ - R_true, 'fro');

                % Ours
                [p_ours_, R_ours_, ~] = estimate_pose_ours( ...
                    b_total, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p );
                ours_pos_errs(r_idx,a_idx,pose_idx,init_idx) = norm(p_ours_ - p_true);
                ours_rot_errs(r_idx,a_idx,pose_idx,init_idx) = norm(R_ours_ - R_true, 'fro');
            end
        end
    end
end

% 横轴：扰动旋量模长（angle_list），纵轴：估计误差（取所有半径平均）
mean_lm_pos = squeeze(mean(mean(lm_pos_errs,1),3));
mean_elm_pos = squeeze(mean(mean(elm_pos_errs,1),3));
mean_Rlm_pos = squeeze(mean(mean(Rlm_pos_errs,1),3));
mean_ours_pos = squeeze(mean(mean(ours_pos_errs,1),3));

mean_lm_rot = squeeze(mean(mean(lm_rot_errs,1),3));
mean_elm_rot = squeeze(mean(mean(elm_rot_errs,1),3));
mean_Rlm_rot = squeeze(mean(mean(Rlm_rot_errs,1),3));
mean_ours_rot = squeeze(mean(mean(ours_rot_errs,1),3));

figure('Name','扰动收敛域二维误差图','Color','white');
subplot(1,2,1);
plot(angle_list, mean_lm_pos, 'k-o', 'LineWidth',2); hold on;
plot(angle_list, mean_elm_pos, 'b-o', 'LineWidth',2);
plot(angle_list, mean_Rlm_pos, 'g-o', 'LineWidth',2);
plot(angle_list, mean_ours_pos, 'r-o', 'LineWidth',2);
xlabel('扰动旋量模长 (deg)');
ylabel('位置估计误差 (m)');
title('不同算法位置误差对比');
legend('LM','ELM','Rlm','Ours');
grid on;

subplot(1,2,2);
plot(angle_list, mean_lm_rot, 'k-o', 'LineWidth',2); hold on;
plot(angle_list, mean_elm_rot, 'b-o', 'LineWidth',2);
plot(angle_list, mean_Rlm_rot, 'g-o', 'LineWidth',2);
plot(angle_list, mean_ours_rot, 'r-o', 'LineWidth',2);
xlabel('扰动旋量模长 (deg)');
ylabel('旋转估计误差 (Frobenius)');
title('不同算法旋转误差对比');
legend('LM','ELM','Rlm','Ours');
grid on;


