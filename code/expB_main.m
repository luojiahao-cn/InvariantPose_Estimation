%--------------------------------------------------------------------------
% beta参数实验：测试不同beta对所提算法定位性能的影响
% 输出：beta-误差关系图（均值±标准差）
%--------------------------------------------------------------------------

clc; clear; close all;
addpath('utils')
addpath('Functions')

rng(2025);

%% 参数设置
num_points = 50; % 每个beta采样数量
beta_list = logspace(-6, 6, 30); % beta取值范围
labels = {'位置误差','旋转误差'};
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

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off', ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6, ...
    'MaxIter', 1000, ...
    'MaxFunctionEvaluations', 10000);

workspace_center = [0; 0; 0];
workspace_radius = 0.5;
lb_p = workspace_center - workspace_radius;
ub_p = workspace_center + workspace_radius;

noise_level = 1e-6;
disturbance = [0; 0; 1e-6];

num_beta = numel(beta_list);

% 误差存储 [beta, 样本]
err_pos = zeros(num_beta, num_points);
err_rot = zeros(num_beta, num_points);
loss_ours = zeros(num_beta, num_points);
loss_rlm  = zeros(num_beta, num_points);

tic;
for beta_idx = 1:num_beta
    beta = beta_list(beta_idx);
    t_beta_start = tic;
    for pt_idx = 1:num_points
        % 随机生成真值
        theta_true = -pi + 2*pi*rand(3,1);
        p_true = -1e-1 + 2e-1*rand(3,1);
        R_true = MatrixExp3(VecToso3(theta_true));
        sensor_positions = p_true + R_true * d_list;

        % 生成磁场
        B_total = zeros(3, num_sensors);
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
        end
        b_total = R_true' * B_total;

        % 生成初值
        init_error = -1 + 2 * rand(3,1);
        p_init = p_true + 0.02 * init_error;
        theta_init = -pi + 2*pi*init_error;
        R_init = MatrixExp3(VecToso3(theta_init));

        % 所提算法，传入beta参数
        [p_ours, R_ours, ~] = estimate_pose_ours( ...
            b_total, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p, beta);

        % RLM算法
        [p_rlm, R_rlm, ~] = estimate_R_lm( ...
            b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options );

        err_pos(beta_idx, pt_idx) = norm(p_ours - p_true);
        err_rot(beta_idx, pt_idx) = norm(R_ours - R_true, 'fro');

        % 损失函数：R^T R与单位阵的距离
        loss_ours(beta_idx, pt_idx) = norm(R_ours' * R_ours - eye(3), 'fro');
        loss_rlm(beta_idx, pt_idx)  = norm(R_rlm'  * R_rlm  - eye(3), 'fro');
    end
    elapsed = toc(t_beta_start);
    percent = beta_idx / num_beta * 100;
    est_left = (num_beta-beta_idx)*elapsed;
    fprintf('beta %d/%d (%.1f%%) 用时 %.2fs，预计剩余 %.2fs\n', ...
        beta_idx, num_beta, percent, elapsed, est_left);
end
fprintf('全部beta实验完成，总耗时 %.2fs\n', toc);

%% 可视化
figure('Name','beta参数-位置误差','Color','white');
mu = mean(err_pos,2);
sigma = std(err_pos,0,2);
fill([beta_list, fliplr(beta_list)], [mu+sigma; flipud(mu-sigma)]', ...
    [0.2 0.6 1], 'FaceAlpha',0.15, 'EdgeColor','none'); hold on;
plot(beta_list, mu, 'b-o', 'LineWidth',2);
set(gca,'XScale','log');
xlabel('\beta');
ylabel('位置误差 (m)');
title('beta参数对位置误差的影响');
grid on;

figure('Name','beta参数-旋转误差','Color','white');
mu = mean(err_rot,2);
sigma = std(err_rot,0,2);
fill([beta_list, fliplr(beta_list)], [mu+sigma; flipud(mu-sigma)]', ...
    [1 0.5 0.2], 'FaceAlpha',0.15, 'EdgeColor','none'); hold on;
plot(beta_list, mu, 'r-o', 'LineWidth',2);
set(gca,'XScale','log');
xlabel('\beta');
ylabel('旋转误差 (Frobenius)');
title('beta参数对旋转误差的影响');
grid on;

figure('Name','beta参数-R^TR损失','Color','white');
mu_ours = mean(loss_ours,2);
sigma_ours = std(loss_ours,0,2);
mu_rlm = mean(loss_rlm,2);
sigma_rlm = std(loss_rlm,0,2);

fill([beta_list, fliplr(beta_list)], [mu_ours+sigma_ours; flipud(mu_ours-sigma_ours)]', ...
    [0.2 0.6 1], 'FaceAlpha',0.15, 'EdgeColor','none'); hold on;
plot(beta_list, mu_ours, 'b-o', 'LineWidth',2, 'DisplayName','Ours');
fill([beta_list, fliplr(beta_list)], [mu_rlm+sigma_rlm; flipud(mu_rlm-sigma_rlm)]', ...
    [1 0.5 0.2], 'FaceAlpha',0.15, 'EdgeColor','none');
plot(beta_list, mu_rlm, 'r-o', 'LineWidth',2, 'DisplayName','RLM');
set(gca,'XScale','log');
xlabel('\beta');
ylabel('$\|R^T R - I\|_F$', 'Interpreter','latex');
title('beta参数对R^T R损失的影响');
legend('Ours','RLM');
grid on;