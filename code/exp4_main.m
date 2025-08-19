% 磁传感器阵列沿直线移动下的位置估计实验（每次用真实位置附近扰动作为初值）

clc; clear; close all;
addpath('utils')
addpath('Functions')

rng(2025);

%% 磁铁参数定义
m_pos = [
    [-0.3; 0; 0], ...
    [0.2; 0; 0]
    ]; % 磁铁位置 [m]
m_hat = [
    [1; 0; 0], ...
    [0; 0; 1]
    ]; % 磁化方向（未归一化）
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

%% 沿直线移动的轨迹生成
num_points = 50;
start_p = [-0.1; 0; -0.05];
end_p   = [0.1; 0; -0.05];
p_traj = repmat(start_p, 1, num_points) + (0:num_points-1) .* repmat((end_p-start_p)/(num_points-1), 1, num_points);

theta_true = pi*rand(3,1); % 固定旋转向量
R_true = MatrixExp3(VecToso3(theta_true));

num_sensors = size(d_list, 2);
num_magnets = size(m_pos, 2);

results = struct();

% 优化器参数和约束提前定义
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

for idx = 1:num_points
    p_true = p_traj(:, idx);

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

    % 每次用真实位置附近扰动作为初值
    p_init = p_true + 0.01*randn(3,1);
    theta_init = theta_true + 0.1*randn(3,1);
    R_init = MatrixExp3(VecToso3(theta_init));

    % LM
    t_lm_start = tic;
    [p_lm, R_lm, stats_lm] = estimate_pose_lm( ...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
    t_lm = toc(t_lm_start);

    % ELM
    t_elm_start = tic;
    [p_elm, R_elm, stats_elm] = estimate_pose_elm( ...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
    t_elm = toc(t_elm_start);

    % 所提算法
    t_ours_start = tic;
    [p_ours, R_ours, stats_ours] = estimate_pose_ours( ...
        b_total, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p );
    t_ours = toc(t_ours_start);

    % Rlm
    t_Rlm_start = tic;
    [p_Rlm, R_Rlm, stats_Rlm] = estimate_R_lm( ...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );
    t_Rlm = toc(t_Rlm_start);

    % 结果存储
    results(idx).p_true = p_true;

    results(idx).p_lm    = p_lm;
    results(idx).R_lm    = R_lm;
    results(idx).lm_pos_error = norm(p_lm - p_true);
    results(idx).lm_rot_error = norm(R_lm - R_true, 'fro');
    results(idx).lm_t = t_lm;

    results(idx).p_elm   = p_elm;
    results(idx).R_elm   = R_elm;
    results(idx).elm_pos_error = norm(p_elm - p_true);
    results(idx).elm_rot_error = norm(R_elm - R_true, 'fro');
    results(idx).elm_t = t_elm;

    results(idx).p_Rlm   = p_Rlm;
    results(idx).R_Rlm   = R_Rlm;
    results(idx).Rlm_pos_error = norm(p_Rlm - p_true);
    results(idx).Rlm_rot_error = norm(R_Rlm - R_true, 'fro');
    results(idx).Rlm_t = t_Rlm;

    results(idx).p_ours  = p_ours;
    results(idx).R_ours  = R_ours;
    results(idx).ours_pos_error = norm(p_ours - p_true);
    results(idx).ours_rot_error = norm(R_ours - R_true, 'fro');
    results(idx).ours_t = t_ours;
end

% 可视化所有算法轨迹和误差
figure('Color','white','Name','沿直线移动位置估计误差');
p_true_mat = cell2mat({results.p_true});
p_lm_mat   = cell2mat({results.p_lm});
p_elm_mat  = cell2mat({results.p_elm});
p_Rlm_mat  = cell2mat({results.p_Rlm});
p_ours_mat = cell2mat({results.p_ours});

subplot(1,2,1);
plot3(p_true_mat(1,:), p_true_mat(2,:), p_true_mat(3,:), 'k.-', 'LineWidth', 2, 'MarkerSize', 16); hold on;
plot3(p_lm_mat(1,:),   p_lm_mat(2,:),   p_lm_mat(3,:),   'b.-', 'LineWidth', 1, 'MarkerSize', 10);
plot3(p_elm_mat(1,:),  p_elm_mat(2,:),  p_elm_mat(3,:),  'g.-', 'LineWidth', 1, 'MarkerSize', 10);
plot3(p_Rlm_mat(1,:),  p_Rlm_mat(2,:),  p_Rlm_mat(3,:),  'm.-', 'LineWidth', 1, 'MarkerSize', 10);
plot3(p_ours_mat(1,:), p_ours_mat(2,:), p_ours_mat(3,:), 'r.-', 'LineWidth', 2, 'MarkerSize', 10);
legend('真实轨迹','LM','ELM','Rlm','Ours');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('传感器阵列沿直线移动轨迹');
grid on; axis equal;

subplot(1,2,2);
plot([results.lm_pos_error],   'b.-', 'LineWidth', 1); hold on;
plot([results.elm_pos_error],  'g.-', 'LineWidth', 1);
plot([results.Rlm_pos_error],  'm.-', 'LineWidth', 1);
plot([results.ours_pos_error], 'r.-', 'LineWidth', 2);
xlabel('轨迹点'); ylabel('位置误差 (m)');
title('位置估计误差');
legend('LM','ELM','Rlm','Ours');
grid on;

sgtitle('磁传感器阵列沿直线移动下的位置估计实验（真实位置扰动初值）');
