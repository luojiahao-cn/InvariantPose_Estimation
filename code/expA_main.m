%--------------------------------------------------------------------------
% 初值实验：测试不同干扰/噪声条件下的定位能力
% 条件：无干扰噪声、仅干扰、仅噪声、干扰和噪声
% 输出：误差分布、误差带图
%--------------------------------------------------------------------------

clc; clear; close all;
addpath('utils')
addpath('Functions')

rng(2025);

%% 参数设置
num_points = 10; % 每个偏离度采样数量
offset_list = linspace(0, 0.1, 51); % 初值偏离度（米）
labels = {'无干扰噪声','仅干扰','仅噪声','干扰和噪声'};
alg_labels = {'LM','ELM','Rlm','Ours'};
alg_colors = lines(numel(alg_labels));

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

% options = optimoptions('lsqnonlin', ...
%     'Algorithm', 'levenberg-marquardt', ...
%     'Display', 'off', ...
%     'TolFun', 1e-6, ...
%     'TolX', 1e-6, ...
%     'MaxIter', 1000, ...
%     'MaxFunctionEvaluations', 10000);

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off');

workspace_center = [0; 0; 0];
workspace_radius = 0.5;
lb_p = workspace_center - workspace_radius;
ub_p = workspace_center + workspace_radius;

noise_level = 3e-5;
disturbance = [0; 0; 2e-3];

num_alg = numel(alg_labels);
num_cond = numel(labels);
num_offset = numel(offset_list);

% 误差存储 [偏离度, 条件, 算法, 样本]
err_pos = zeros(num_offset, num_cond, num_alg, num_points);
err_rot = zeros(num_offset, num_cond, num_alg, num_points);

tic; % 总计时
hwb = waitbar(0, '正在运行实验，请稍候...');
for off_idx = 1:num_offset
    offset = offset_list(off_idx);
    t_offset_start = tic;
    for pt_idx = 1:num_points
        % waitbar进度更新，增加已用时和剩余时间
        total_sample = (off_idx-1)*num_points + pt_idx;
        total_samples = num_offset * num_points;
        percent_sample = total_sample / total_samples;
        elapsed = toc;
        if total_sample > 1
            avg_time = elapsed / total_sample;
            est_left = (total_samples - total_sample) * avg_time;
        else
            est_left = 0;
        end
        waitbar(percent_sample, hwb, ...
            sprintf('进度 %.2f%% (%d/%d)\n已用时 %.1fs，预计剩余 %.1fs', ...
                percent_sample*100, total_sample, total_samples, elapsed, est_left));

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
        p_init = p_true + offset * init_error; 
        theta_init = -pi + 2*pi*init_error;
        theta_init = theta_true + offset * 10 * theta_init;
        R_init = MatrixExp3(VecToso3(theta_init));

        % 四种情况
        for mode = 1:num_cond
            b_meas = b_total;
            switch mode
                case 1 % 无干扰噪声
                    % 不做处理
                case 2 % 仅干扰
                    b_meas = b_meas + disturbance;
                case 3 % 仅噪声
                    b_meas = b_meas + randn(size(b_meas)) * noise_level;
                case 4 % 干扰和噪声
                    b_meas = b_meas + disturbance + randn(size(b_meas)) * noise_level;
            end

            % 多算法比较
            % LM
            [p_lm, R_lm, ~] = estimate_pose_lm( ...
                b_meas, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
            err_pos(off_idx, mode, 1, pt_idx) = norm(p_lm - p_true);
            err_rot(off_idx, mode, 1, pt_idx) = norm(R_lm - R_true, 'fro');

            % ELM
            [p_elm, R_elm, ~] = estimate_pose_elm( ...
                b_meas, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
            err_pos(off_idx, mode, 2, pt_idx) = norm(p_elm - p_true);
            err_rot(off_idx, mode, 2, pt_idx) = norm(R_elm - R_true, 'fro');

            % Ours
            [p_ours, R_ours, ~] = estimate_pose_ours( ...
                b_meas, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p );
            err_pos(off_idx, mode, 4, pt_idx) = norm(p_ours - p_true);
            err_rot(off_idx, mode, 4, pt_idx) = norm(R_ours - R_true, 'fro');
            
            % Rlm
            [p_Rlm, R_Rlm, ~] = estimate_R_lm( ...
                b_meas, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );
            err_pos(off_idx, mode, 3, pt_idx) = norm(p_Rlm - p_true);
            err_rot(off_idx, mode, 3, pt_idx) = norm(R_Rlm - R_true, 'fro');

        end
    end
    elapsed_offset = toc(t_offset_start);
    total_elapsed = toc;
    percent = off_idx / num_offset * 100;
    est_left = (num_offset-off_idx)*elapsed_offset;
    fprintf('偏离度 %d/%d (%.1f%%) 用时 %.2fs，预计剩余 %.2fs\n', ...
        off_idx, num_offset, percent, elapsed_offset, est_left);
end
close(hwb);
fprintf('全部实验完成，总耗时 %.2fs\n', toc);

save('../results/expA_main.mat', 'err_pos', 'err_rot', 'offset_list', 'labels', 'alg_labels', 'alg_colors', 'num_alg', 'num_cond', 'num_offset', 'num_points');