function batch_results = run_batch_mu_experiments(params, test_points, num_trials)
% RUN_BATCH_EXPERIMENTS 批量执行实验（多个测试点，每个测试点多次实验）
% 输入：
%   params            - 实验参数结构体
%   test_points       - 测试点结构体数组（由 generate_test_points 生成）
%   num_trials_per_point - 每个测试点的实验次数（用于测试不同初始值）
% 输出：
%   batch_results     - 结构体数组，每个元素对应一个测试点：
%                       - test_point: 测试点信息（p_true, theta_true）
%                       - results: 该测试点的所有实验结果数组
%                       - summary: 统计摘要

num_test_points = length(test_points);
batch_results = struct('results', {});

fprintf('\n========== 开始批量实验 ==========\n');
fprintf('测试点数量: %d\n', num_test_points);
fprintf('每个测试点实验次数: %d\n', num_trials);
fprintf('总实验次数: %d\n\n', num_test_points * num_trials);

% 提取常用参数
d_list = params.sensor.d_list;
workspace_center = params.workspace.center;
workspace_radius = params.workspace.radius;
lb_p = workspace_center - workspace_radius;
lb_p(3) = 0;
ub_p = workspace_center + workspace_radius;

% 遍历每个测试点
mu_vec = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3];
W_vec = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1];
for trial_idx = 1:num_trials
    mu = mu_vec(trial_idx);
    params.optimization.mu = mu;

    % W = W_vec(trial_idx);
    % params.optimization.W(1, 1) = W;
    for point_idx = 1:num_test_points
        fprintf('\n========== 测试点 %d/%d ==========\n', point_idx, num_test_points);
        fprintf('mu = %.4f\n', mu);
        % fprintf('W = %.4f\n', W);
        p_t = test_points(point_idx).p_true;
        theta_t = test_points(point_idx).theta_true;
        fprintf('p_true = [%.4f, %.4f, %.4f]\n', p_t(1), p_t(2), p_t(3));
        fprintf('theta_true = [%.4f, %.4f, %.4f]\n', theta_t(1), theta_t(2), theta_t(3));
        
        % 更新参数中的真实姿态
        params_current = params;
        params_current.ground_truth.p_true = test_points(point_idx).p_true;
        params_current.ground_truth.theta_true = test_points(point_idx).theta_true;
        
        % 为当前测试点生成磁场数据
        [b_total, ~, ~, ~, ~] = generate_magnetic_data(params_current);
        b_total = b_total + kron(params_current.uncertainty.disturbance_strength, ones(1, size(params_current.sensor.d_list, 2))); % add noise to the magnetic field data
        
        % 计算真实旋转矩阵
        R_true = MatrixExp3(VecToso3(test_points(point_idx).theta_true));
        
        % 对该测试点执行多次实验
            exp_idx = (point_idx - 1) * num_trials + trial_idx;
            batch_results(trial_idx).results(point_idx) = run_single_experiment(exp_idx, ...
                num_test_points * num_trials, params_current, ...
                b_total, d_list, test_points(point_idx).p_true, R_true, lb_p, ub_p);
    end
end

fprintf('\n========== 批量实验完成 ==========\n');
