function batch_results = run_beta_batch_experiment(params, test_points, num_trials_per_point)
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
batch_results = struct('test_point', {}, 'results', {});

fprintf('\n========== 开始批量实验 ==========\n');
fprintf('测试点数量: %d\n', num_test_points);
fprintf('每个测试点实验次数: %d\n', num_trials_per_point);
fprintf('总实验次数: %d\n\n', num_test_points * num_trials_per_point);

% 提取常用参数
d_list = params.sensor.d_list;

% 遍历每个测试点
for point_idx = 1:num_test_points
    fprintf('\n========== 测试点 %d/%d ==========\n', point_idx, num_test_points);
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
    
    % 初始化结果元胞数组
    results = cell(num_trials_per_point, 1);
    
    unc = linspace(0.02, 0.1, num_trials_per_point);
    for trial_idx = 1:num_trials_per_point
        % params_current.experiment.p_uncertainty = unc(trial_idx); % 位置不确定度
        params_current.experiment.p_uncertainty = 0.04;
        % params_current.experiment.r_uncertainty = unc(trial_idx);  % 旋转不确定度
        params_current.experiment.r_uncertainty = 0.1;
        exp_idx = (point_idx - 1) * num_trials_per_point + trial_idx;
        results{trial_idx} = run_beta_experiment(exp_idx, ...
            num_test_points * num_trials_per_point, params_current, ...
            b_total, d_list, test_points(point_idx).p_true);
    end
    
    % 存储结果
    batch_results(point_idx).test_point = test_points(point_idx);
    batch_results(point_idx).results = results;
end

fprintf('\n========== 批量实验完成 ==========\n');

end
