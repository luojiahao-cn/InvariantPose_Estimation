function batch_results = run_batch_experiments(params, test_points, num_trials_per_point)
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
batch_results = struct('test_point', {}, 'results', {}, 'summary', {});

fprintf('\n========== 开始批量实验 ==========\n');
fprintf('测试点数量: %d\n', num_test_points);
fprintf('每个测试点实验次数: %d\n', num_trials_per_point);
fprintf('总实验次数: %d\n\n', num_test_points * num_trials_per_point);

% 提取常用参数
d_list = params.sensor.d_list;
workspace_center = params.workspace.center;
workspace_radius = params.workspace.radius;
lb_p = workspace_center - workspace_radius;
lb_p(3) = 0;
ub_p = workspace_center + workspace_radius;

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
    % 静磁场下LM就炸了
    b_total = b_total + kron(params_current.uncertainty.disturbance_strength, ones(1, size(params_current.sensor.d_list, 2))); % add noise to the magnetic field data
    
    % 计算真实旋转矩阵
    R_true = MatrixExp3(VecToso3(test_points(point_idx).theta_true));
    
    % 对该测试点执行多次实验
    % 预先初始化结果结构体数组，确保字段一致
    results = struct('p_lm', [], 'R_lm', [], 'p_elm', [], 'R_elm', [],...
     'p_ours', [], 'R_ours', [], 'p_Rlm', [], 'R_Rlm', [],...
      'lm_pos_error', [], 'lm_rot_error', [],...
       'elm_pos_error', [], 'elm_rot_error', [], ... 
       'ours_pos_error', [], 'ours_rot_error', [],...
        'Rlm_pos_error', [], 'Rlm_rot_error', [],...
        'p_init', [], 'theta_init', [],...
        'p_fischer', [], 'R_fischer', [], 'fischer_pos_error', [], 'fischer_rot_error', [],...
        'time_lm', [], 'time_elm', [], 'time_ours', [], 'time_fischer', [], 'time_Rlm', []);
    
    for trial_idx = 1:num_trials_per_point
        exp_idx = (point_idx - 1) * num_trials_per_point + trial_idx;
        results(trial_idx) = run_single_experiment(exp_idx, ...
            num_test_points * num_trials_per_point, params_current, ...
            b_total, d_list, test_points(point_idx).p_true, R_true, lb_p, ub_p);
    end

    % 沿着num_trials_per_point计算平均
    summary = calculate_point_summary(results);
        
    % 存储结果
    batch_results(point_idx).test_point = test_points(point_idx);
    batch_results(point_idx).results = results;
    batch_results(point_idx).summary = summary;
end

fprintf('\n========== 批量实验完成 ==========\n');

end

function summary = calculate_point_summary(results)
    % 计算单个测试点的统计摘要
    % 提取所有误差
    lm_pos_errors = [results.lm_pos_error];
    lm_rot_errors = [results.lm_rot_error];
    elm_pos_errors = [results.elm_pos_error];
    elm_rot_errors = [results.elm_rot_error];
    ours_pos_errors = [results.ours_pos_error];
    ours_rot_errors = [results.ours_rot_error];
    fischer_pos_errors = [results.fischer_pos_error];
    fischer_rot_errors = [results.fischer_rot_error];
    Rlm_pos_errors = [results.Rlm_pos_error];
    Rlm_rot_errors = [results.Rlm_rot_error];
    time_lm = [results.time_lm];
    time_elm = [results.time_elm];
    time_ours = [results.time_ours];
    time_fischer = [results.time_fischer];
    time_Rlm = [results.time_Rlm];

    % 计算统计量
    summary.lm.pos_mean = mean(lm_pos_errors);
    summary.lm.pos_std = std(lm_pos_errors);
    summary.lm.pos_max = max(lm_pos_errors);
    summary.lm.rot_mean = mean(lm_rot_errors);
    summary.lm.rot_std = std(lm_rot_errors);
    summary.lm.rot_max = max(lm_rot_errors);
    summary.lm.time_mean = mean(time_lm);
    summary.lm.time_std = std(time_lm);
    summary.lm.time_max = max(time_lm);

    summary.elm.pos_mean = mean(elm_pos_errors);
    summary.elm.pos_std = std(elm_pos_errors);
    summary.elm.pos_max = max(elm_pos_errors);
    summary.elm.rot_mean = mean(elm_rot_errors);
    summary.elm.rot_std = std(elm_rot_errors);
    summary.elm.rot_max = max(elm_rot_errors);
    summary.elm.time_mean = mean(time_elm);
    summary.elm.time_std = std(time_elm);
    summary.elm.time_max = max(time_elm);

    summary.ours.pos_mean = mean(ours_pos_errors);
    summary.ours.pos_std = std(ours_pos_errors);
    summary.ours.pos_max = max(ours_pos_errors);
    summary.ours.rot_mean = mean(ours_rot_errors);
    summary.ours.rot_std = std(ours_rot_errors);
    summary.ours.rot_max = max(ours_rot_errors);
    summary.ours.time_mean = mean(time_ours);
    summary.ours.time_std = std(time_ours);
    summary.ours.time_max = max(time_ours);

    summary.fischer.pos_mean = mean(fischer_pos_errors);
    summary.fischer.pos_std = std(fischer_pos_errors);
    summary.fischer.pos_max = max(fischer_pos_errors);
    summary.fischer.rot_mean = mean(fischer_rot_errors);
    summary.fischer.rot_std = std(fischer_rot_errors);
    summary.fischer.rot_max = max(fischer_rot_errors);
    summary.fischer.time_mean = mean(time_fischer);
    summary.fischer.time_std = std(time_fischer);
    summary.fischer.time_max = max(time_fischer);

    summary.Rlm.pos_mean = mean(Rlm_pos_errors);
    summary.Rlm.pos_std = std(Rlm_pos_errors);
    summary.Rlm.pos_max = max(Rlm_pos_errors);
    summary.Rlm.rot_mean = mean(Rlm_rot_errors);
    summary.Rlm.rot_std = std(Rlm_rot_errors);
    summary.Rlm.rot_max = max(Rlm_rot_errors);
    summary.Rlm.time_mean = mean(time_Rlm);
    summary.Rlm.time_std = std(time_Rlm);
    summary.Rlm.time_max = max(time_Rlm);

end

