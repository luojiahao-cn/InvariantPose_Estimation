function batch_results = run_batch_experiments(params, test_points, num_trials_per_point, b_total)
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
workspace_radius = params.workspace.radius;

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
    
    % 计算真实旋转矩阵
    R_true = MatrixExp3(VecToso3(test_points(point_idx).theta_true));
    
    % 使用元胞数组收集结果，避开繁琐且易错的结构体模板
    results_cell = cell(1, num_trials_per_point);
    
    % 提取常量变量
    total_experiments_count = num_test_points * num_trials_per_point;
    p_true_current = test_points(point_idx).p_true;

    for trial_idx = 1:num_trials_per_point
        current_exp_idx = (point_idx - 1) * num_trials_per_point + trial_idx;
        lb_p = p_true_current - workspace_radius;
        ub_p = p_true_current + workspace_radius;
        
        % 磁场索引逻辑
        if ndims(b_total) == 4
            b_input = b_total(:, :, point_idx, trial_idx);
        elseif size(b_total, 3) == total_experiments_count
            b_input = b_total(:, :, current_exp_idx);
        else
            b_input = b_total(:, :, point_idx);
        end

        results_cell{trial_idx} = run_single_experiment(current_exp_idx, ...
            total_experiments_count, params_current, ...
            b_input, d_list, p_true_current, R_true, lb_p, ub_p);
    end

    % 将 cell 转换回结构体数组
    results_struct_array = [results_cell{:}];
    
    % 将结构体数组转换为 Table，处理不同维度的字段（如 3x1 向量或 3x3 矩阵）
    % 我们使用这种方式强制将所有内容视为单行条目
    results_table = struct2table(results_struct_array, 'AsArray', isscalar(results_struct_array));
    
    % 沿着 trials 计算平均值和汇总
    summary = calculate_point_summary_from_table(results_table);
        
    % 存储结果（转换为 struct 数组以保持与原有绘图代码兼容）
    batch_results(point_idx).test_point = test_points(point_idx);
    batch_results(point_idx).results = results_struct_array;
    batch_results(point_idx).summary = summary;
end

fprintf('\n========== 批量实验完成 ==========\n');

end

function summary = calculate_point_summary_from_table(T)
    % 使用 Table 的向量化操作进行极简统计
    methods = {'lm', 'elm', 'ours', 'fischer', 'Rlm'};
    summary = struct();
    var_names = T.Properties.VariableNames;

    for i = 1:length(methods)
        m = methods{i};
        
        % 自动匹配字段名并计算均值/标准差/最大值
        f_pos = [m '_pos_error'];
        f_rot = [m '_rot_error'];
        f_r   = [m '_r_error'];
        f_time = ['time_' m];
        
        if ismember(f_pos, var_names)
            summary.(m).pos_mean = mean(T.(f_pos));
            summary.(m).pos_std  = std(T.(f_pos));
            summary.(m).pos_max  = max(T.(f_pos));
        end
        if ismember(f_rot, var_names)
            summary.(m).rot_mean = mean(T.(f_rot));
            summary.(m).rot_std  = std(T.(f_rot));
            summary.(m).rot_max  = max(T.(f_rot));
        end
        if ismember(f_r, var_names)
            summary.(m).r_mean = mean(T.(f_r));
            summary.(m).r_std  = std(T.(f_r));
        end
        if ismember(f_time, var_names)
            summary.(m).time_mean = mean(T.(f_time));
            summary.(m).time_std  = std(T.(f_time));
        end
    end
    
    % 特殊字段处理 (Direct R)
    if ismember('direct_r_error_MM', var_names)
        summary.ours.direct_r_mean_MM = mean(T.direct_r_error_MM);
    end
    if ismember('direct_r_error_RO', var_names)
        summary.ours.direct_r_mean_RO = mean(T.direct_r_error_RO);
    end
    if ismember('direct_r_error_SDP', var_names)
        summary.ours.direct_r_mean_SDP = mean(T.direct_r_error_SDP);
    end
end

