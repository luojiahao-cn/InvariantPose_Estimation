%% 批量测试框架 - 测试多个 p_true 和 theta_true 组合
clc; clear; close all;
addpath('./utils')
addpath('./Functions')
addpath('./tools')
addpath('./exp')
%% ========== 参数配置 ==========
% params = get_experiment_params();
% 设置随机种子
% rng(params.experiment.random_seed);
% load('./exp/mat_data/optimized_params.mat', 'params');
load('./exp/mat_data/optimized_params_cone.mat', 'test_points', 'b_total', 'b_total_bg', 'params');
% b_total = b_total + b_total_bg; % 加回背景场，制造噪声

params_current = params;
params_current.uncertainty.p_uncertainty = 0.005;
params_current.uncertainty.r_uncertainty = 0.5;
params_current.workspace.radius = params_current.uncertainty.p_uncertainty;
params_current.optimization.W = eye(3);
params_current.optimization.options.FunctionTolerance = 1e-8;
params_current.optimization.options.StepTolerance = 1e-8;
params_current.optimization.mu = 1e-3;
params_current.sensor.d_list = params.sensor.d_list(:, [2,3,5]);
%% ========== 实验设置 ==========
num_trials_per_point = 1;  % 每个测试点的实验次数（用于测试不同初始值）
% 生成方法选项：
%   'random': 在工作空间内随机生成
%   'grid': 在工作空间内网格采样
%   'uniform_sphere': 在球面上均匀采样位置
% test_points = generate_test_points(params, num_test_points, 'random');

% fprintf('已生成 %d 个测试点\n', num_test_points);

%% ========== 批量执行实验 ==========
batch_results = run_batch_experiments(params_current, test_points(:, :), num_trials_per_point, b_total(:, [2,3,5], :));

%% ========== 结果分析 ==========
% 计算所有测试点的总平均值
all_summary = struct();
methods = {'lm', 'elm', 'ours', 'fischer', 'Rlm'};
for i = 1:length(methods)
    m = methods{i};
    all_summary.(m).pos_error = mean(arrayfun(@(x) x.summary.(m).pos_mean, batch_results));
    all_summary.(m).rot_error = mean(arrayfun(@(x) x.summary.(m).rot_mean, batch_results));
    all_summary.(m).r_error = mean(arrayfun(@(x) x.summary.(m).r_mean, batch_results));
end
all_summary.ours.direct_r_error_SDP = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_SDP, batch_results));
all_summary.ours.direct_r_error_RO  = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_RO, batch_results));
all_summary.ours.direct_r_error_spec = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_spec, batch_results));
all_summary.ours.direct_r_error_MM  = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_MM, batch_results));
all_summary.ours.direct_r_error_SSL = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_SSL, batch_results));

fprintf('\n========== 总平均误差 (1-abs(inner)) ==========\n');
fprintf('Method    | Pos Error (m) | R-axis Error\n');
fprintf('----------|---------------|--------------\n');
for i = 1:length(methods)
    m = methods{i};
    fprintf('%-9s | %13.6f | %12.6f\n', upper(m), all_summary.(m).pos_error, all_summary.(m).r_error);
end
fprintf('%-9s | %13s | %12.6f (SCA)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_MM);
fprintf('%-9s | %13s | %12.6f (SSL)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_SSL);
fprintf('%-9s | %13s | %12.6f (RO)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_RO);
fprintf('%-9s | %13s | %12.6f (SDP)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_SDP);
fprintf('%-9s | %13s | %12.6f (SPEC)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_spec);

% plot_path_comparison(batch_results, params_current);

% 绘图函数：
% 1. 标准版：生成全景图
% plot_batch_results(batch_results, [test_points.p_true], params);

% 2. 主轴对比版：重点展示 SCA, RO, SDP, SDP RED (ours), SPEC 的指向估计性能
% plot_principal_axis_results(batch_results, params_current);

plot_principal_axis_results_cone(batch_results, params_current);

% 3. 文本标注版：将误差统计信息 (mean ± std) 直接显示在子图标题中
% plot_batch_results_text(batch_results, [test_points.p_true], params);

% 4. 2D 投影版：仅绘制当前层 (highlight)，保留 ep/er 的 boxchart
% plot_batch_results_2d(batch_results, [test_points.p_true], params_current);

% %% ========== 保存结果 ==========
% save('results/exp1_batch_results.mat', 'batch_results', 'test_points', 'params');

% fprintf('\n结果已保存到 results/exp1_batch_results.mat\n');