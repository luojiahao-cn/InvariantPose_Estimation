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
params_current.uncertainty.p_uncertainty = 0.01;
params_current.uncertainty.r_uncertainty = 0.5;
params_current.workspace.radius = params_current.uncertainty.p_uncertainty;
params_current.optimization.W = eye(3);
params_current.optimization.options.FunctionTolerance = 1e-8;
params_current.optimization.options.StepTolerance = 1e-8;
params_current.optimization.mu = 1e-3;
% d_list_ind = [1:12];
d_list_ind = [1:6];
% d_list_ind = [2,3,6];
params_current.sensor.d_list = params.sensor.d_list(:, d_list_ind);
% params_current.sensor.d_list = params.sensor.d_list(:, [1,2,3,7,8,9]);
%% ========== 实验设置 ==========
num_trials_per_point = 1;  % 每个测试点的实验次数（用于测试不同初始值）
% 生成方法选项：
%   'random': 在工作空间内随机生成
%   'grid': 在工作空间内网格采样
%   'uniform_sphere': 在球面上均匀采样位置
% test_points = generate_test_points(params, num_test_points, 'random');

% fprintf('已生成 %d 个测试点\n', num_test_points);

%% ========== 批量执行实验 ==========
batch_results = run_batch_experiments(params_current, test_points(:, :), num_trials_per_point, b_total(:, d_list_ind, :));

%% ========== 结果分析 ==========
% 计算所有测试点的总平均值
all_summary = struct();
methods = {'lm', 'elm', 'ours', 'fischer', 'Rlm'};
for i = 1:length(methods)
    m = methods{i};
    all_summary.(m).pos_error = mean(arrayfun(@(x) x.summary.(m).pos_mean, batch_results));
    all_summary.(m).rot_error = mean(arrayfun(@(x) x.summary.(m).rot_mean, batch_results));
    all_summary.(m).r_error = mean(arrayfun(@(x) x.summary.(m).r_mean, batch_results));
    all_summary.(m).time = mean(arrayfun(@(x) x.summary.(m).time_mean, batch_results));
end
all_summary.ours.direct_r_error_SDP = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_SDP, batch_results));
all_summary.ours.direct_r_error_RGD = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_RGD, batch_results));
all_summary.ours.direct_r_error_spec = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_spec, batch_results));
all_summary.ours.direct_r_error_MM  = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_MM, batch_results));
all_summary.ours.direct_r_error_SSL = mean(arrayfun(@(x) x.summary.ours.direct_r_mean_SSL, batch_results));

% 时间统计：基于所有 trial 的原始时间（更准确）
all_trials = [batch_results.results];

all_summary.lm.time = mean([all_trials.time_lm]);
all_summary.lm.time_std = std([all_trials.time_lm]);
all_summary.elm.time = mean([all_trials.time_elm]);
all_summary.elm.time_std = std([all_trials.time_elm]);
all_summary.ours.time = mean([all_trials.time_ours]);
all_summary.ours.time_std = std([all_trials.time_ours]);
all_summary.fischer.time = mean([all_trials.time_fischer]);
all_summary.fischer.time_std = std([all_trials.time_fischer]);
all_summary.Rlm.time = mean([all_trials.time_Rlm]);
all_summary.Rlm.time_std = std([all_trials.time_Rlm]);

all_summary.ours.direct_time_SCA = mean([all_trials.time_SCA]);
all_summary.ours.direct_time_SCA_std = std([all_trials.time_SCA]);
all_summary.ours.direct_time_SSL = mean([all_trials.time_SSL]);
all_summary.ours.direct_time_SSL_std = std([all_trials.time_SSL]);
all_summary.ours.direct_time_RGD = mean([all_trials.time_RO]);
all_summary.ours.direct_time_RGD_std = std([all_trials.time_RO]);

fprintf('\n========== 总平均误差 (1-abs(inner)) ==========\n');
fprintf('Method    | Pos Error (m) | R-axis Error\n');
fprintf('----------|---------------|--------------\n');
for i = 1:length(methods)
    m = methods{i};
    fprintf('%-9s | %13.6f | %12.6f\n', upper(m), all_summary.(m).pos_error, all_summary.(m).r_error);
end
fprintf('%-9s | %13s | %12.6f (SCA)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_MM);
fprintf('%-9s | %13s | %12.6f (SSL)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_SSL);
fprintf('%-9s | %13s | %12.6f (RGD)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_RGD);
fprintf('%-9s | %13s | %12.6f (SDP)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_SDP);
fprintf('%-9s | %13s | %12.6f (SPEC)\n', 'OURS_DIR', '-', all_summary.ours.direct_r_error_spec);

fprintf('\n========== Direct Principal Axis 平均耗时 (mean ± std, s) ==========\n' );
fprintf('SSL  : %.6e ± %.6e\n', all_summary.ours.direct_time_SSL, all_summary.ours.direct_time_SSL_std);
fprintf('SCA  : %.6e ± %.6e\n', all_summary.ours.direct_time_SCA, all_summary.ours.direct_time_SCA_std);
fprintf('RGD  : %.6e ± %.6e\n', all_summary.ours.direct_time_RGD, all_summary.ours.direct_time_RGD_std);

fprintf('\n========== 算法平均耗时对比 (mean ± std, s) ==========\n' );
fprintf('OURS    : %.6e ± %.6e\n', all_summary.ours.time, all_summary.ours.time_std);
fprintf('LM      : %.6e ± %.6e\n', all_summary.lm.time, all_summary.lm.time_std);
fprintf('ELM     : %.6e ± %.6e\n', all_summary.elm.time, all_summary.elm.time_std);
fprintf('FISCHER : %.6e ± %.6e\n', all_summary.fischer.time, all_summary.fischer.time_std);
fprintf('RLM     : %.6e ± %.6e\n', all_summary.Rlm.time, all_summary.Rlm.time_std);

%% ========== 保存时间统计到 CSV ==========
if ~exist('results', 'dir')
    mkdir('results');
end

num_sensors = length(d_list_ind);
num_points = length(test_points);
timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

algorithm_names = {'LM'; 'ELM'; 'OURS'; 'FISCHER'; 'RLM'; 'SCA'; 'SSL'; 'RGD'};
time_mean_s = [ ...
    all_summary.lm.time; ...
    all_summary.elm.time; ...
    all_summary.ours.time; ...
    all_summary.fischer.time; ...
    all_summary.Rlm.time; ...
    all_summary.ours.direct_time_SCA; ...
    all_summary.ours.direct_time_SSL; ...
    all_summary.ours.direct_time_RGD];
time_std_s = [ ...
    all_summary.lm.time_std; ...
    all_summary.elm.time_std; ...
    all_summary.ours.time_std; ...
    all_summary.fischer.time_std; ...
    all_summary.Rlm.time_std; ...
    all_summary.ours.direct_time_SCA_std; ...
    all_summary.ours.direct_time_SSL_std; ...
    all_summary.ours.direct_time_RGD_std];

T_time = table( ...
    repmat(timestamp, numel(algorithm_names), 1), ...
    repmat(num_sensors, numel(algorithm_names), 1), ...
    repmat(num_points, numel(algorithm_names), 1), ...
    repmat(num_trials_per_point, numel(algorithm_names), 1), ...
    algorithm_names, ...
    time_mean_s, ...
    time_std_s, ...
    'VariableNames', {'timestamp', 'num_sensors', 'num_points', 'num_trials_per_point', 'algorithm', 'time_mean_s', 'time_std_s'});

csv_path = sprintf('results/time_summary_N%d.csv', num_sensors);
if exist(csv_path, 'file')
    writetable(T_time, csv_path, 'WriteMode', 'append', 'WriteVariableNames', false);
else
    writetable(T_time, csv_path);
end

fprintf('\n时间统计已保存到: %s\n', csv_path);


% plot_path_comparison(batch_results, params_current);

% 绘图函数：
% 1. 标准版：生成全景图
% plot_batch_results(batch_results, [test_points.p_true], params);

% 2. 主轴对比版：重点展示 SCA, RGD, SDP, SDP RED (ours), SPEC 的指向估计性能
% plot_principal_axis_results(batch_results, params_current);

% plot_principal_axis_results_cone(batch_results, params_current);

% 3. 文本标注版：将误差统计信息 (mean ± std) 直接显示在子图标题中
% plot_batch_results_text(batch_results, [test_points.p_true], params);

% 4. 2D 投影版：仅绘制当前层 (highlight)，保留 ep/er 的 boxchart
% plot_batch_results_2d(batch_results, [test_points.p_true], params_current);

% %% ========== 保存结果 ==========
% save('results/exp1_batch_results.mat', 'batch_results', 'test_points', 'params');

% fprintf('\n结果已保存到 results/exp1_batch_results.mat\n');