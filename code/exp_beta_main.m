%% 批量测试框架 - 测试多个 p_true 和 theta_true 组合
clc; clear; close all;
addpath('./utils')
addpath('./Functions')
addpath('./tools')

%% ========== 参数配置 ==========
params = get_experiment_params();

% 设置随机种子
rng(params.experiment.random_seed);

%% ========== 生成测试点 ==========
% 生成100个测试点（p_true 和 theta_true 的组合）
num_test_points = 9;  % 测试点数量
num_trials_per_point = 1;  % 每个测试点的实验次数（用于测试不同初始值）

% 生成方法选项：
%   'random': 在工作空间内随机生成
%   'grid': 在工作空间内网格采样
%   'uniform_sphere': 在球面上均匀采样位置
test_point_fixed = generate_test_points(params, 1, 'fixed');
test_points = generate_test_points(params, num_test_points, 'random');
test_points = [test_point_fixed, test_points];

fprintf('已生成 %d 个测试点\n', num_test_points);

%% ========== 批量执行实验 ==========
batch_results = run_beta_batch_experiment(params, test_points, num_trials_per_point);

%% ========== 结果分析 ==========
analyze_beta_results(batch_results);

%% ========== 结果可视化 ==========
% plot_beta_results(batch_results);

