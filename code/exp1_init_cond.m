%% 批量测试框架 - 测试多个 p_true 和 theta_true 组合
clc; clear; close all;
addpath('./utils')
addpath('./Functions')
addpath('./tools')

d0 = 1.5e-3;
d_list1 = [[0; 0; 0],...
    [d0; 0; 0],...
    [0; d0; 0],...
    [d0; d0; 0],...
    [0; 2*d0; 0],...
    [d0; 2*d0; 0],...
    [0; 0; -d0],...
    [d0; 0; -d0],...
    [0; d0; -d0],...
    [d0; d0; -d0],...
    [0; 2*d0; -d0],...
    [d0; 2*d0; -d0]]; % redundant config 1

% 计算几何中心
dc1 = mean(d_list1, 2);
d_list1 = d_list1 - dc1;

d_list2 = [[0; 0; 0],...
    [0; d0; 0],...
    [0; 0; -d0],...
    [0; d0; -d0]]; % redundant config 2
dc2 = mean(d_list2, 2);
d_list2 = d_list2 - dc2;

d_list3 = [[0; 0; 0],...
    [0; d0; 0],...
    [0; d0; -d0]]; % minimal config 3
dc3 = [0; d0/2; -d0/2];
d_list3 = d_list3 - dc3;

d_list4 = [[-d0; 0; 0],...
    [0; d0; 0],...
    [0; 0; d0]];
dc4 = [0; 0; 0];
d_list4 = d_list4 - dc4;

%% ========== 参数配置 ==========
params = get_experiment_params();
params.sensor.d_list = d_list1;
% 设置随机种子
rng(params.experiment.random_seed);

%% ========== 生成测试点 ==========
% 生成100个测试点（p_true 和 theta_true 的组合）
num_test_points = 300;  % 测试点数量
num_trials_per_point = 1;  % 每个测试点的实验次数（用于测试不同初始值）

% 生成方法选项：
%   'random': 在工作空间内随机生成
%   'grid': 在工作空间内网格采样
%   'uniform_sphere': 在球面上均匀采样位置
test_points = generate_test_points(params, num_test_points, 'fixed');

fprintf('已生成 %d 个测试点\n', num_test_points);

%% ========== 批量执行实验 ==========
batch_results = run_batch_experiments(params, test_points, num_trials_per_point);

%% ========== 结果分析 ==========
analyze_batch_results(batch_results);

%% ========== 结果可视化 ==========
plot_batch_results(batch_results);

% %% ========== 保存结果 ==========
% save('results/exp1_batch_results.mat', 'batch_results', 'test_points', 'params');

% fprintf('\n结果已保存到 results/exp1_batch_results.mat\n');