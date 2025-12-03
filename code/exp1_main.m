%% 基本优化框架
clc; clear; close all;
addpath('./utils')
addpath('./Functions')

%% ========== 参数配置 ==========
% 所有实验参数统一在 get_experiment_params() 函数中配置
% 如需修改参数，请编辑 code/utils/get_experiment_params.m 文件
params = get_experiment_params();

% 设置随机种子
rng(params.experiment.random_seed);

%% ========== 生成磁铁测量数据 ==========
[b_total, B_total, gradb_total, gradB_total, sensor_positions] = ...
    generate_magnetic_data(params);

% 提取常用参数
m_pos = params.magnet.m_pos;
m_hat = params.magnet.m_hat;
m_norm = params.magnet.m_norm;
d_list = params.sensor.d_list;
p_true = params.ground_truth.p_true;
theta_true = params.ground_truth.theta_true;
R_true = MatrixExp3(VecToso3(theta_true));

%% ========== 工作空间约束设置 ==========
workspace_center = params.workspace.center;
workspace_radius = params.workspace.radius;
lb_p = workspace_center - workspace_radius; % 下界
ub_p = workspace_center + workspace_radius; % 上界

%% ========== 多次实验执行 ==========
num_experiments = params.experiment.num_experiments;
results = struct();

for exp_idx = 1:num_experiments
    results(exp_idx) = run_single_experiment(exp_idx, num_experiments, params, ...
        b_total, d_list, p_true, R_true, lb_p, ub_p);
end

%% ========== 结果可视化与分析 ==========
visualize_pose(m_pos, m_hat, m_norm, p_true, R_true, results, d_list);
display_statistical_summary(results, num_experiments);
plot_error_distributions(results);


