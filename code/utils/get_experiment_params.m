function params = get_experiment_params()
% GET_EXPERIMENT_PARAMS 获取实验参数配置
% 输出：
%   params - 包含所有实验参数的结构体，分为以下子结构：
%     - magnet: 磁铁相关参数
%     - sensor: 传感器相关参数
%     - ground_truth: 真实姿态参数
%     - uncertainty: 不确定性参数
%     - optimization: 优化算法参数
%     - workspace: 工作空间约束参数
%     - experiment: 实验设置参数

%% 默认磁铁参数
params.magnet.m_pos = [
    [-0.1125;0.0002;0.3099], ...
    [0.1125;0.0002;0.3099]
    ];
params.magnet.m_hat = [
    [-0.2588;0;0.9659], ...
    [0.2588;0;0.9659]
    ];
params.magnet.m_hat = params.magnet.m_hat ./ vecnorm(params.magnet.m_hat);
params.magnet.m_norm = [300, 300];

%% 默认传感器参数
params.sensor.d_list = [
    [0; 0; 0], ...
    [1e-3; 0; 0], ...
    [-1e-3; 0; 0], ...
    [0; 1e-3; 0], ...
    [0; -1e-3; 0], ...
    [0; 0; 1e-3], ...
    [0; 0; -1e-3]
    ];
%% 默认真实位姿
params.ground_truth.theta_true = [0; 0; 1]; % 真实旋转向量 [rad]
params.ground_truth.p_true = [0; 0; 0]; % 传感器阵列参考点真实位置 [m]
%% 不确定性参数
params.uncertainty.p_uncertainty = 0.05; % 位置不确定性
params.uncertainty.r_uncertainty = 0.1;  % 旋转不确定性

%% 优化算法参数
params.optimization.options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off');
params.optimization.mu = 1;   % 所提算法的mu参数
params.optimization.beta = 1e-3;   % 所提算法的beta参数

%% 工作空间约束参数
params.workspace.center = [0; 0; 0];
params.workspace.radius = 0.15;

%% 实验设置参数
params.experiment.num_experiments = 10; % 实验次数
params.experiment.random_seed = 2025;  % 随机种子

end

