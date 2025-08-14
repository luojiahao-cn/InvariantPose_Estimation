clc; clear; close all;
addpath('utils')
addpath('Functions')

% 设置随机种子
rng(2025);

%% 磁铁参数定义
m_pos = [
    [-0.3; 0; 0] ,...
    [0.2; 0; 0];
    ]; % 磁铁位置 [m]
m_hat = [
    [1; 0; 0],...
    [0; 0; 1]
    ];  % 磁化方向（未归一化）

m_hat = m_hat ./ vecnorm(m_hat); % 归一化磁化方向
m_norm = [8e2 , 8e2];            % 磁矩幅值 [A·m²]

d_list_e = [
    [0; 0; 0],...
    [1e-3; 0; 0],...
    [-1e-3; 0; 0],...
    [0; 1e-3; 0],...
    [0; -1e-3; 0],...
    [0; 0; 1e-3],...
    [0; 0; -1e-3],...
    ];
row_means = mean(d_list_e, 2);
d_list = d_list_e - row_means;

% 真实姿态
theta_true = pi*rand(3,1); % 真实旋转向量 [rad]
p_true = -1e-2 + 2e-2 * rand(3,1); %[0; 0; -0.05]; % 传感器阵列参考点真实位置 [m]
R_true = MatrixExp3(VecToso3(theta_true));

% 计算传感器全局位置
sensor_positions = p_true + R_true * d_list;

num_sensors = size(sensor_positions, 2);
num_magnets = size(m_pos, 2);

%% 生成磁铁测量数据

% 初始化磁场存储
B_total = zeros(3, num_sensors);        % 全局坐标系磁场
b_total = zeros(3, num_sensors);        % 局部坐标系磁场
gradB_total = zeros(3, 3, num_sensors); % 全局坐标系梯度
gradb_total = zeros(3, 3, num_sensors); % 局部坐标系梯度

for sensor_idx = 1:num_sensors
    % 临时存储全局坐标系下的总和
    B_t = zeros(3, 1);
    gradB_t = zeros(3, 3);

    for magnet_idx = 1:num_magnets
        % 计算相对位置（全局坐标系）
        r = sensor_positions(:, sensor_idx) - m_pos(:, magnet_idx);

        % 获取当前磁铁参数
        moment_unit = m_hat(:, magnet_idx);
        moment_mag = m_norm(magnet_idx);

        % 计算磁场和梯度（全局坐标系）
        [B_single, gradB_single] = dipole_b_and_gradb(r, moment_unit, moment_mag);

        % 累加全局磁场和梯度
        B_t = B_t + B_single;
        gradB_t = gradB_t + gradB_single;
    end

    % 存储全局坐标系结果
    B_total(:, sensor_idx) = B_t;
    gradB_total(:, :, sensor_idx) = gradB_t;

    % 转换到局部坐标系并存储
    b_total(:, sensor_idx) = R_true' * B_t;
    gradb_total(:, :, sensor_idx) = R_true' * gradB_t * R_true;
end

%% 增加地磁干扰和噪声
% % 添加高斯噪声
B_total = B_total + randn(size(B_total)) * 1e-6;
b_total = b_total + randn(size(b_total)) * 1e-6;

% % 增加地磁干扰
B_total = B_total + [0; 0; 1e-6];
b_total = b_total + [0; 0; 1e-6];

%% 实验设置
num_experiments = 50;
results = struct();
% 优化器参数

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'off', ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6, ...
    'MaxIter', 1000, ...
    'MaxFunctionEvaluations', 10000);

% options2 = optimoptions('lsqnonlin', ...
%     'Algorithm', 'levenberg-marquardt', ...
%     'Display', 'off');

% 工作空间约束参数
workspace_center = [0; 0; 0];
workspace_radius = 0.5;

% 增加位置约束
lb_p = workspace_center - workspace_radius; % 下界
ub_p = workspace_center + workspace_radius; % 上界

for exp_idx = 1:num_experiments
    fprintf('\n===== 实验 %d/%d =====\n', exp_idx, num_experiments);

    % ==== 生成初始扰动 ====
    init_error = -1 + 2 * rand(3,1);
    p_init = p_true + 1e-3 * init_error;
    theta_init = theta_true + pi/2 * init_error;
    R_init = MatrixExp3(VecToso3(theta_init));

    % ==== 算法调用 ====
    % LM
    % t_lm_start = tic;
    % [p_lm, R_lm, stats_lm] = estimate_pose_lm( ...
    %     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
    % t_lm = toc(t_lm_start);

    % ELM
    % t_elm_start = tic;
    % [p_elm, R_elm, stats_elm] = estimate_pose_elm( ...
    %     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
    % t_elm = toc(t_elm_start);

    % 所提算法
    t_ours_start = tic;
    [p_ours, R_ours, stats_ours] = estimate_pose_ours( ...
        b_total, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p );
    t_ours = toc(t_ours_start);

    % Rlm
    t_Rlm_start = tic;
    [p_Rlm, R_Rlm, stats_Rlm] = estimate_R_lm( ...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );
    t_Rlm = toc(t_Rlm_start);

    % ==== 结果存储 ====
    % results(exp_idx).p_lm    = p_lm;
    % results(exp_idx).R_lm    = R_lm;
    % results(exp_idx).t_lm    = t_lm;
    % results(exp_idx).p_elm   = p_elm;
    % results(exp_idx).R_elm   = R_elm;
    % results(exp_idx).t_elm   = t_elm;
    results(exp_idx).p_Rlm   = p_Rlm;
    results(exp_idx).R_Rlm   = R_Rlm;
    results(exp_idx).t_Rlm   = t_Rlm;
    results(exp_idx).p_ours  = p_ours;
    results(exp_idx).R_ours  = R_ours;
    results(exp_idx).t_ours  = t_ours;

    % ==== 误差分析 ====

    % LM
    % results(exp_idx).lm_pos_error = norm(p_lm - p_true);
    % results(exp_idx).lm_rot_error = norm(R_lm - R_true, 'fro');
    % results(exp_idx).lm_field_error = calculate_field_errors(p_lm, R_lm, b_total, d_list, m_pos, m_hat, m_norm);

    % ELM
    % results(exp_idx).elm_pos_error = norm(p_elm - p_true);
    % results(exp_idx).elm_rot_error = norm(R_elm - R_true, 'fro');
    % results(exp_idx).elm_field_error = calculate_field_errors(p_elm, R_elm, b_total, d_list, m_pos, m_hat, m_norm);

    % Rlm
    results(exp_idx).Rlm_pos_error = norm(p_Rlm - p_true);
    results(exp_idx).Rlm_rot_error = norm(R_Rlm - R_true, 'fro');
    results(exp_idx).Rlm_field_error = calculate_field_errors(p_Rlm, R_Rlm, b_total, d_list, m_pos, m_hat, m_norm);

    % ==== 其它算法（如union）可按需补充 ====

    % 所提算法
    results(exp_idx).ours_pos_error = norm(p_ours - p_true);
    results(exp_idx).ours_rot_error = norm(R_ours - R_true, 'fro');
    results(exp_idx).ours_field_error = calculate_field_errors(p_ours, R_ours, b_total, d_list, m_pos, m_hat, m_norm);

end

% 统一可视化所有实验结果
visualize_pose(m_pos, m_hat, m_norm, p_true, R_true, results, d_list);

display_statistical_summary(results, num_experiments);

plot_error_distributions(results);



