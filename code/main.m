clc; clear; close all;

%% 磁铁参数定义
m_pos = [ 0 , 0.05;
          0 , 0.02;
          0 , 0.10]; % 磁铁位置 [m]

m_hat = [ 0 , 0.7;
          0 , 0.1;
          1 , 0.7];  % 磁化方向（未归一化）

m_hat = m_hat ./ vecnorm(m_hat); % 归一化磁化方向
m_norm = [8e2 , 8e2];            % 磁矩幅值 [A·m²]

d_list_e = [
    [0; 0; 0],...
    [1e-3; 0; 0],...
    [2e-3; 0; 0],...
    [0; 1e-3; 0],...
    [0; 0; 1e-3],...
    % [1e-3; 0; 1e-3]
];

row_means = mean(d_list_e, 2);
d_list = d_list_e - row_means;

theta_true = [0.1; 0.2; 0.3]; % 真实旋转向量 [rad]
p_true = [0.05; -0.03; 0.04]; % 传感器阵列参考点真实位置 [m]
R_true = MatrixExp3(VecToso3(theta_true));

% 计算传感器全局位置
sensor_positions = p_true + R_true * d_list;

num_sensors = size(sensor_positions, 2);
num_magnets = size(m_pos, 2);

% 设置随机种子
rng(2025);

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

%% 多次实验设置
num_experiments = 10; % 实验次数
results = struct();

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off');

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
    p_init = p_true + 0.05 * init_error;
    theta_init = theta_true + 0.1 * init_error;
    R_init = MatrixExp3(VecToso3(theta_init));

    % ==== 算法调用 ====
    % LM
    % [p_lm, R_lm, stats_lm] = estimate_pose_lm( ...
    %     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
    % ELM
    % [p_elm, R_elm, stats_elm] = estimate_pose_elm( ...
    %     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
    % 所提算法
    [p_ours, R_ours, stats_ours] = estimate_pose_ours( ...
        b_total, d_list, m_pos, m_hat, m_norm, p_init, R_true, p_true, options, lb_p, ub_p );
    % Rlm
    [p_Rlm, R_Rlm, stats_Rlm] = estimate_R_lm( ...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );

    % ==== 结果存储 ====
    % results(exp_idx).p_lm    = p_lm;
    % results(exp_idx).R_lm    = R_lm;
    % results(exp_idx).p_elm   = p_elm;
    % results(exp_idx).R_elm   = R_elm;
    results(exp_idx).p_ours  = p_ours;
    results(exp_idx).R_ours  = R_ours;
    results(exp_idx).p_Rlm   = p_Rlm;
    results(exp_idx).R_Rlm   = R_Rlm;

    % ==== 误差分析 ====
    % 初始误差
    results(exp_idx).init_pos_error = norm(p_init - p_true);
    results(exp_idx).init_rot_error = norm(R_init - R_true, 'fro');

    % LM
    % results(exp_idx).lm_pos_error = norm(p_lm - p_true);
    % results(exp_idx).lm_rot_error = norm(R_lm - R_true, 'fro');
    % results(exp_idx).lm_field_error = calculate_field_errors(p_lm, R_lm, b_total, d_list, m_pos, m_hat, m_norm);

    % ELM
    % results(exp_idx).elm_pos_error = norm(p_elm - p_true);
    % results(exp_idx).elm_rot_error = norm(R_elm - R_true, 'fro');
    % results(exp_idx).elm_field_error = calculate_field_errors(p_elm, R_elm, b_total, d_list, m_pos, m_hat, m_norm);

    % 所提算法
    results(exp_idx).ours_pos_error = norm(p_ours - p_true);
    results(exp_idx).ours_rot_error = norm(R_ours - R_true, 'fro');
    results(exp_idx).ours_field_error = calculate_field_errors(p_ours, R_ours, b_total, d_list, m_pos, m_hat, m_norm);

    % Rlm
    results(exp_idx).Rlm_pos_error = norm(p_Rlm - p_true);
    results(exp_idx).Rlm_rot_error = norm(R_Rlm - R_true, 'fro');
    results(exp_idx).Rlm_field_error = calculate_field_errors(p_Rlm, R_Rlm, b_total, d_list, m_pos, m_hat, m_norm);

    % ==== 其它算法（如union）可按需补充 ====
end

% 统一可视化所有实验结果
visualize_pose(m_pos, m_hat, m_norm, p_true, R_true, results, d_list);

display_statistical_summary(results, num_experiments);

plot_error_distributions(results);


