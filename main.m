clc; clear; close all;
%% 磁铁参数定义
m_pos = [ 0 , 0.05;
          0 , 0.02;
          0 , 0.10]; % 磁铁位置 [m]

m_hat = [ 0 , 0.7;
          0 , 0.1;
          1 , 0.7];  % 磁化方向（未归一化）

m_hat = m_hat ./ vecnorm(m_hat); % 归一化磁化方向
m_norm = [1e3 , 8e2];            % 磁矩幅值 [A·m²]

d_list_e = [   0,  1e-3, 2e-3,     0;
               0,     0,    0,  1e-3;
               0,     0,    0,     0]; % 传感器偏移 [m]
row_means = mean(d_list_e, 2);
d_list = d_list_e - row_means;


theta_true = [0.1; 0.2; 0.3]; % 真实旋转向量 [rad]
p_true = [0.05; -0.03; 0.04];    % 传感器阵列参考点真实位置 [m]
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

%% 多次实验设置
num_experiments = 50; % 实验次数
results = struct();

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt',...
    'Display', 'off');

for exp_idx = 1:num_experiments
    fprintf('\n===== 实验 %d/%d =====', exp_idx, num_experiments);
    
    % 生成随机初始误差
    init_error =  -1 + 2 * rand(3,1);
    % init_error = 0;
    
    % 初始位置和旋转（添加扰动）
    p_init = p_true + 0.03 * init_error; 
    % theta_init = theta_true + 1 * init_error;
    theta_init = 3 * init_error;
    R_init = MatrixExp3(VecToso3(theta_init));

    % 调用LM算法
    [p_lm, R_lm, stats_lm] = estimate_pose_lm(...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options);

    % 调用ELM算法
    [p_elm, R_elm, stats_elm] = estimate_pose_elm(...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options);

    % 调用融合算法
    [p_union, R_union, stats_union] = estimate_pose_union(...
        b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options);

    % 调用所提算法
    [p_prop, R_prop, stats_prop] = estimate_pose_ours(...
        b_total, d_list, m_pos, m_hat, m_norm, p_init, R_true, p_true, options);
    
    % 计算初始误差
    init_pos_error = norm(p_init - p_true);
    init_rot_error = norm(R_init - R_true,'fro');
    
    % 计算LM算法结果误差
    lm_pos_error = norm(p_lm - p_true);
    lm_rot_error = norm(R_lm - R_true,'fro');

    % 计算ELM算法结果误差
    elm_pos_error = norm(p_elm - p_true);
    elm_rot_error = norm(R_elm - R_true,'fro');

    % 计算融合算法结果误差
    union_pos_error = norm(p_union - p_true);
    union_rot_error = norm(R_union - R_true,'fro');

    % 计算所提算法结果误差
    prop_pos_error = norm(p_prop - p_true);
    prop_rot_error = norm(R_prop - R_true,'fro');

    % 存储结果
    results(exp_idx).init_pos_error = init_pos_error;
    results(exp_idx).init_rot_error = init_rot_error;
    results(exp_idx).lm_pos_error = lm_pos_error;
    results(exp_idx).lm_rot_error = lm_rot_error;
    results(exp_idx).elm_pos_error = elm_pos_error;
    results(exp_idx).elm_rot_error = elm_rot_error;
    results(exp_idx).union_pos_error = union_pos_error;
    results(exp_idx).union_rot_error = union_rot_error;
    results(exp_idx).prop_pos_error = prop_pos_error;
    results(exp_idx).prop_rot_error = prop_rot_error;
end

%% 调用文本输出函数
display_statistical_summary(results, num_experiments);

%% 调用图像输出函数
plot_error_distributions(results);
