function [p_est, R_est, stats] = estimate_pose_elm(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options)
% ESTIMATE_POSE_ELM 通过优化传感器间的磁场差估计姿态（同时优化位置和旋转）
%
% 输入参数：
%   b_total    - 3×N传感器测量矩阵（局部坐标系），单位：T
%   d_list     - 3×N位移矩阵，传感器在参考坐标系中的偏移，单位：m
%   m_pos      - 3×K磁铁位置矩阵，单位：m
%   m_hat      - 3×K磁化方向单位向量（归一化）
%   m_norm     - 1×K磁矩幅值向量，单位：A·m²
%   theta_init - 3×1初始旋转向量估计 [rad]
%   p_init     - 3×1初始位置估计 [m]
%
% 输出参数：
%   p_est    - 3×1估计的位置向量 [m]
%   R_est    - 3×3估计的旋转矩阵
%   stats    - 包含优化统计信息的结构体

% 参数检查
num_sensors = size(b_total, 2);
if num_sensors ~= size(d_list, 2)
    error('传感器测量数量(%d)和偏移数量(%d)不匹配', num_sensors, size(d_list, 2));
end

% 获取传感器对组合
pairs = nchoosek(1:num_sensors, 2);
num_pairs = size(pairs, 1);

% 构建测量差矩阵
delta_b_meas = zeros(3, num_pairs);
delta_d = zeros(3, num_pairs);

for idx = 1:num_pairs
    i = pairs(idx, 1);
    j = pairs(idx, 2);
    delta_b_meas(:, idx) = b_total(:, j) - b_total(:, i);
    delta_d(:, idx) = d_list(:, j) - d_list(:, i);
end

% 设置优化选项
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off');

% 组合初始变量（旋转向量和位置）
x0 = [theta_init; p_init];

% 运行优化
[x_opt, resnorm, residual, exitflag, output] = lsqnonlin(...
    @(x)elm_objective(x, m_pos, m_hat, m_norm, d_list, delta_d, delta_b_meas), ...
    x0, [], [], options);

% 提取结果
theta_opt = x_opt(1:3);
p_est = x_opt(4:6);
R_est = MatrixExp3(VecToso3(theta_opt));

% 保存统计信息
stats.resnorm = resnorm;
stats.residual = residual;
stats.exitflag = exitflag;
stats.output = output;

% ===== 内部目标函数 =====
function residuals = elm_objective(x, m_pos, m_hat, m_norm, d_list, delta_d, delta_b_meas)
    % 提取参数
    theta = x(1:3);
    p = x(4:6);
    
    % 将旋转向量转换为旋转矩阵
    R = MatrixExp3(VecToso3(theta));
    
    num_pairs = size(delta_d, 2);
    residuals = zeros(3*num_pairs, 1);
    
    % 对每个传感器对
    for idx = 1:num_pairs
        i = pairs(idx, 1);
        j = pairs(idx, 2);
        
        % 计算传感器i和j的全局位置
        p_i = p + R * d_list(:, i);
        p_j = p + R * d_list(:, j);
        
        % 计算在p_i处的磁场（全局坐标系）
        B_i = zeros(3,1);
        for k = 1:size(m_pos,2)
            r_i = p_i - m_pos(:, k);
            [B_single, ~] = dipole_b_and_gradb(r_i, m_hat(:,k), m_norm(k));
            B_i = B_i + B_single;
        end
        
        % 计算在p_j处的磁场（全局坐标系）
        B_j = zeros(3,1);
        for k = 1:size(m_pos,2)
            r_j = p_j - m_pos(:, k);
            [B_single, ~] = dipole_b_and_gradb(r_j, m_hat(:,k), m_norm(k));
            B_j = B_j + B_single;
        end
        
        % 转换到局部坐标系
        b_i = R' * B_i;
        b_j = R' * B_j;
        
        % 计算预测的磁场差
        delta_b_pred = b_j - b_i;
        
        % 残差
        residuals(3*(idx-1)+1:3*idx) = delta_b_pred - delta_b_meas(:, idx);
    end
end

end