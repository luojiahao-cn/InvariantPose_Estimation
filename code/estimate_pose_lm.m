function [p_opt, R_opt, stats] = estimate_pose_lm(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options)
% ESTIMATE_POSE_LM 使用Levenberg-Marquardt算法估计传感器姿态
%
% 输入参数：
%   b_total    - 3×N传感器测量矩阵（局部坐标系），单位：T
%   d_list     - 3×N位移矩阵，传感器在参考坐标系中的偏移，单位：m
%   m_pos      - 3×K磁铁位置矩阵，单位：m
%   m_hat      - 3×K磁化方向单位向量（归一化）
%   m_norm     - 1×K磁矩幅值向量，单位：A·m²
%   theta_init - 3×1初始旋转向量估计 [rad]
%   p_init     - 3×1初始位置估计 [m]
%   options    - （可选）优化选项
%
% 输出参数：
%   R_opt      - 3×3估计的旋转矩阵（传感器到全局坐标系的旋转）
%   p_opt      - 3×1估计的位置向量 [m]
%   stats      - 包含优化统计信息的结构体

% 参数检查
if nargin < 8
    options = [];
end

% 组合初始变量
x0 = [theta_init; p_init];

% 检查传感器数量
num_sensors = size(b_total, 2);
if num_sensors ~= size(d_list, 2)
    error('传感器测量数量(%d)和偏移数量(%d)不匹配', num_sensors, size(d_list, 2));
end

% 设置默认优化选项
if isempty(options)
    options = optimoptions('lsqnonlin', ...
        'Algorithm', 'levenberg-marquardt', ...
        'Display', 'off');
end

% 清空历史记录
history.iter = [];
history.resnorm = [];

% 优化计算
[x_opt, resnorm, residual, exitflag, output] = lsqnonlin(...
    @(x)lm_objective(x, m_pos, m_hat, m_norm, d_list, b_total), ...
    x0, [], [], options);

% 提取结果
omega_opt = x_opt(1:3);
p_opt = x_opt(4:6);
R_opt = MatrixExp3(VecToso3(omega_opt));

% 保存统计信息
stats.resnorm = resnorm;
stats.residual = residual;
stats.exitflag = exitflag;
stats.output = output;
stats.history = history;

% ===== 内部函数：优化目标函数 =====
function residuals = lm_objective(x, m_pos, m_hat, m_norm, d_list, b_list)
    % 提取参数
    omega = x(1:3);
    p = x(4:6);
    
    % 将旋转向量转换为旋转矩阵
    R = MatrixExp3(VecToso3(omega));
    
    num_sensors = size(d_list, 2);
    residuals = zeros(3*num_sensors, 1); % 每个传感器3个分量
    
    % 对每个传感器
    for i = 1:num_sensors
        % 传感器全局位置
        p_sensor_i = p + R * d_list(:, i);
        
        % 计算该位置的总磁场（全局坐标系）
        B_global = zeros(3,1);
        num_magnets = size(m_pos,2);
        for j = 1:num_magnets
            r = p_sensor_i - m_pos(:, j);
            [B_single, ~] = dipole_b_and_gradb(r, m_hat(:,j), m_norm(j));
            B_global = B_global + B_single;
        end
        
        % 将全局磁场转换到传感器局部坐标系
        b_pred = R' * B_global;
        
        % 实际测量值（局部坐标系）
        b_meas = b_list(:, i);
        
        % 残差
        residuals(3*(i-1)+1:3*i) = b_pred - b_meas;
    end
end

end % 主函数结束