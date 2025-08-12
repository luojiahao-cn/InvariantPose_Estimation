function [p_opt, R_opt, stats] = estimate_R_lm(b_total, d_list, m_pos, m_hat, m_norm, omega_init, p_fixed, options)
% ESTIMATE_R_LM 固定位置p，优化旋转矩阵R
%
% 输入参数：
%   b_total    - 3×N传感器测量矩阵（局部坐标系），单位：T
%   d_list     - 3×N位移矩阵，传感器在参考坐标系中的偏移，单位：m
%   m_pos      - 3×K磁铁位置矩阵，单位：m
%   m_hat      - 3×K磁化方向单位向量（归一化）
%   m_norm     - 1×K磁矩幅值向量，单位：A·m²
%   omega_init - 3×1初始旋转向量估计 [rad]
%   p_fixed    - 3×1固定位置向量 [m]
%   options    - （可选）优化选项
%
% 输出参数：
%   p_opt      - 3×1估计的位置向量
%   R_opt      - 3×3估计的旋转矩阵
%   stats      - 优化统计信息

p_opt = p_fixed;

if nargin < 8
    options = [];
end

% 设置初始变量（仅旋转向量）
x0 = omega_init;

% 设置变量上下界
lb = -inf(3,1);
ub = inf(3,1);

num_sensors = size(b_total, 2);
if num_sensors ~= size(d_list, 2)
    error('传感器测量数量(%d)和偏移数量(%d)不匹配', num_sensors, size(d_list, 2));
end

if isempty(options)
    options = optimoptions('lsqnonlin', ...
        'Algorithm', 'levenberg-marquardt', ...
        'Display', 'off');
end

history.iter = [];
history.resnorm = [];

[x_opt, resnorm, residual, exitflag, output] = lsqnonlin( ...
    @(x)lm_objective(x, m_pos, m_hat, m_norm, d_list, b_total, p_fixed), ...
    x0, lb, ub, options);

R_opt = MatrixExp3(VecToso3(x_opt));

stats.resnorm = resnorm;
stats.residual = residual;
stats.exitflag = exitflag;
stats.output = output;
stats.history = history;

% ===== 内部函数：优化目标函数 =====
function residuals = lm_objective(x, m_pos, m_hat, m_norm, d_list, b_list, p)
    omega = x(1:3);
    R = MatrixExp3(VecToso3(omega));
    num_sensors = size(d_list, 2);
    residuals = zeros(3*num_sensors, 1);
    for i = 1:num_sensors
        p_sensor_i = p + R * d_list(:, i);
        B_global = zeros(3,1);
        num_magnets = size(m_pos,2);
        for j = 1:num_magnets
            r = p_sensor_i - m_pos(:, j);
            [B_single, ~] = dipole_b_and_gradb(r, m_hat(:,j), m_norm(j));
            B_global = B_global + B_single;
        end
        b_pred = R' * B_global;
        b_meas = b_list(:, i);
        residuals(3*(i-1)+1:3*i) = b_pred - b_meas;
    end
end

end
