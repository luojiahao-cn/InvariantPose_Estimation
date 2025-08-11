function field_error = calculate_field_errors(p_est, R_est, b_total_actual, d_list, m_pos, m_hat, m_norm)
% CALCULATE_FIELD_ERRORS 计算估计位置和旋转与实际磁场测量值之间的误差
%
% 输入参数：
%   p_est - 估计的位置 [3×1]
%   R_est - 估计的旋转矩阵 [3×3]
%   b_total_actual - 实际测量的磁场数据 [3×N]
%   d_list - 传感器在参考坐标系中的偏移 [3×N]
%   m_pos - 磁铁位置矩阵 [3×K]
%   m_hat - 磁化方向单位向量 [3×K]
%   m_norm - 磁矩幅值向量 [1×K]
%
% 输出参数：
%   field_error - 磁场误差总和 (标量)

num_sensors = size(d_list, 2);
num_magnets = size(m_pos, 2);

% 初始化误差存储
field_error_total = 0;

% 计算每个传感器的预测值和误差
for sensor_idx = 1:num_sensors
    % 计算传感器在全局坐标系中的位置
    sensor_pos_global = p_est + R_est * d_list(:, sensor_idx);
    
    % 初始化预测磁场
    B_pred = zeros(3, 1);
    
    % 计算所有磁铁对该传感器的影响
    for magnet_idx = 1:num_magnets
        r = sensor_pos_global - m_pos(:, magnet_idx);
        B = dipole_b_and_gradb(r, m_hat(:, magnet_idx), m_norm(magnet_idx));
        B_pred = B_pred + B;
    end
    
    % 将预测值转换到局部坐标系
    b_pred = R_est' * B_pred;
    
    % 计算磁场误差
    field_error = norm(b_pred - b_total_actual(:, sensor_idx));
    field_error_total = field_error_total + field_error;
end

field_error = field_error_total;
end