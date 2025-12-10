function [b_total, B_total, gradb_total, gradB_total, sensor_positions] = ...
    generate_magnetic_data(params)
% GENERATE_MAGNETIC_DATA 生成磁铁测量数据
% 输入：
%   params - 实验参数结构体
% 输出：
%   b_total      - 3×N传感器测量矩阵（局部坐标系），单位：T
%   B_total      - 3×N传感器测量矩阵（全局坐标系），单位：T
%   gradb_total  - 3×3×N局部坐标系梯度张量
%   gradB_total  - 3×3×N全局坐标系梯度张量
%   sensor_positions - 3×N传感器全局位置矩阵

% 提取参数
m_pos = params.magnet.m_pos;
m_hat = params.magnet.m_hat;
m_norm = params.magnet.m_norm;
d_list = params.sensor.d_list;
theta_true = params.ground_truth.theta_true;
p_true = params.ground_truth.p_true;

% 计算真实旋转矩阵和传感器全局位置
R_true = MatrixExp3(VecToso3(theta_true));
sensor_positions = p_true + R_true * d_list;

num_sensors = size(d_list, 2);

% 初始化磁场存储
B_total = zeros(3, num_sensors);        % 全局坐标系磁场
b_total = zeros(3, num_sensors);        % 局部坐标系磁场
gradB_total = zeros(3, 3, num_sensors); % 全局坐标系梯度
gradb_total = zeros(3, 3, num_sensors); % 局部坐标系梯度

% 计算每个传感器的磁场和梯度
for sensor_idx = 1:num_sensors
    % 临时存储全局坐标系下的总和

    [B_t, gradB_t] = calcFieldAndGradient(sensor_positions(:, sensor_idx), m_pos, m_hat, m_norm);

    % 存储全局坐标系结果
    B_total(:, sensor_idx) = B_t;
    gradB_total(:, :, sensor_idx) = gradB_t;

    % 转换到局部坐标系并存储
    b_total(:, sensor_idx) = R_true' * B_t;
    gradb_total(:, :, sensor_idx) = R_true' * gradB_t * R_true;
end
end

