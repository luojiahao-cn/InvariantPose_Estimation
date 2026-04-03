function [b_total, sensor_positions] = generate_magnetic_data(params)
% GENERATE_MAGNETIC_DATA Generate magnetic measurement data
% Input:
%   params - Experiment parameter structure
% Output:
%   b_total      - 3xN sensor measurement matrix (local frame)
%   sensor_positions - 3xN sensor global position matrix

m_pos = params.magnet.m_pos;
m_hat = params.magnet.m_hat;
m_norm = params.magnet.m_norm;
d_list = params.sensor.d_list;
theta_true = params.ground_truth.theta_true;
p_true = params.ground_truth.p_true;

R_true = MatrixExp3(VecToso3(theta_true));
sensor_positions = p_true + R_true * d_list;

num_sensors = size(d_list, 2);
b_total = zeros(3, num_sensors);

for sensor_idx = 1:num_sensors
    [B_t, ~] = calcFieldAndGradient(sensor_positions(:, sensor_idx), m_pos, m_hat, m_norm);
    b_total(:, sensor_idx) = R_true' * B_t;
end
end
