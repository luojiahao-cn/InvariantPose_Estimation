clc,clear,close all
load('./exp/mat_data/optimized_params.mat', 'params');

d_list = params.sensor.d_list(:, [1:12]);
% d_list = d_list - mean(d_list, 2);
num_sensors = size(d_list, 2);
pairs = nchoosek(1:num_sensors, 2);
D_delta = zeros(3, size(pairs,1));
% B_delta = zeros(3, size(pairs,1));

for idx = 1:size(pairs,1)
    i = pairs(idx, 1);
    j = pairs(idx, 2);
    D_delta(:, idx) = d_list(:, j) - d_list(:, i);
end

%% 估计局部梯度张量与球面采样可视化
% 构建选择矩阵S (5个独立分量到9个矢量的映射)
S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
    0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
    0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

C_matrix = kron(D_delta', eye(3)) * S;

% 预计算采样核心矩阵 (Information Matrix part)
P = D_delta * D_delta';
M = S * (inv(S' * kron(P, eye(3)) * S)) * S';

[V, D] = eig(P);

% 球面采样 (降低密度)
[theta, phi] = meshgrid(linspace(0, pi, 18), linspace(0, 2*pi, 36));
ux = sin(theta) .* cos(phi);
uy = sin(theta) .* sin(phi);
uz = cos(theta);

J_vals = zeros(size(theta));
for i = 1:numel(theta)
    u = [ux(i); uy(i); uz(i)];
    uu = kron(u, u);
    J_vals(i) = uu' * M * uu;
end

% 绘制“长度为J(u)”的图标 (使用线条)
X = J_vals .* ux;
Y = J_vals .* uy;
Z = J_vals .* uz;

figure('Color', 'w', 'Name', 'Optimal Sensitivity Icon');
quiver3(zeros(size(X)), zeros(size(Y)), zeros(size(Z)), X, Y, Z, 0, 'Color', 'b', 'LineWidth', 0.8);
hold on;
plot3(X, Y, Z, 'r.', 'MarkerSize', 8); % 添加末端点供参考
axis equal;
grid on;
view(45, 30);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(['Directional Sensitivity J(u) (Lines) for ' num2str(num_sensors) ' Sensors']);
hold off;
