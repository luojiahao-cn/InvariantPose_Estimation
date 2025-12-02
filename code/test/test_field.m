clc,clear,close all;
%% 测试磁场强度满足+- 32Gs

% 设置两个偶极子的朝向和位置

nu1 = [-0.2588         0   0.9659   -0.1125    0.0002    0.3099]';
nu2 = [0.2588         0   0.9659    0.1125    0.0002    0.3099]';

% figure; hold on; axis equal; grid on;
% % 绘制第一个磁偶极子的朝向和位置
% quiver3(nu1(4), nu1(5), nu1(6), nu1(1), nu1(2), nu1(3), 0.05, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% % 绘制第二个磁偶极子的朝向和位置
% quiver3(nu2(4), nu2(5), nu2(6), nu2(1), nu2(2), nu2(3), 0.05, 'b', 'LineWidth', 2, 'MaxHeadSize', 2);
% xlabel('x'); ylabel('y'); zlabel('z');
% legend('Dipole 1', 'Dipole 2');
% title('磁偶极子的朝向与位置');

%% 计算磁场
% 在R = 15e-3的半球内随机采样1000个点并计算磁场
R = 15e-3;
num_points = 1000;
points = rand(3, num_points);
points = points ./ vecnorm(points);
points = R * points;

% 计算这些点的平均磁场模值
b_total = zeros(3, num_points);
nabla_b_total = zeros(3, 3, num_points);
k = 300;
for i = 1:num_points
    p = points(:, i);
    [b1, nabla_b1] = dipole_magnetic_field(nu1, p, k);
    [b2, nabla_b2] = dipole_magnetic_field(nu2, p, k);
    b = b1 + b2;
    b_total(:, i) = b;
    nabla_b_total(:, :, i) = nabla_b1 + nabla_b2;
end
b_total = mean(b_total, 2);
nabla_b_total = mean(nabla_b_total, 3);

figure; hold on; axis equal; grid on;
plot(b_total);
plot(nabla_b_total(3, :, :));
xlabel('x'); ylabel('y'); zlabel('z');
legend('b', 'nabla_b');
title('磁场和梯度');

%% 计算磁场模值
b_norm = vecnorm(b_total);
b_norm = mean(b_norm);

disp(['磁场模值: ', num2str(b_norm)]);