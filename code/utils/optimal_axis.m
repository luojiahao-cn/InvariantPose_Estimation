clc,clear,close all
load('./exp/mat_data/optimized_params.mat', 'params');

% 定义要对比的传感器索引子集
sensor_subsets = {1:12, [2,3,5,6,8,9,11,12], 1:6, [2,3,11,12], [5,6,3], [2,6,9]};
subset_labels = {'Sensors 1-12', 'Sensors 2,3,5,6,8,9,11,12', 'Sensors 1-6', 'Sensors 2,3,11,12', 'Sensors 5,6,3', 'Sensors 2,6,9'};

% 构建基础矩阵与网格
S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
    0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
    0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

[theta, phi] = meshgrid(linspace(0, pi, 50), linspace(0, 2*pi, 100));
ux = sin(theta) .* cos(phi);
uy = sin(theta) .* sin(phi);
uz = cos(theta);

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.2, 0.8, 0.6]);
set(gcf, 'color', 'w', 'units', 'centimeters', 'position', [0, 0, 19.8, 12]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:length(sensor_subsets)
    idx_list = sensor_subsets{k};
    d_list = params.sensor.d_list(:, idx_list);
    num_sensors = size(d_list, 2);
    
    % 计算全组合微分向量
    pairs = nchoosek(1:num_sensors, 2);
    D_delta = zeros(3, size(pairs, 1));
    for idx = 1:size(pairs, 1)
        D_delta(:, idx) = d_list(:, pairs(idx, 2)) - d_list(:, pairs(idx, 1));
    end
    
    % 预计算核心矩阵
    P = D_delta * D_delta';
    M = S * (inv(S' * kron(P, eye(3)) * S)) * S';
    
    % 计算球面灵敏度
    J_vals = zeros(size(theta));
    for i = 1:numel(theta)
        u = [ux(i); uy(i); uz(i)];
        uu = kron(u, u);
        J_vals(i) = uu' * M * uu;
    end
    J_vals = J_vals ./ max(J_vals(:)); % 归一化到 1
    
    % 绘制图标
    nexttile;
    X = J_vals .* ux;
    Y = J_vals .* uy;
    Z = J_vals .* uz;
    
    surf(X, Y, Z, J_vals, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on
    draw_axes(eye(4), '', 1.5);
    shading interp;
    
    % 上色方案
    try
        colormap(colorwarm);
    catch
        c1 = [0.2, 0.2, 0.8]; c2 = [0.9, 0.9, 0.9]; c3 = [0.8, 0.2, 0.2];
        colormap([interp1([0, 0.5, 1], [c1; c2; c3], linspace(0, 1, 256))]);
    end
    
    % axis equal; grid on; view(45, 30);
    % title(subset_labels{k}, 'Interpreter', 'latex', 'FontSize', 14);
    % xlabel('$x$', 'Interpreter', 'latex');
    % ylabel('$y$', 'Interpreter', 'latex');
    % zlabel('$z$', 'Interpreter', 'latex');
    % xlim([-1.2, 1.2]); ylim([-1.2, 1.2]); zlim([-1.2, 1.2]);
    % set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
    axis equal
    axis off
    view([1,1,1])
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 16;
ylabel(cb, 'Normalized $J(\boldsymbol{u})$', 'Interpreter', 'latex');
