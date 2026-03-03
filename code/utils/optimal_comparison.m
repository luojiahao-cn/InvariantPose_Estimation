clc,clear,close all
load('./exp/mat_data/optimized_params.mat', 'params');

% 重构后的绘图逻辑：始终绘制全部12个传感器作为背景，根据传入的 active_indices 高亮显示
% 所有的坐标都基于 d_list1 (12个传感器)
d_base = params.sensor.d_list; % 3×12
dc_base = mean(d_base, 2);
d_base = d_base - dc_base; % 统一中心化

% 配置1: 全部12个
idx1 = 1:12;

% 配置2: 
idx2 = [2,3,5,6,8,9,11,12];

% 配置3: 1-6个
idx4 = 1:6;

% 配置4:
idx3 = [3,6,9,12];

% 配置5:
idx5 = [5,6,3];

% 配置6:
idx6 = [2,6,9];

sensor_subsets = {idx1, idx2, idx3, idx4, idx5, idx6};
subset_labels = {'Redundant config. 1', 'Redundant config. 2', 'Redundant config. 3', 'Redundant config. 4', 'Minimal config. 1', 'Minimal config. 2'};

% 构建选择矩阵S (5个独立分量到9个矢量的映射)
S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
    0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
    0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

% 球面均匀采样
[theta, phi] = meshgrid(linspace(0, pi, 50), linspace(0, 2*pi, 100));
ux = sin(theta) .* cos(phi);
uy = sin(theta) .* sin(phi);
uz = cos(theta);

% 绘制传感器阵列与灵敏度图标
figure('Color', 'w', 'Units', 'centimeters', 'Position', [0, 0, 36, 12]);
t = tiledlayout(2, 6, 'TileSpacing', 'tight', 'Padding', 'tight');

% 第一行：绘制传感器布局
scalefactor = 2.5e-3;
for k = 1:6
    nexttile;
    active_idx = sensor_subsets{k};
    draw_sensor_new(d_base, active_idx, [0.7 0.7 0.7], 1e-3, 1e-3, 0.8e-3);
    draw_axes(eye(4), '', scalefactor); 
    axis equal;
    xlim([-2.5e-3 ,2.5e-3]);
    ylim([-2.5e-3 ,2.5e-3]);
    zlim([-2.5e-3 ,2.5e-3]);
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'LineWidth', 1, 'box', 'off');
    axis off
    view(125, 13);
    title(subset_labels{k}, 'Interpreter', 'latex', 'FontSize', 16);
end

% 第二行：绘制灵敏度图标 J(u)
for k = 1:6
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
    
    surf(J_vals.*ux, J_vals.*uy, J_vals.*uz, J_vals, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
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
    
    axis equal; axis off;
    view([1,1,1]);
end

% 设置共享 Colorbar
cb = colorbar;
cb.Layout.Tile = 'south';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;
cb.Label.String = 'Normalized Sensitivity $J(\mbox{\boldmath{$u$}})$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

export_fig('./results/sensitivityJ.png','-r600','-nocrop')
%% 辅助函数
function draw_sensor_new(d_list, active_indices, active_fc, L, W, H)
    % 绘制所有传感器，active_indices 指定的传感器使用 active_fc，其他的用浅色
    hold on;
    inactive_fc = [0.95, 0.95, 0.95]; % 浅灰色
    
    for i = 1:size(d_list, 2)
        center = d_list(:, i);
        if ismember(i, active_indices)
            % Active sensor
            draw_single_box(center, L, W, H, active_fc, 0.8); 
        else
            % Inactive sensor
            draw_single_box(center, L, W, H, inactive_fc, 0.1); % alpha = 0.1
        end
    end
end

function draw_single_box(center, L, W, H, fc, alpha)
        % 计算长方体的8个顶点（以中心为中心）
        cx = center(1); cy = center(2); cz = center(3);
        vertices = [
            cx - L/2, cy - W/2, cz - H/2;  % 1
            cx + L/2, cy - W/2, cz - H/2;  % 2
            cx + L/2, cy + W/2, cz - H/2;  % 3
            cx - L/2, cy + W/2, cz - H/2;  % 4
            cx - L/2, cy - W/2, cz + H/2;  % 5
            cx + W/2, cy - W/2, cz + H/2;  % 6
            cx + W/2, cy + W/2, cz + H/2;  % 7
            cx - W/2, cy + W/2, cz + H/2   % 8
        ];
        faces = [1, 2, 3, 4; 8, 7, 6, 5; 1, 5, 6, 2; 3, 7, 8, 4; 1, 4, 8, 5; 2, 6, 7, 3];
        patch('Faces', faces, 'Vertices', vertices, ...
            'FaceColor', fc, 'FaceAlpha', alpha, 'EdgeColor', 'k', 'LineWidth', 1);
end
