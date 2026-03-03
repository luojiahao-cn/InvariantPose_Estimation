figure(2)
close all
d0 = 1.5e-3;
d_list1 = [[0; 0; 0],...
    [d0; 0; 0],...
    [0; d0; 0],...
    [d0; d0; 0],...
    [0; 2*d0; 0],...
    [d0; 2*d0; 0],...
    [0; 0; -d0],...
    [d0; 0; -d0],...
    [0; d0; -d0],...
    [d0; d0; -d0],...
    [0; 2*d0; -d0],...
    [d0; 2*d0; -d0]]; % redundant config 1

% 计算几何中心
dc1 = mean(d_list1, 2);
d_list1 = d_list1 - dc1;

d_list2 = [[0; 0; 0],...
    [0; d0; 0],...
    [0; 0; -d0],...
    [0; d0; -d0]]; % redundant config 2
dc2 = mean(d_list2, 2);
d_list2 = d_list2 - dc2;

d_list3 = [[0; 0; 0],...
    [0; d0; 0],...
    [0; d0; -d0]]; % minimal config 3
dc3 = [0; d0/2; -d0/2];
d_list3 = d_list3 - dc3;

d_list4 = [[-d0; 0; 0],...
    [0; d0; 0],...
    [0; 0; d0]];
dc4 = [0; 0; 0];
d_list4 = d_list4 - dc4;

rankC1 = calcRankC(d_list1);
rankC2 = calcRankC(d_list2);
rankC3 = calcRankC(d_list3);
rankC4 = calcRankC(d_list4);
fprintf('rankC1: %d\n', rankC1);
fprintf('rankC2: %d\n', rankC2);
fprintf('rankC3: %d\n', rankC3);
fprintf('rankC4: %d\n', rankC4);


% 重构后的绘图逻辑：始终绘制全部12个传感器作为背景，根据传入的 active_indices 高亮显示
% 所有的坐标都基于 d_list1 (12个传感器)
load('./exp/mat_data/optimized_params.mat', 'params');
d_base = params.sensor.d_list; % 3×12
dc_base = mean(d_base, 2);
d_base = d_base - dc_base; % 统一中心化

% 配置1: 全部12个
idx1 = 1:12;

% 配置2: 
idx2 = [2,3,5,6,8,9,11,12];

% 配置3: 1-6个
idx3 = 1:6;

% 配置4:
idx4 = [3,6,9,12];

% 配置5:
idx5 = [5,6,3];

% 配置6:
idx6 = [2,6,9];

indices_cell = {idx1, idx2, idx3, idx4, idx5, idx6};
subset_labels = {'Sensors 1-12', 'Sensors 2,3,5,6,8,9,11,12', 'Sensors 1-6', 'Sensors 3,6,9,12', 'Sensors 5,6,3', 'Sensors 2,6,9'};

% 绘制传感器阵列
figure
tl = tiledlayout(2, 3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
set(gcf, 'color', 'w', 'units', 'centimeters', 'position', 2*[0, 0, 8.9, 6]);
scalefactor = 2.5e-3;

for k = 1:6
    nexttile
    active_idx = indices_cell{k};
    draw_sensor_new(d_base, active_idx, [0.7 0.7 0.7], 1e-3, 1e-3, 0.8e-3);
    % draw_axes(RpToTrans(rotx(-90), [0; 0; 0]), '', scalefactor); % 使用原来的绘制坐标轴
    draw_axes(eye(4), '', scalefactor); % 换回更标准的坐标轴
    axis equal;
    xlim([-2.5e-3 ,2.5e-3]);
    ylim([-2.5e-3 ,2.5e-3]);
    zlim([-2.5e-3 ,2.5e-3]);
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'LineWidth', 1, 'box', 'off'); % box off
    axis off
    view(115, 10);
    % view([0,0,1]);
    % title(subset_labels{k}, 'Interpreter', 'latex', 'FontSize', 14);
end

%%
export_fig('./results/magneto_sensor.png', '-r600');

function draw_sensor(d_list, fc, L, W, H)
    % 遍历d_list中的每一列，以该位置为几何中心绘制传感器
    hold on;
    for i = 1:size(d_list, 2)
        center = d_list(:, i);  % 传感器中心位置
        draw_single_box(center, L, W, H, fc, 0.5);
    end
end

function draw_sensor_new(d_list, active_indices, active_fc, L, W, H)
    % 绘制所有传感器，active_indices 指定的传感器使用 active_fc，其他的用浅色
    hold on;
    inactive_fc = [0.95, 0.95, 0.95]; % 浅灰色
    
    for i = 1:size(d_list, 2)
        center = d_list(:, i);
        if ismember(i, active_indices)
            % Active sensor
            % 使用传入的 active_fc 并设置较高的不透明度
            draw_single_box(center, L, W, H, active_fc, 0.8); 
        else
            % Inactive sensor
            % 设置较低的不透明度
            draw_single_box(center, L, W, H, inactive_fc, 0.1);
        end
    end
end

function draw_single_box(center, L, W, H, fc, alpha)
        % 计算长方体的8个顶点（以中心为中心）
        % 顶点定义：8×3矩阵，每行是一个顶点的[x, y, z]坐标
        cx = center(1);
        cy = center(2);
        cz = center(3);
    
        vertices = [
            cx - L/2, cy - W/2, cz - H/2;  % 1: 左下前
            cx + L/2, cy - W/2, cz - H/2;  % 2: 右下前
            cx + L/2, cy + W/2, cz - H/2;  % 3: 右上前
            cx - L/2, cy + W/2, cz - H/2;  % 4: 左上前
            cx - L/2, cy - W/2, cz + H/2;  % 5: 左下后
            cx + W/2, cy - W/2, cz + H/2;  % 6: 右下后
            cx + W/2, cy + W/2, cz + H/2;  % 7: 右上后
            cx - W/2, cy + W/2, cz + H/2   % 8: 左上后
        ];
        
        % 定义6个面的顶点索引（每个面4个顶点，按逆时针顺序）
        faces = [
            1, 2, 3, 4;  % 底面（z = -H/2）
            8, 7, 6, 5;  % 顶面（z = H/2）
            1, 5, 6, 2;  % 前面（y = -W/2）
            3, 7, 8, 4;  % 后面（y = W/2）
            1, 4, 8, 5;  % 左面（x = -W/2）
            2, 6, 7, 3   % 右面（x = W/2）
        ];

        % 使用patch绘制长方体
        patch('Faces', faces, 'Vertices', vertices, ...
            'FaceColor', fc, ...
            'FaceAlpha', alpha, ...
            'EdgeColor', 'k', ...
            'LineWidth', 1);
end

function rankC = calcRankC(d_list)
    pairs = nchoosek(1:size(d_list, 2), 2);
    D_matrix = zeros(3, size(pairs,1));

    for idx = 1:size(pairs,1)
        i = pairs(idx, 1);
        j = pairs(idx, 2);
        D_matrix(:, idx) = d_list(:, j) - d_list(:, i);
    end
    %% 估计局部梯度张量
    % 构建选择矩阵S
    S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
        0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
        0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

    % 构建完整约束矩阵C
    C_matrix = kron(D_matrix', eye(3)) * S;
    rankC = rank(C_matrix);
end