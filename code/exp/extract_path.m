clc,clear,close all

% 1. 加载数据
load('../exp/mat_data/optimized_params_1.mat', 'test_points');
p_all = cat(2, test_points.p_true);

p_ind = 1:1000;

for i = 1:2:9
    range = (i*100 + 1) : (i+1)*100;
    p_ind(range) = fliplr(p_ind(range));
end

% 蛇形索引处理：从第11个开始，每隔10个翻转一次 (11-20, 31-40, ...)
for i = 1:2:99 
    range = (i*10 + 1) : (i+1)*10;
    p_ind(range) = fliplr(p_ind(range));
end

p_all = p_all(:, p_ind); % 重新排序点

% --- 3. 定义嵌套立方体路径 (4层) ---
layers_lo = [5, 4, 3, 2]; % 每层最小索引 (内到外)
layers_hi = [6, 7, 8, 9]; % 每层最大索引 (内到外)
x_coords = linspace(-0.1, 0.1, 10); 
grid_sz = [10, 10, 10];

all_grid_path = [];
cubic_path_ideal = [];

for L = 1:length(layers_lo)
    i_min = layers_lo(L); i_max = layers_hi(L);
    j_min = i_min; j_max = i_max;
    k_min = i_min; k_max = i_max;
    
    % 定义当前层立方体的 8 个顶点顺序
    nodes = [
        i_min, j_min, k_min; % V1
        i_min, j_max, k_min; % V2
        i_min, j_max, k_max; % V3
        i_min, j_min, k_max; % V4
        i_max, j_min, k_max; % V5
        i_max, j_max, k_max; % V6
        i_max, j_max, k_min; % V7
        i_max, j_min, k_min; % V8
        i_min, j_min, k_min  % 回到层起点
    ];
    
    % 记录理想物理轨迹 (用于绘图)
    cubic_path_ideal = [cubic_path_ideal, x_coords(nodes)'];
    
    % 4. 逐段生成当前层的格点
    for s = 1:size(nodes, 1)-1
        p1 = nodes(s, :); p2 = nodes(s+1, :);
        diff_dim = find(p1 ~= p2);
        if isempty(diff_dim)
            segment = p1';
        else
            v1 = p1(diff_dim); v2 = p2(diff_dim);
            step = sign(v2 - v1);
            vals = v1:step:v2;
            segment = repmat(p1', 1, length(vals));
            segment(diff_dim, :) = vals;
        end
        all_grid_path = [all_grid_path, segment];
    end
    
    % --- 重新应用：层间直角连接桥 (Axial Bridge) ---
    if L < length(layers_lo)
        p_current_end = nodes(end, :); 
        p_next_start_val = layers_lo(L+1);
        
        % 生成直角转弯的三个关键支撑点
        bridge_steps = [
            p_next_start_val, p_current_end(2), p_current_end(3); % 移到下一层的 X
            p_next_start_val, p_next_start_val, p_current_end(3); % 移到下一层的 Y
            p_next_start_val, p_next_start_val, p_next_start_val  % 移到下一层的 Z
        ];
        
        % 将桥接点加入逻辑路径
        all_grid_path = [all_grid_path, bridge_steps'];
        
        % 将桥接点加入理想物理轨迹 (不放 NaN，使其在视觉上连通)
        bridge_phys = x_coords(bridge_steps)';
        cubic_path_ideal = [cubic_path_ideal, bridge_phys];
    end
end

% 5. 转换为线性索引并去重
path_indices = sub2ind(grid_sz, all_grid_path(1,:), all_grid_path(2,:), all_grid_path(3,:));
path_indices = unique(path_indices, 'stable'); 

% 6. 还原回原始与纠偏后的索引
cubic_path_raw = p_ind(path_indices)'; 
cubic_path = path_indices'; 

% --- 理想路径已在循环中生成 ---

save('./mat_data/path_ind.mat', 'cubic_path_raw', 'cubic_path', 'cubic_path_ideal', '-append');
fprintf('生成成功！数据已保存：\n - cubic_path: 矫正后的几何索引\n - cubic_path_ideal: 理想物理路径 (8点顶点顺序)\n');

% 7. 可视化：按顺序绘制点并标记轨迹
figure('Color', 'w', 'Position', [100 100 800 600]);
% 绘制参考背景
scatter3(p_all(1,:), p_all(2,:), p_all(3,:), 15, [0.85 0.85 0.85], 'filled', 'MarkerEdgeAlpha', 0.1);
hold on;

% 提取坐标
pts = p_all(:, path_indices);

% 绘制主轨迹线
plot3(pts(1,:), pts(2,:), pts(3,:), 'b-', 'LineWidth', 1.5);

% 绘制带编号的点
colors = jet(length(path_indices));
for i = 1:length(path_indices)
    plot3(pts(1,i), pts(2,i), pts(3,i), 'o', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'k');
    % 每隔几个点标一下序号，避免太拥挤
    if mod(i, 5) == 1 || i == length(path_indices)
        text(pts(1,i), pts(2,i), pts(3,i), sprintf('  #%d', i), ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end

% 特别标记起点和终点
plot3(pts(1,1), pts(2,1), pts(3,1), 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'START');
plot3(pts(1,end), pts(2,end), pts(3,end), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'END');

grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Ordered Path Trajectory (Vertex Sequence Trace)');
view(125, 25);
% colorbar('Ticks', [0, 1], 'TickLabels', {'Start', 'End'}, 'Label', 'Path Progress');

%% --- 新增路径 A：3D 螺旋轨迹 (Spiral Path) ---
t = linspace(0, 4*pi, 60); % 绕两圈
x_sp = round(5 + 4*cos(t)); % 中心在5, 半径4
y_sp = round(5 + 4*sin(t));
z_sp = round(linspace(1, 8, 60)); % 从1层平滑上升到8层

% 理想连续螺旋路径
t_fine = linspace(0, 4*pi, 200);
spiral_path_ideal = [
    -0.1 + (5 + 4*cos(t_fine) - 1) * (0.2/9);
    -0.1 + (5 + 4*sin(t_fine) - 1) * (0.2/9);
    -0.1 + (linspace(1, 8, 200) - 1) * (0.2/9)
];

path_indices_sp = sub2ind(grid_sz, x_sp, y_sp, z_sp);
path_indices_sp = unique(path_indices_sp, 'stable');
spiral_path = path_indices_sp';
spiral_path_raw = p_ind(path_indices_sp)';

save('path_ind.mat', 'spiral_path', 'spiral_path_raw', 'spiral_path_ideal', '-append');
fprintf('Spiral Path 生成成功！\n');

%% --- 新增路径 B：3D 空间星形轨迹 (Star Path) ---
% 选取 8 个极值顶点，按照交叉步进的方式排列，形成类似 3D 星形的连线
i_lo = 2; i_hi = 9;
j_lo = 2; j_hi = 9;
k_lo = 1; k_hi = 8;

star_nodes = [
    5, 5, k_hi;    % 顶点
    i_lo, j_lo, k_lo;
    i_hi, j_hi, k_lo;
    i_lo, j_hi, k_hi;
    i_hi, j_lo, k_hi;
    5, 5, k_lo;    % 底心
    5, 5, k_hi     % 回到顶点
];

% 理想星形路径 (折线)
star_path_ideal = x_coords(star_nodes)';

grid_path_star = [];
for s = 1:size(star_nodes, 1)-1
    p1 = star_nodes(s, :); p2 = star_nodes(s+1, :);
    steps = 15; % 细化连线
    x_seg = round(linspace(p1(1), p2(1), steps));
    y_seg = round(linspace(p1(2), p2(2), steps));
    z_seg = round(linspace(p1(3), p2(3), steps));
    grid_path_star = [grid_path_star, [x_seg; y_seg; z_seg]];
end

path_indices_star = sub2ind(grid_sz, grid_path_star(1,:), grid_path_star(2,:), grid_path_star(3,:));
path_indices_star = unique(path_indices_star, 'stable');
star_path = path_indices_star';
star_path_raw = p_ind(path_indices_star)';

save('path_ind.mat', 'star_path', 'star_path_raw', 'star_path_ideal', '-append');
fprintf('Star Path 生成成功！\n');

%% --- 新增路径 C：随高度收缩的六边形棱台轨迹 (Hexagon Frustum Path) ---
% 该路径由若干层六边形组成，每一层的高度增加而半径减小，层与层之间连通
z_steps = [1, 3, 5, 7, 9]; % 选取的 Z 层索引
r_steps = [4.5, 3.5, 2.5, 1.5, 0.5]; % 每层对应的半径（网格单位）
center = [5.5, 5.5]; % 网格中心

hex_grid_path = [];
angles = linspace(0, 2*pi, 7); % 六边形 6 个顶点 + 闭合点

for i = 1:length(z_steps)
    z = z_steps(i);
    r = r_steps(i);
    
    % 当前层的六边形顶点坐标
    x_hex = center(1) + r * cos(angles);
    y_hex = center(2) + r * sin(angles);
    
    % 将顶点连接成离散步进路径
    for v = 1:length(angles)-1
        p1 = [x_hex(v), y_hex(v), z];
        p2 = [x_hex(v+1), y_hex(v+1), z];
        
        pts_count = 10; % 边上的细分点数
        seg_x = round(linspace(p1(1), p2(1), pts_count));
        seg_y = round(linspace(p1(2), p2(2), pts_count));
        seg_z = round(linspace(p1(3), p2(3), pts_count));
        
        hex_grid_path = [hex_grid_path, [seg_x; seg_y; seg_z]];
    end
    
    % 如果不是最后一层，添加连接到下一层起点的路径
    if i < length(z_steps)
        r_next = r_steps(i+1);
        z_next = z_steps(i+1);
        p_current_end = [x_hex(end), y_hex(end), z];
        p_next_start = [center(1) + r_next * cos(angles(1)), center(2) + r_next * sin(angles(1)), z_next];
        
        pts_count_up = 8;
        up_x = round(linspace(p_current_end(1), p_next_start(1), pts_count_up));
        up_y = round(linspace(p_current_end(2), p_next_start(2), pts_count_up));
        up_z = round(linspace(p_current_end(3), p_next_start(3), pts_count_up));
        hex_grid_path = [hex_grid_path, [up_x; up_y; up_z]];
    end
end

% 转换理想路径坐标（物理单位）
hex_path_ideal = [];
for i = 1:length(z_steps)
    z = z_steps(i);
    r = r_steps(i);
    x_hex_phys = -0.1 + (center(1) + r * cos(angles) - 1) * (0.2/9);
    y_hex_phys = -0.1 + (center(2) + r * sin(angles) - 1) * (0.2/9);
    z_hex_phys = -0.1 + (z - 1) * (0.2/9);
    hex_path_ideal = [hex_path_ideal, [x_hex_phys; y_hex_phys; repmat(z_hex_phys, 1, length(angles))]];
end

% 限制在 [1, 10] 范围内
hex_grid_path = max(1, min(10, hex_grid_path));

path_indices_hex = sub2ind(grid_sz, hex_grid_path(1,:), hex_grid_path(2,:), hex_grid_path(3,:));
path_indices_hex = unique(path_indices_hex, 'stable');
hex_path = path_indices_hex';
hex_path_raw = p_ind(path_indices_hex)';

save('path_ind.mat', 'hex_path', 'hex_path_raw', 'hex_path_ideal', '-append');
fprintf('Hexagon Frustum Path 生成成功！\n');

%% --- 新增路径 D：MAGNETIC 字母艺术轨迹 (MAGNETIC Art Path) ---
% 该路径将 "MAGNETIC" 八个字母分布在 x=1..8 的 8 个切层中
letters_str = 'MAGNETIC'; 
magnetic_indices = [];
magnetic_char_id_map = []; % 用于标记每个点所属的字母索引
magnetic_path_ideal = [];

% 检查 p_all 的真实边界
fprintf('p_all 坐标范围: X[%.3f, %.3f], Y[%.3f, %.3f], Z[%.3f, %.3f]\n', ...
    min(p_all(1,:)), max(p_all(1,:)), min(p_all(2,:)), max(p_all(2,:)), min(p_all(3,:)), max(p_all(3,:)));
fprintf('p_all(:,1) = [%.4f, %.4f, %.4f]\n', p_all(1,1), p_all(2,1), p_all(3,1));

% 定义网格到物理空间的映射函数 (保持与 p_all 一致)
grid_to_phys = @(g) -0.1 + (g - 1) * (0.2/9);

for i = 1:numel(letters_str)
    ch = letters_str(i);
    ix = i; 
    strokes = getLetterStrokes(ch);
    
    for s = 1:numel(strokes)
        P = strokes{s}; % Nx2 [y z]
        stroke_x = []; stroke_y = []; stroke_z = [];
        
        for k = 1:size(P, 1)-1
            p1 = P(k, :); p2 = P(k+1, :);
            dist = max(abs(p2 - p1));
            num_pts = max(5, round(dist * 10)); 
            
            y_pts_grid = linspace(p1(1), p2(1), num_pts);
            z_pts_grid = linspace(p1(2), p2(2), num_pts);
            
            % 1. 采样用于实验的离散网格索引
            iy_grid = round(y_pts_grid); 
            iz_grid = round(z_pts_grid);
            iy_grid = max(1, min(10, iy_grid)); 
            iz_grid = max(1, min(10, iz_grid));
            
            % 重要：ix, iy, iz 的对应关系必须与 sub2ind 时一致
            % 这里的 lin_idxs 必须与 grid_sz = [10,10,10] 匹配
            lin_idxs = ix + (iy_grid - 1) * 10 + (iz_grid - 1) * 100;
            magnetic_indices = [magnetic_indices, lin_idxs];
            magnetic_char_id_map = [magnetic_char_id_map, repmat(i, 1, numel(lin_idxs))];
            
            % 2. 构造连续的理想物理坐标
            stroke_y = [stroke_y, grid_to_phys(y_pts_grid)];
            stroke_z = [stroke_z, grid_to_phys(z_pts_grid)];
            stroke_x = [stroke_x, repmat(grid_to_phys(ix), 1, numel(y_pts_grid))];
        end
        magnetic_path_ideal = [magnetic_path_ideal, [stroke_x; stroke_y; stroke_z], [NaN; NaN; NaN]];
    end
end

% 去除 magnetic_indices 中相邻的重复点，保持路径精简
diff_idx = [true, diff(magnetic_indices) ~= 0];
magnetic_indices = magnetic_indices(diff_idx);
magnetic_char_id_map = magnetic_char_id_map(diff_idx);

% 计算每个字符占据的点的长度向量
magnetic_char_lengths = zeros(1, numel(letters_str));
for i = 1:numel(letters_str)
    magnetic_char_lengths(i) = sum(magnetic_char_id_map == i);
end

magnetic_path_raw = p_ind(magnetic_indices)';

save('path_ind.mat', 'magnetic_path_raw', 'magnetic_path_ideal', 'magnetic_char_lengths', '-append');
fprintf('MAGNETIC Art Path 生成成功！采样点数: %s\n', mat2str(magnetic_char_lengths));

%% --- 增强版可视化 (仅展示 MAGNETIC 路径) ---
figure('Color', 'w', 'Name', 'MAGNETIC 3D Art Path', 'Position', [100 100 800 600]);

% 绘制参考全网格（淡灰色点阵）
scatter3(p_all(1,:), p_all(2,:), p_all(3,:), 5, [0.9 0.9 0.9], 'filled', 'MarkerEdgeAlpha', 0.1);
hold on;

% 绘制理想路径 (红色虚线，用于展示设计的形状)
plot3(magnetic_path_ideal(1,:), magnetic_path_ideal(2,:), magnetic_path_ideal(3,:), ...
      'r--', 'LineWidth', 1.5, 'DisplayName', 'Ideal Path');

% 绘制实际选中的网格点 (蓝色圆点)
pts = p_all(:, magnetic_indices);
plot3(pts(1,:), pts(2,:), pts(3,:), 'b.', 'MarkerSize', 12, 'DisplayName', 'Sampled Grid Points'); 

% 标记起点
plot3(pts(1,1), pts(2,1), pts(3,1), 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
% 标记终点
plot3(pts(1,end), pts(2,end), pts(3,end), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

grid on; axis equal; box on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('MAGNETIC 3D Art Path (Fixed Discontinuity & Offset)');
legend('show', 'Location', 'northeastoutside');
view(135, 30);

function plot_path(p_all, idx, ideal_pts, title_str)
    % 此函数目前已弃用，保留以维持脚本结构完整性
end

function strokes = getLetterStrokes(ch)
% Return a cell array of polylines. Each polyline is Nx2 [y z], integer grid.
% Using an "Outline/Ribbon" style based on the 10x10 grid.
% 改用整数或 .0 坐标以确保点线重合
switch upper(ch)
    case 'M'
        % 重新设计 M，确保其坐标尽量在网格线上
        strokes = {[2 2; 2 10; 4 10; 6 8; 8 10; 10 10; 10 2; 8 2; 8 7; 6 5; 4 7; 4 2; 2 2]};
    case 'A'
        strokes = {[2 2; 2 7; 5 10; 7 10; 10 7; 10 2; 8 2; 8 4; 5 4; 5 6; 8 6; 8 7; 7 8; 5 8; 4 7; 4 2; 2 2]};
    case 'G'
        strokes = {[1 1; 1 9; 8 9; 8 7; 3 7; 3 3; 7 3; 7 4; 5 4; 5 6; 9 6; 9 1; 1 1]+1};
    case 'N'
        strokes = {[2 2; 2 10; 4 10; 8 6; 8 10; 10 10; 10 2; 8 2; 4 6; 4 2; 2 2]};
    case 'E'
        strokes = {[2 2; 2 10; 10 10; 10 8; 4 8; 4 7; 9 7; 9 5; 4 5; 4 4; 10 4; 10 2; 2 2]};
    case 'T'
        strokes = {[2 8; 2 10; 10 10; 10 8; 7 8; 7 2; 5 2; 5 8; 2 8]};
    case 'I'
        strokes = {[2 2; 2 4; 5 4; 5 8; 2 8; 2 10; 10 10; 10 8; 7 8; 7 4; 10 4; 10 2; 2 2]};
    case 'C'
        strokes = {[2 2; 2 10; 10 10; 10 8; 4 8; 4 4; 10 4; 10 2; 2 2]};
    otherwise
        strokes = {};
end
end

