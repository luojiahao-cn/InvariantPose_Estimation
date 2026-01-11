%% MAGNETIC in Y-Z plane (vector-stroke, angular style)
% One figure, 8 letters tiled along y direction.
clear; clc; close all;

% ---- Layout ----
cellW = 10;           % letter box width in y
cellH = 10;           % letter box height in z
gap   = 2;            % spacing between letters (in y)
letters = 'MAGNETIC';
nL = numel(letters);

% Total canvas range in y
Ymax = nL*cellW + (nL-1)*gap;
Zmax = cellH;

figure('Color','w'); hold on; axis equal;
xlim([0, Ymax]); ylim([0, Zmax+1]);
xlabel('y'); ylabel('z');
title('MAGNETIC on the y-z plane (angular vector strokes)');

% ---- Draw background grid (optional but helpful) ----
drawGrid = true;
if drawGrid
    % vertical grid lines
    for yy = 0:Ymax
        plot([yy yy],[0 Zmax], 'Color',[0.85 0.85 0.85], 'LineWidth',0.5);
    end
    % horizontal grid lines
    for zz = 0:Zmax
        plot([0 Ymax],[zz zz], 'Color',[0.85 0.85 0.85], 'LineWidth',0.5);
    end
end

% ---- Stroke style ----
strokeColor = [0.05 0.35 0.45];
LW = 2;          % For outline style, thinner lines are better
LW2 = 2;         
capStyle = 'square';

% ---- Draw each letter ----
yOffset = 0;
for idx = 1:nL
    ch = letters(idx);
    baseY = yOffset;   % left edge of this letter box in global y

    % draw bounding box (optional)
    plot([baseY baseY+cellW baseY+cellW baseY baseY], ...
         [0 0 cellH cellH 0], 'Color',[0.75 0.75 0.75], 'LineWidth',1);

    % get strokes in local coords (y,z) within 1..10
    strokes = getLetterStrokes(ch);

    % plot strokes (translate by baseY)
    meta = strokesMeta(ch);
    for s = 1:numel(strokes)
        P = strokes{s};              % Nx2 array [y z]
        P(:,1) = P(:,1) + baseY;     % translate y
        
        % Check if this stroke is marked as secondary
        isSecondary = false;
        if ~isempty(meta.secondary) && numel(meta.secondary) >= s
            isSecondary = meta.secondary(s);
        end
        
        if isSecondary
            plot(P(:,1), P(:,2), '-', 'Color', strokeColor, ...
                'LineWidth', LW2);
        else
            plot(P(:,1), P(:,2), '-', 'Color', strokeColor, ...
                'LineWidth', LW);
        end
    end

    % letter label
    text(baseY+cellW/2, cellH+0.6, ch, 'HorizontalAlignment','center', ...
        'FontName','Arial','FontWeight','bold','FontSize',12);

    yOffset = yOffset + cellW + gap;
end

set(gca,'YDir','normal');  % z upward
box on; grid off;

%% -------- 路径解析与索引转换 (MAGNETIC 3D Path) --------
% 我们将 8 个字母分布在 x = 1, 2, ..., 8 的 8 个切层中
% 每个字母在各自的 x 层内占据一个 10x10 的 y-z 区域

fprintf('\n正在计算 MAGNETIC 路径索引...\n');

% 加载原始索引映射 (用于转换逻辑索引到原始扫描顺序)
if exist('path_ind.mat', 'file')
    S = load('path_ind.mat');
    if isfield(S, 'p_ind')
        p_ind = S.p_ind;
    else
        p_ind = 1:1000;
    end
else
    p_ind = 1:1000;
end

magnetic_indices = [];
magnetic_path_ideal = []; % 用于保存物理坐标 (x, y, z)

for i = 1:nL
    ch = letters(i);
    x_val = i; % 字母 A-H 分别对应 x = 1-8
    
    strokes = getLetterStrokes(ch);
    
    for s = 1:numel(strokes)
        P = strokes{s}; % Nx2 [y z]
        
        % 生成理想路径坐标 (保持浮点，用于绘图)
        ideal_seg = [repmat(x_val, size(P, 1), 1), P];
        magnetic_path_ideal = [magnetic_path_ideal; ideal_seg];
        
        % 离散化笔画为网格索引
        for k = 1:size(P, 1)-1
            p1 = P(k, :);
            p2 = P(k+1, :);
            
            % 根据线段长度确定插值点数，确保路径连续
            dist = max(abs(p2 - p1));
            num_pts = max(2, round(dist * 5)); % 适当增加采样密度
            
            y_pts = linspace(p1(1), p2(1), num_pts);
            z_pts = linspace(p1(2), p2(2), num_pts);
            
            % 转换为整数网格坐标 (1-10)
            ix = x_val;
            iy = round(y_pts);
            iz = round(z_pts);
            
            % 限制在 1-10 范围内
            iy = max(1, min(10, iy));
            iz = max(1, min(10, iz));
            
            % 计算线性索引: (x-1) + (y-1)*10 + (z-1)*100 + 1
            lin_idx = (ix - 1) + (iy - 1) * 10 + (iz - 1) * 100 + 1;
            magnetic_indices = [magnetic_indices, lin_idx];
        end
    end
end

% 去重并保持顺序 (保持路径的连续性)
% 先去除相邻的重复点
diff_idx = [true, diff(magnetic_indices) ~= 0];
magnetic_indices = magnetic_indices(diff_idx);

% 映射回 raw 原始索引
magnetic_path_raw = p_ind(magnetic_indices)';

% 保存到文件
save('path_ind.mat', 'magnetic_path_raw', 'p_ind', '-append');
fprintf('Success! MAGNETIC 路径已保存到 path_ind.mat\n');
fprintf(' - magnetic_path_raw: 共 %d 个点\n', length(magnetic_path_raw));

%% -------- Local functions --------
function strokes = getLetterStrokes(ch)
% Return a cell array of polylines. Each polyline is Nx2 [y z], integer grid.
% Using an "Outline/Ribbon" style based on the 10x10 grid.

Y0=1.5; Y1=8.5; Z0=1.5; Z1=8.5;
YM=5; ZM=5;
W = 1.2; % Stroke width for ribbon

switch upper(ch)
    case 'M'
        % Outline path for M
        strokes = {
            [Y0 Z0; Y0 Z1; Y0+W+0.3 Z1; YM ZM; Y1-W-0.3 Z1; Y1 Z1; Y1 Z0; ...
             Y1-W Z0; Y1-W Z1-1.3; YM ZM-1.8; Y0+W Z1-1.3; Y0+W Z0; Y0 Z0]
        };

    case 'A'
        % Refined A based on the latest stencil image (Trapezoid arch with bottom notch and inner triangle)
        strokes = {
            [4 8.5; 6 8.5; 9 1.5; 7 1.5; 6 3.5; 4 3.5; 3 1.5; 1 1.5; 4 8.5], ... % Outer Frame
            [4.1 5.2; 5.9 5.2; 5 7.7; 4.1 5.2]                                 % Inner Triangle
        };

    case 'G'
        % Blocky G outline
        strokes = {
            [Y1 Z1; Y0 Z1; Y0 Z0; Y1 Z0; Y1 ZM; YM ZM; YM ZM-W; Y1-W ZM-W; Y1-W Z0+W; Y0+W Z0+W; Y0+W Z1-W; Y1 Z1-W; Y1 Z1]
        };

    case 'N'
        % N outline
        strokes = {
            [Y0 Z0; Y0 Z1; Y0+W Z1; Y1-W Z0+W+1; Y1-W Z1; Y1 Z1; Y1 Z0; Y1-W Z0; Y0+W Z1-W-1; Y0+W Z0; Y0 Z0]
        };

    case 'E'
        % E outline
        strokes = {
            [Y0 Z0; Y0 Z1; Y1 Z1; Y1 Z1-W; Y0+W Z1-W; Y0+W ZM+W/2; Y1-1 ZM+W/2; Y1-1 ZM-W/2; Y0+W ZM-W/2; Y0+W Z0+W; Y1 Z0+W; Y1 Z0; Y0 Z0]
        };

    case 'T'
        % T outline
        strokes = {
            [Y0 Z1; Y1 Z1; Y1 Z1-W; YM+W/2 Z1-W; YM+W/2 Z0; YM-W/2 Z0; YM-W/2 Z1-W; Y0 Z1-W; Y0 Z1]
        };

    case 'I'
        % I outline
        strokes = {
            [Y0 Z1; Y1 Z1; Y1 Z1-W; YM+W/2 Z1-W; YM+W/2 Z0+W; Y1 Z0+W; Y1 Z0; Y0 Z0; Y0 Z0+W; YM-W/2 Z0+W; YM-W/2 Z1-W; Y0 Z1-W; Y0 Z1]
        };

    case 'C'
        % C outline
        strokes = {
            [Y1 Z1; Y0 Z1; Y0 Z0; Y1 Z0; Y1 Z0+W; Y0+W Z0+W; Y0+W Z1-W; Y1 Z1-W; Y1 Z1]
        };

    otherwise
        strokes = {};
end
end

function meta = strokesMeta(ch)
% Mark which strokes are secondary (thinner) for each letter.
% Helps make "thickening" strokes slightly thinner than main strokes.
meta.secondary = [];
switch upper(ch)
    case 'M'
        meta.secondary = [false false false true];
    case 'N'
        meta.secondary = [false false false true];
    otherwise
        % all main
        % (we will handle missing fields in caller)
        meta.secondary = [];
end
end
