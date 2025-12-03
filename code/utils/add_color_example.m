% 示例：如何添加新颜色到 slanCM_Data
% 
% 使用步骤：
% 1. 先运行 tmp.m 加载数据并查看结构
% 2. 然后运行此脚本添加新颜色
% 3. 最后保存文件

%% 确保已加载数据
if ~exist('slanCM_Data', 'var')
    fprintf('正在加载数据...\n');
    load('D:\Softwares\MATLAB\toolbox\slan_CM\slanCM_Data.mat');
end

%% 示例1: 添加一个简单的蓝色渐变
% 从深蓝色到白色的渐变
color_name = 'myblue';
rgb_start = [0.2 0.4 0.8];  % 起始RGB值
n_points = 256;  % 颜色映射点数

% 生成颜色映射（从指定颜色到白色）
color_map = zeros(n_points, 3);
for i = 1:n_points
    t = (i - 1) / (n_points - 1);
    % 线性插值：从rgb_start到白色[1,1,1]
    color_map(i, :) = rgb_start + (1 - rgb_start) * t;
end

% 添加到结构体
slanCM_Data.(color_name) = color_map;
fprintf('已添加颜色: %s\n', color_name);

%% 示例2: 添加多色渐变（可选）
% 如果您想添加多色渐变，可以这样做：
% color_name2 = 'sunset';
% control_points = [0.8 0.2 0.2; 1 0.8 0.2; 1 1 1];  % 从红到黄到白
% color_map2 = zeros(n_points, 3);
% for i = 1:n_points
%     t = (i - 1) / (n_points - 1);
%     if t < 0.5
%         % 前一半：从第一个颜色到第二个颜色
%         color_map2(i, :) = (1 - t*2) * control_points(1, :) + (t*2) * control_points(2, :);
%     else
%         % 后一半：从第二个颜色到第三个颜色
%         color_map2(i, :) = (1 - (t-0.5)*2) * control_points(2, :) + ((t-0.5)*2) * control_points(3, :);
%     end
% end
% slanCM_Data.(color_name2) = color_map2;
% fprintf('已添加颜色: %s\n', color_name2);

%% 显示添加后的颜色列表
fprintf('\n当前所有颜色映射:\n');
field_names = fieldnames(slanCM_Data);
for i = 1:length(field_names)
    fprintf('  %d. %s\n', i, field_names{i});
end

%% 保存文件
fprintf('\n正在保存文件...\n');
save('D:\Softwares\MATLAB\toolbox\slan_CM\slanCM_Data.mat', 'slanCM_Data', '-v7.3');
fprintf('保存完成！\n');
fprintf('现在可以使用: slanCM(''%s'', n) 来调用新颜色\n', color_name);

