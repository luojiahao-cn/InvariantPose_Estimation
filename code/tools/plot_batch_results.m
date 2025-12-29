function plot_batch_results(batch_results)
% PLOT_BATCH_RESULTS 绘制批量实验结果
% 输入：
%   batch_results - 由 run_batch_experiments 返回的批量结果结构体数组

fontsize = 14;

num_points = length(batch_results);
num_methods = 4;  % LM, ELM, Ours, Rlm

% 提取所有测试点的坐标
all_coords = [batch_results.test_point];
all_p_true = [all_coords.p_true];  % 3×num_points 矩阵
% all_p_init = [all_coords.p_init];
% 提取所有batch_results(N).results.p_init

% 提取所有测试点的误差（使用summary中的均值）
pos_error_matrix = zeros(num_points, num_methods);
rot_error_matrix = zeros(num_points, num_methods);

method_names = {'ours', 'lm', 'elm', 'fischer'};
method_labels = {'Ours', '(1)', '(4)', 'Fischer'};

for i = 1:num_points
    for j = 1:num_methods
        method = method_names{j};
        % 获取位置误差和旋转误差的均值
        pos_error_matrix(i, j) = batch_results(i).summary.(method).pos_mean;
        r_error_matrix(i, j) = batch_results(i).summary.(method).r_mean;
        rot_error_matrix(i, j) = 2*asin(batch_results(i).summary.(method).rot_mean/(2*sqrt(2)));
        time_matrix(i, j) = batch_results(i).summary.(method).time_mean;
    end
end

% 先将位置误差和旋转误差直接相加，再全局归一化
rot_error_matrix_ = rot_error_matrix;
rot_error_matrix_(rot_error_matrix_ > 0.5) = 0.5;
error_matrix = pos_error_matrix + rot_error_matrix_;

% 归一化综合误差到[0,1]用于颜色映射（全局归一化）
err_max_all = max(error_matrix(:));
err_min_all = min(error_matrix(:));
if err_max_all > err_min_all
    error_norm = (error_matrix - err_min_all) / (err_max_all - err_min_all);
else
    error_norm = ones(size(error_matrix)) * 0.5;
end

% 创建颜色映射
n_colors = 256;
custom_colormap = slanCM('gor', n_colors);

% 创建图形
figure('Color', 'white', 'units', 'centimeters', 'Position', [0, 0, 30, 15]);

tl = tiledlayout(4, num_methods, 'TileSpacing', 'compact', 'Padding', 'compact');
% x-y 平面 scatter，占据前两行, 每个方法合并两行（2行1列）
for j = 1:num_methods
    % 使用 nexttile 的 span 参数合并两行
    % nexttile(tile_number, [row_span, col_span])
    % 第j列，占据第1-2行（合并两行）
    nexttile(j, [2, 1]);
    hold on;
    
    % 根据误差大小映射颜色
    colors = zeros(num_points, 3);
    for i = 1:num_points
        color_idx = max(1, min(n_colors, round(error_norm(i, j) * (n_colors - 1)) + 1));
        colors(i, :) = custom_colormap(color_idx, :);
    end
    
    % 绘制散点图
    scatter(all_p_true(1, :), all_p_true(2, :), 50, colors, 'filled');
    % scatter(all_p_init(1, :), all_p_init(2, :), 50, colors, 'filled');
    
    xlabel('$x$ [mm]', 'FontSize', fontsize, 'Interpreter', 'latex');
    if j == 1
        ylabel('$y$ [mm]', 'FontSize', fontsize, 'Interpreter', 'latex');
        yticks([-0.1, 0, 0.1]);
        yticklabels({'-10', '0', '10'});
    else
        yticklabels('');
    end
    xticks([-0.1, 0, 0.1]);
    xticklabels({'-10', '0', '10'});
    xlim([-0.15, 0.15]);
    ylim([-0.15, 0.15]);
    title(sprintf(method_labels{j}), 'FontSize', fontsize, 'Interpreter', 'latex');
    grid off;
    % axis equal;
    box on;

    set(gca, 'FontSize', fontsize, 'TickLabelInterpreter', 'latex', 'TickDir', 'out', 'LineWidth', 1);
end

% 第三行：x-z 平面 scatter（第9~12列, 只用一行，每方法一列）
for j = 1:num_methods
    nexttile(j + 8); % 第3行：9,10,11,12
    hold on;
    
    % 根据误差大小映射颜色
    colors = zeros(num_points, 3);
    for i = 1:num_points
        color_idx = max(1, min(n_colors, round(error_norm(i, j) * (n_colors - 1)) + 1));
        colors(i, :) = custom_colormap(color_idx, :);
    end
    
    % 绘制散点图
    scatter(all_p_true(1, :), all_p_true(3, :), 50, colors, 'filled');
    
    xlabel('$x$ [mm]', 'FontSize', fontsize, 'Interpreter', 'latex');
    if j == 1
        ylabel('$z$ [mm]', 'FontSize', fontsize, 'Interpreter', 'latex');
        yticks([0, 0.1]);
        yticklabels({0, '10'});
    else
        yticklabels('');
    end
    xlim([-0.15, 0.15]);
    xticks([-0.1, 0, 0.1]);
    xticklabels({'-10', '0', '10'});
    ylim([0, 0.15]);
    % title(sprintf('%s - X-Z平面', method_labels{j}), 'FontSize', 12, 'FontWeight', 'bold');
    grid off;
    % axis equal;
    box on;
    set(gca, 'FontSize', fontsize, 'TickLabelInterpreter', 'latex', 'TickDir', 'out', 'LineWidth', 1);
end

% 第四行：每个方法的误差统计（第13~16列, 每方法一列）
for j = 1:num_methods

    subtl = tiledlayout(tl, 2, 1);
    subtl.Layout.Tile = 12 + j;

    % 位置+旋转误差统计（第13-16列，第j列为第12+j个tile）
    nexttile(subtl, 1);

    rot_errors = rot_error_matrix(:, j);
    
    % 使用 boxchart（横向）
    gc = categorical(repmat("Rotation", num_points, 1));
    bc = boxchart(gc, rot_errors, 'BoxFaceColor', colors(1, :), 'Orientation', 'horizontal');
    % bc.MarkerStyle = 'none';  % 不显示离群点
    
    % % 颜色设置
    % bc.BoxFaceColor = colors(1, :);
    % bc.BoxEdgeColor = 'k';
    % bc.WhiskerLineColor = 'k';
    % bc.MedianLineColor = colors(2, :);
    
    if j == 1
        yticklabels('$e_{\mbox{\boldmath ${R}$}}$ [rad]');
    else
        yticklabels('');
    end
    
    set(gca, 'FontSize', fontsize, 'TickDir', 'out', ...
             'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    box on;
    

    nexttile(subtl, 2);

    pos_errors = 1e3 * pos_error_matrix(:, j);  % 换成 mm
    
    gc = categorical(repmat("Position", num_points, 1));
    bc = boxchart(gc, pos_errors, 'BoxFaceColor', colors(1, :), 'Orientation', 'horizontal');
    % bc.MarkerStyle = 'none';
    
    % 颜色设置
    % bc.BoxFaceColor = colors(1, :);
    % bc.BoxEdgeColor = 'k';
    % bc.WhiskerLineColor = 'k';
    % bc.MedianLineColor = colors(2, :);
    % bc.MedianLineWidth = 1.5;
    
    if j == 1
        yticklabels('$e_{\mbox{\boldmath ${p}$}}$ [mm]');
    else
        yticklabels('');
    end
    
    set(gca, 'FontSize', fontsize, 'TickDir', 'out', ...
             'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    box on;

end

%%
mean(time_matrix)
% save('../results/initial_cond_config1.mat');
% export_fig('../figures/initial_cond_config1.png', '-r600', '-nocrop');
end

