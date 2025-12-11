function plot_batch_results(batch_results)
% PLOT_BATCH_RESULTS 绘制批量实验结果
% 输入：
%   batch_results - 由 run_batch_experiments 返回的批量结果结构体数组

fontsize = 16;

num_points = length(batch_results);
num_methods = 4;  % LM, ELM, Ours, Rlm

% 提取所有测试点的坐标
all_coords = [batch_results.test_point];
all_p_true = [all_coords.p_true];  % 3×num_points 矩阵
% all_p_init = [all_coords.p_init];
% 提取所有batch_results(N).results.p_init
all_p_init = [];
for i = 1:num_points
    all_p_init = [all_p_init, batch_results(i).results.p_init];
end

% 提取所有测试点的误差（使用summary中的均值）
pos_error_matrix = zeros(num_points, num_methods);
rot_error_matrix = zeros(num_points, num_methods);

method_names = {'ours', 'fischer', 'elm', 'lm'};
method_labels = {'Ours', 'Fischer', 'ELM', 'LM'};

for i = 1:num_points
    for j = 1:num_methods
        method = method_names{j};
        % 获取位置误差和旋转误差的均值
        pos_error_matrix(i, j) = batch_results(i).summary.(method).pos_mean;
        rot_error_matrix(i, j) = 2*asin(batch_results(i).summary.(method).rot_mean/(2*sqrt(2)));
        time_matrix(i, j) = batch_results(i).summary.(method).time_mean;
    end
end

% 先将位置误差和旋转误差直接相加，再全局归一化
rot_error_matrix_ = rot_error_matrix;
rot_error_matrix_(rot_error_matrix_ > 0.1) = 0.1;
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
figure('Color', 'white', 'units', 'centimeters', 'Position', [0, 0, 30, 16]);

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
    % scatter(all_p_true(1, :), all_p_true(2, :), 50, colors, 'filled');
    scatter(all_p_init(1, :), all_p_init(2, :), 50, colors, 'filled');
    
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
    title(sprintf(method_labels{j}), 'FontSize', fontsize, 'FontWeight', 'bold');
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
    % scatter(all_p_true(1, :), all_p_true(3, :), 50, colors, 'filled');
    scatter(all_p_init(1, :), all_p_init(3, :), 50, colors, 'filled');
    
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
    nexttile(subtl, 1); % 13,14,15,16
    % hold on;

    % 取出该方法的所有测试点位置和旋转误差
    rot_errors = rot_error_matrix(:, j);

    % 构造数据用于boxchart（横向）
    % y: [位置误差; 旋转误差]
    % group: 1=位置, 2=旋转
    y = rot_errors;
    
    boxplot(y, 'Labels', {'Rotation'}, 'Orientation', 'horizontal', 'Widths', 2);

    % 设置颜色：灰色，黑色边框
    % boxchart 会为每个类别创建一个 BoxChart 对象
    colors = slanCM('gor', 2);

    boxes = findobj(gca, 'Tag', 'Box');
    set(boxes, 'Color', colors(1, :));
    medians = findobj(gca, 'Tag', 'Median');
    set(medians, 'Color', colors(2, :));  % 中位线颜色设为黑色

    if j == 1
        yticklabels({'$e_{\mbox{\boldmath ${R}$}}$ [rad]'});
    else
        yticklabels('');
    end

    % 设置坐标轴（使用 categorical 值）
    set(gca, 'FontSize', fontsize, 'TickDir', 'out', 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    % xlabel('RMSE [-]', 'FontSize', fontsize, 'Interpreter', 'latex');
    xlim([0, 0.15]);
    % axis auto;
    box on;

    nexttile(subtl, 2); % 15,16
    pos_errors = pos_error_matrix(:, j);
    y = pos_errors;
    
    boxplot(1e3 * y, 'Labels', {'Position'}, 'Orientation', 'horizontal');

    % 设置颜色：灰色，黑色边框
    % boxchart 会为每个类别创建一个 BoxChart 对象
    colors = slanCM('gor', 2);

    boxes = findobj(gca, 'Tag', 'Box');
    set(boxes, 'Color', colors(1, :));
    medians = findobj(gca, 'Tag', 'Median');
    set(medians, 'Color', colors(2, :), 'LineWidth', 1.5);  % 中位线颜色设为黑色

    if j == 1
        yticklabels({'$e_{\mbox{\boldmath ${p}$}}$ [mm]'});
    else
        yticklabels('');
    end

    % 设置坐标轴（使用 categorical 值）
    set(gca, 'FontSize', fontsize, 'TickDir', 'out', 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    % xlabel('RMSE', 'FontSize', fontsize, 'Interpreter', 'latex');
    % xlim([0, 40]);
    axis auto;
    box on;

end

%%
mean(time_matrix)
% save('../results/initial_cond_config1.mat');
% export_fig('../figures/initial_cond_config1.png', '-r600', '-nocrop');
end

