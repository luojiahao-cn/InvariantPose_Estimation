function plot_init_results(batch_results)
% PLOT_BATCH_RESULTS 绘制批量实验结果
% 输入：
%   batch_results - 由 run_batch_experiments 返回的批量结果结构体数组
% load('../../results/initial_cond_config1_20251211.mat', 'batch_results');
num_points = length(batch_results);
num_methods = 4;  % LM, ELM, Ours, Rlm

% 提取所有测试点的坐标
% all_coords = [batch_results.test_point];
% all_p_true = [all_coords.p_true];  % 3×num_points 矩阵
% all_p_init = [all_coords.p_init];
% 提取所有batch_results(N).results.p_init
all_p_init = [];
for i = 1:num_points
    all_p_init = [all_p_init, batch_results(i).results.p_init];
end

% 提取所有测试点的误差（使用summary中的均值）
pos_error_matrix = zeros(num_points, num_methods);
rot_error_matrix = zeros(num_points, num_methods);
time_matrix = zeros(num_points, num_methods);

method_names = {'ours', 'lm', 'elm', 'fischer'};

for i = 1:num_points
    for j = 1:num_methods
        method = method_names{j};
        % 获取位置误差和旋转误差的均值
        pos_error_matrix(i, j) = batch_results(i).summary.(method).pos_mean;
        rot_error_matrix(i, j) = 2*asin(batch_results(i).summary.(method).rot_mean/(2*sqrt(2)));
        time_matrix(i, j) = batch_results(i).summary.(method).time_mean;
    end
end

figure;
set(gcf, 'color', 'w', 'units', 'centimeters', ...
    'Position', [0, 0, 19.8, 6.3]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

fontsize = 16;
method_labels = {'Ours', '(1)', '(4)', 'Fischer'};

%% 1. Pos
nexttile;
hold on;
for j = 1:num_methods
    boxchart(j * ones(num_points, 1), 1e3 * pos_error_matrix(:, j)); 
end
hold off;
set(gca, 'XTick', 1:num_methods, 'XTickLabel', method_labels, ...
         'TickLabelInterpreter', 'latex');
ylabel({'$e_{\mbox{\boldmath ${p}$}}$ [mm]'}, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', fontsize, 'TickDir', 'out', 'LineWidth', 1);

%% 2. Rot
nexttile;
hold on;
for j = 1:num_methods
    boxchart(j * ones(num_points, 1), rot_error_matrix(:, j)); 
end
hold off;
set(gca, 'XTick', 1:num_methods, 'XTickLabel', method_labels, ...
         'TickLabelInterpreter', 'latex');
ylabel({'$e_{\mbox{\boldmath ${R}$}}$ [rad]'}, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', fontsize, 'TickDir', 'out', 'LineWidth', 1);

%% 3. Time
nexttile;
hold on;
for j = 1:num_methods
    boxchart(j * ones(num_points, 1), 1e3 * time_matrix(:, j)); 
end
hold off;
set(gca, 'XTick', 1:num_methods, 'XTickLabel', method_labels, ...
         'TickLabelInterpreter', 'latex');
ylabel('$t$ [ms]', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', fontsize, 'TickDir', 'out', 'LineWidth', 1);
%%
% export_fig('../../figures/init_results_config1.png', '-png', '-r600', '-nocrop');
end

