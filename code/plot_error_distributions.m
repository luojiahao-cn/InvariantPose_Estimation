function plot_error_distributions(results)
% 改进的误差分布可视化函数

fields = fieldnames(results);
algorithms = {};
pos_errors = {};
rot_errors = {};
field_errors = {};

for i = 1:numel(fields)
    fname = fields{i};
    if endsWith(fname, '_pos_error')
        alg = extractBefore(fname, '_pos_error');
        algorithms{end+1} = alg;
        pos_errors{end+1} = [results.(fname)];
        rot_field = [alg '_rot_error'];
        if isfield(results, rot_field)
            rot_errors{end+1} = [results.(rot_field)];
        else
            rot_errors{end+1} = [];
        end
        field_field = [alg '_field_error'];
        if isfield(results, field_field)
            field_errors{end+1} = [results.(field_field)];
        else
            field_errors{end+1} = [];
        end
    end
end

figure('Color', 'white', 'Position', [100, 100, 1600, 800], 'Name', '误差分布');

% 位置误差分布 - KDE
subplot(2,3,1); hold on;
colors = lines(numel(algorithms));
for i = 1:numel(algorithms)
    [f, x] = ksdensity(pos_errors{i});
    plot(x, f, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
title('位置误差分布 (KDE)', 'FontSize', 14);
xlabel('位置误差 (m)', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend(strrep(algorithms, 'prop', '所提算法'), 'FontSize', 10, 'Location', 'best');
grid on; box on;

% 旋转误差分布 - KDE
subplot(2,3,2); hold on;
for i = 1:numel(algorithms)
    [f, x] = ksdensity(rot_errors{i});
    plot(x, f, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
title('旋转误差分布 (KDE)', 'FontSize', 14);
xlabel('旋转误差', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend(strrep(algorithms, 'prop', '所提算法'), 'FontSize', 10, 'Location', 'best');
grid on; box on;

% 磁场误差分布 - KDE（只显示有数据的算法，自动排除init）
subplot(2,3,3); hold on;
valid_idx = find(~cellfun(@isempty, field_errors) & cellfun(@(x) numel(x)>0, field_errors));
field_errors_valid = field_errors(valid_idx);
algorithms_valid = algorithms(valid_idx);
for i = 1:numel(field_errors_valid)
    [f, x] = ksdensity(field_errors_valid{i});
    plot(x, f, '-', 'LineWidth', 2, 'Color', colors(valid_idx(i),:));
end
title('磁场误差分布 (KDE)', 'FontSize', 14);
xlabel('磁场误差', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend(strrep(algorithms_valid, 'prop', '所提算法'), 'FontSize', 10, 'Location', 'best');
grid on; box on;

% 位置误差箱线图
subplot(2,3,4);
boxplot(cell2mat(cellfun(@(x) x(:), pos_errors, 'UniformOutput', false)), ...
    'Labels', strrep(algorithms, 'prop', 'Proposed'));
title('位置误差分布比较', 'FontSize', 14);
ylabel('位置误差 (m)', 'FontSize', 12);
grid on;

% 旋转误差箱线图
subplot(2,3,5);
boxplot(cell2mat(cellfun(@(x) x(:), rot_errors, 'UniformOutput', false)), ...
    'Labels', strrep(algorithms, 'prop', 'Proposed'));
title('旋转误差分布比较', 'FontSize', 14);
ylabel('旋转误差', 'FontSize', 12);
grid on;

% 磁场误差箱线图（只显示有数据的算法，自动排除init）
subplot(2,3,6);
boxplot(cell2mat(cellfun(@(x) x(:), field_errors_valid, 'UniformOutput', false)), ...
    'Labels', strrep(algorithms_valid, 'prop', 'Proposed'));
title('磁场误差分布比较', 'FontSize', 14);
ylabel('磁场误差', 'FontSize', 12);
grid on;

sgtitle('算法误差分布比较', 'FontSize', 16, 'FontWeight', 'bold');
end