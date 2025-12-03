function plot_error_distributions(results)
% 多算法自适应误差分布可视化

fields = fieldnames(results);
algorithms = {};
pos_errors = {};
rot_errors = {};
field_errors = {};
time_errors = {};

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

% 自动收集耗时数据（支持t_开头字段）
time_fields = fields(startsWith(fields, 't_'));
time_algorithms = cellfun(@(x) extractAfter(x, 't_'), time_fields, 'UniformOutput', false);
for i = 1:numel(algorithms)
    idx = find(strcmp(time_algorithms, algorithms{i}), 1);
    if ~isempty(idx)
        time_errors{end+1} = [results.(time_fields{idx})];
    else
        time_errors{end+1} = [];
    end
end

% 算法显示名映射表
alg_map = containers.Map({'init','lm','elm','Rlm','ours'}, ...
    {'初始','LM','ELM','Rlm','Ours'});
alg_labels = cell(size(algorithms));
for i = 1:numel(algorithms)
    x = algorithms{i};
    if isKey(alg_map, x)
        alg_labels{i} = alg_map(x);
    else
        alg_labels{i} = x;
    end
end

figure('Color', 'white', 'Position', [100, 100, 1800, 900], 'Name', '误差分布');

% 位置误差分布 - KDE
subplot(2,4,1); hold on;
colors = lines(numel(algorithms));
for i = 1:numel(algorithms)
    [f, x] = ksdensity(pos_errors{i});
    plot(x, f, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
title('位置误差分布 (KDE)', 'FontSize', 14);
xlabel('位置误差 (m)', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend(alg_labels, 'FontSize', 10, 'Location', 'best');
grid on; box on;

% 旋转误差分布 - KDE
subplot(2,4,2); hold on;
for i = 1:numel(algorithms)
    [f, x] = ksdensity(rot_errors{i});
    plot(x, f, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
title('旋转误差分布 (KDE)', 'FontSize', 14);
xlabel('旋转误差', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend(alg_labels, 'FontSize', 10, 'Location', 'best');
grid on; box on;

% 磁场误差分布 - KDE（只显示有数据的算法，自动排除init）
subplot(2,4,3); hold on;
valid_idx = find(~cellfun(@isempty, field_errors) & cellfun(@(x) numel(x)>0, field_errors));
field_errors_valid = field_errors(valid_idx);
alg_labels_valid = alg_labels(valid_idx);
for i = 1:numel(field_errors_valid)
    [f, x] = ksdensity(field_errors_valid{i});
    plot(x, f, '-', 'LineWidth', 2, 'Color', colors(valid_idx(i),:));
end
title('磁场误差分布 (KDE)', 'FontSize', 14);
xlabel('磁场误差', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend(alg_labels_valid, 'FontSize', 10, 'Location', 'best');
grid on; box on;

% 位置误差箱线图
subplot(2,4,4);
boxplot(cell2mat(cellfun(@(x) x(:), pos_errors, 'UniformOutput', false)), ...
    'Labels', alg_labels);
title('位置误差分布比较', 'FontSize', 14);
ylabel('位置误差 (m)', 'FontSize', 12);
grid on;

% 旋转误差箱线图
subplot(2,4,5);
boxplot(cell2mat(cellfun(@(x) x(:), rot_errors, 'UniformOutput', false)), ...
    'Labels', alg_labels);
title('旋转误差分布比较', 'FontSize', 14);
ylabel('旋转误差', 'FontSize', 12);
grid on;

% 磁场误差箱线图（只显示有数据的算法，自动排除init）
subplot(2,4,6);
boxplot(cell2mat(cellfun(@(x) x(:), field_errors_valid, 'UniformOutput', false)), ...
    'Labels', alg_labels_valid);
title('磁场误差分布比较', 'FontSize', 14);
ylabel('磁场误差', 'FontSize', 12);
grid on;

% 算法耗时箱线图（只显示有数据的算法）
subplot(2,4,[7 8]);
valid_time_idx = find(~cellfun(@isempty, time_errors) & cellfun(@(x) numel(x)>0, time_errors));
time_errors_valid = time_errors(valid_time_idx);
alg_labels_time = alg_labels(valid_time_idx);
if ~isempty(time_errors_valid)
    boxplot(cell2mat(cellfun(@(x) x(:), time_errors_valid, 'UniformOutput', false)), ...
        'Labels', alg_labels_time);
    title('算法耗时分布比较', 'FontSize', 14);
    ylabel('耗时 (秒)', 'FontSize', 12);
    grid on;
end

sgtitle('算法误差与耗时分布比较', 'FontSize', 16, 'FontWeight', 'bold');
end