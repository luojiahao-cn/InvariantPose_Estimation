function display_statistical_summary(results, num_experiments)
% DISPLAY_STATISTICAL_SUMMARY 显示实验结果的统计摘要
%
% 输入参数：
%   results - 结构体数组，包含每次实验的误差数据
%   num_experiments - 实验次数
%
% 输出：无（在命令窗口显示统计摘要）

fields = fieldnames(results);
algorithms = {};
pos_errors = {};
rot_errors = {};

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
    end
end

init_idx = find(strcmp(algorithms, 'init'));
fprintf('\n===== 统计摘要 (%d次实验) =====\n', num_experiments);
fprintf('初始位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
    mean(pos_errors{init_idx}), std(pos_errors{init_idx}));
fprintf('初始旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
    mean(rot_errors{init_idx}), std(rot_errors{init_idx}));

for i = 1:numel(algorithms)
    if strcmp(algorithms{i}, 'init'), continue; end
    fprintf('\n%s算法:\n', strrep(algorithms{i}, 'prop', '所提'));
    fprintf('位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
        mean(pos_errors{i}), std(pos_errors{i}));
    fprintf('旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
        mean(rot_errors{i}), std(rot_errors{i}));
    pos_reduction = (pos_errors{init_idx} - pos_errors{i}) ./ pos_errors{init_idx} * 100;
    rot_reduction = (rot_errors{init_idx} - rot_errors{i}) ./ rot_errors{init_idx} * 100;
    fprintf('位置误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
        mean(pos_reduction), std(pos_reduction));
    fprintf('旋转误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
        mean(rot_reduction), std(rot_reduction));
end
end