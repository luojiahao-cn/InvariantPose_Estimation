function display_statistical_summary(results, num_experiments)
% DISPLAY_STATISTICAL_SUMMARY 显示实验结果的统计摘要
%
% 输入参数：
%   results - 结构体数组，包含每次实验的误差数据
%   num_experiments - 实验次数
%
% 输出：无（在命令窗口显示统计摘要）

% 提取数据（增加Union算法，保持prop在最后）
init_pos_errors = [results.init_pos_error];
init_rot_errors = [results.init_rot_error];
lm_pos_errors = [results.lm_pos_error];
lm_rot_errors = [results.lm_rot_error];
elm_pos_errors = [results.elm_pos_error];
elm_rot_errors = [results.elm_rot_error];
union_pos_errors = [results.union_pos_error];  % 新增Union位置误差
union_rot_errors = [results.union_rot_error];  % 新增Union旋转误差
prop_pos_errors = [results.prop_pos_error];    % prop保持最后
prop_rot_errors = [results.prop_rot_error];    % prop保持最后

% 计算误差减少比例（保持prop在最后）
pos_error_reduction_lm = (init_pos_errors - lm_pos_errors) ./ init_pos_errors * 100;
pos_error_reduction_elm = (init_pos_errors - elm_pos_errors) ./ init_pos_errors * 100;
pos_error_reduction_union = (init_pos_errors - union_pos_errors) ./ init_pos_errors * 100;  % Union
pos_error_reduction_prop = (init_pos_errors - prop_pos_errors) ./ init_pos_errors * 100;   % prop最后

rot_error_reduction_lm = (init_rot_errors - lm_rot_errors) ./ init_rot_errors * 100;
rot_error_reduction_elm = (init_rot_errors - elm_rot_errors) ./ init_rot_errors * 100;
rot_error_reduction_union = (init_rot_errors - union_rot_errors) ./ init_rot_errors * 100;  % Union
rot_error_reduction_prop = (init_rot_errors - prop_rot_errors) ./ init_rot_errors * 100;   % prop最后

% 显示统计摘要（保持prop在最后）
fprintf('\n===== 统计摘要 (%d次实验) =====\n', num_experiments);
fprintf('初始位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
    mean(init_pos_errors), std(init_pos_errors));
fprintf('初始旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
    mean(init_rot_errors), std(init_rot_errors));

fprintf('\nLM算法:\n');
fprintf('位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
    mean(lm_pos_errors), std(lm_pos_errors));
fprintf('旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
    mean(lm_rot_errors), std(lm_rot_errors));
fprintf('位置误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(pos_error_reduction_lm), std(pos_error_reduction_lm));
fprintf('旋转误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(rot_error_reduction_lm), std(rot_error_reduction_lm));

fprintf('\nELM算法:\n');
fprintf('位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
    mean(elm_pos_errors), std(elm_pos_errors));
fprintf('旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
    mean(elm_rot_errors), std(elm_rot_errors));
fprintf('位置误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(pos_error_reduction_elm), std(pos_error_reduction_elm));
fprintf('旋转误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(rot_error_reduction_elm), std(rot_error_reduction_elm));

fprintf('\nUnion算法:\n');  % 新增Union部分
fprintf('位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
    mean(union_pos_errors), std(union_pos_errors));
fprintf('旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
    mean(union_rot_errors), std(union_rot_errors));
fprintf('位置误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(pos_error_reduction_union), std(pos_error_reduction_union));
fprintf('旋转误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(rot_error_reduction_union), std(rot_error_reduction_union));

fprintf('\n所提算法:\n');  % prop保持最后输出
fprintf('位置误差: 平均 = %.6f m, 标准差 = %.6f m\n', ...
    mean(prop_pos_errors), std(prop_pos_errors));
fprintf('旋转误差: 平均 = %.6f, 标准差 = %.6f\n', ...
    mean(prop_rot_errors), std(prop_rot_errors));
fprintf('位置误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(pos_error_reduction_prop), std(pos_error_reduction_prop));
fprintf('旋转误差减少比例: 平均 = %.1f%%, 标准差 = %.1f%%\n', ...
    mean(rot_error_reduction_prop), std(rot_error_reduction_prop));
end