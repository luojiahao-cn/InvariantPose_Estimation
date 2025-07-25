function plot_error_distributions(results)
% 改进的误差分布可视化函数

% 提取数据（增加Union算法，保持prop在最后）
lm_pos = [results.lm_pos_error];
elm_pos = [results.elm_pos_error];
union_pos = [results.union_pos_error];  % 新增Union位置误差
prop_pos = [results.prop_pos_error];   % prop保持最后

lm_rot = [results.lm_rot_error];
elm_rot = [results.elm_rot_error];
union_rot = [results.union_rot_error];  % 新增Union旋转误差
prop_rot = [results.prop_rot_error];    % prop保持最后

% 创建新图表
figure('Color', 'white', 'Position', [100, 100, 1200, 800], 'Name', '改进的误差分布');

% 位置误差分布 - 使用KDE（保持prop在最后）
subplot(2, 2, 1);
hold on;
[f_lm, x_lm] = ksdensity(lm_pos);
[f_elm, x_elm] = ksdensity(elm_pos);
[f_union, x_union] = ksdensity(union_pos);  % Union
[f_prop, x_prop] = ksdensity(prop_pos);   % prop最后

% 绘制顺序：LM -> ELM -> Union -> Proposed
plot(x_lm, f_lm, 'r-', 'LineWidth', 2);
plot(x_elm, f_elm, 'g-', 'LineWidth', 2);
plot(x_union, f_union, 'm-', 'LineWidth', 2);  % 紫色表示Union算法
plot(x_prop, f_prop, 'b-', 'LineWidth', 2);   % 蓝色表示所提算法（最后）

title('位置误差分布 (KDE)', 'FontSize', 14);
xlabel('位置误差 (m)', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend({'LM算法', 'ELM算法', 'Union算法', '所提算法'}, 'FontSize', 10, 'Location', 'best');  % 更新图例
grid on;
box on;

% 旋转误差分布 - 使用KDE（保持prop在最后）
subplot(2, 2, 2);
hold on;
[f_lm, x_lm] = ksdensity(lm_rot);
[f_elm, x_elm] = ksdensity(elm_rot);
[f_union, x_union] = ksdensity(union_rot);  % Union
[f_prop, x_prop] = ksdensity(prop_rot);   % prop最后

% 绘制顺序：LM -> ELM -> Union -> Proposed
plot(x_lm, f_lm, 'r-', 'LineWidth', 2);
plot(x_elm, f_elm, 'g-', 'LineWidth', 2);
plot(x_union, f_union, 'm-', 'LineWidth', 2);  % 紫色表示Union算法
plot(x_prop, f_prop, 'b-', 'LineWidth', 2);   % 蓝色表示所提算法（最后）

title('旋转误差分布 (KDE)', 'FontSize', 14);
xlabel('旋转误差', 'FontSize', 12);
ylabel('概率密度', 'FontSize', 12);
legend({'LM算法', 'ELM算法', 'Union算法', '所提算法'}, 'FontSize', 10, 'Location', 'best');  % 更新图例
grid on;
box on;

% 位置误差箱线图（保持prop在最后）
subplot(2, 2, 3);
boxplot([lm_pos', elm_pos', union_pos', prop_pos'], ...  % prop最后
    'Labels', {'LM', 'ELM', 'Union', 'Proposed'});       % 标签顺序匹配
title('位置误差分布比较', 'FontSize', 14);
ylabel('位置误差 (m)', 'FontSize', 12);
grid on;

% 旋转误差箱线图（保持prop在最后）
subplot(2, 2, 4);
boxplot([lm_rot', elm_rot', union_rot', prop_rot'], ...  % prop最后
    'Labels', {'LM', 'ELM', 'Union', 'Proposed'});       % 标签顺序匹配
title('旋转误差分布比较', 'FontSize', 14);
ylabel('旋转误差', 'FontSize', 12);
grid on;

% 添加整体标题
sgtitle('算法误差分布比较', 'FontSize', 16, 'FontWeight', 'bold');
end