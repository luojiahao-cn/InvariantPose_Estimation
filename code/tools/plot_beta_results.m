clc,clear,close all
load('../results/beta_sensitive.mat');


tl = tiledlayout(1, 2);
set(gcf, 'color', 'w', 'units', 'centimeters', 'position', 2*[0, 0, 8.9, 3]);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
nexttile(tl, 1)

num_beta = batch_results(1).results{1}.num_beta;
colors = slanCM('gor', num_beta);

for i = 1:num_beta
    [~, deltaidx] = find(batch_results(1).results{1}.R_beta_history{i}.delta_history < 1e-5, 1, 'first');
    semilogx(batch_results(1).results{1}.R_beta_history{i}.eR_history, 'LineWidth', 2, 'Color', colors(i, :));
    hold on;
    scatter(deltaidx, batch_results(1).results{1}.R_beta_history{i}.eR_history(deltaidx),...
    100, 'Marker', 'x', 'markerEdgeColor', colors(i, :), 'LineWidth', 2, 'HandleVisibility', 'off');
end
h = legend('$\beta = 0$', '$\beta = 10^{-2}$', '$\beta = 1$', '$\beta = 10^1$', '$\beta = 10^2$', '$\beta = 10^3$',...
'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
h.ItemTokenSize = [10, 10];
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$e_{\mbox{\boldmath{$R$}}}$ [-]', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
ylim([0, 0.25])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'LineWidth', 1);

nexttile(tl, 2)
% 注意：不要使用clear，会清除之前加载的数据
load('../results/beta_robust.mat');

num_beta = batch_results(1).results{1}.num_beta;
num_points = size(batch_results, 2);
beta_vec = batch_results(1).results{1}.beta_vec;

% 将beta=0替换为1e-15以便显示
beta_vec_display = beta_vec;
beta_vec_display(beta_vec == 0) = 1e-15;

% 为每个beta值收集eR_vec数据
eR_matrix = zeros(num_points, num_beta); % 存储所有beta的eR_vec
iter_step_matrix = zeros(num_points, num_beta); % 存储迭代步数

for i = 1:num_beta
    for j = 1:num_points
        [eR_matrix(j, i), ind] = min(batch_results(j).results{1}.R_beta_history{i}.eR_history);
        iter_step_matrix(j, i) = ind;
    end
end

% 创建eR_matrix的箱线图
% figure;
boxplot(eR_matrix, 'Labels', {'1', '2', '3', '4', '5', '6', '7', '8'}, ...
'Labels', {'$0$', '$10^{-2}$', '$1$', '$10^1$', '$50$', '$10^2$', '$500$', '$10^3$'}, 'widths', 0.2);
% 设置boxplot的线条为黑色，内线为红色。box的T字形为实线

colors = slanCM('gor', 2);
% 获取箱子对象（外框）
boxes = findobj(gca, 'Tag', 'Box');
set(boxes, 'Color', colors(1, :));   % 箱体外边颜色改为绿色

% 获取中位数线对象
medians = findobj(gca, 'Tag', 'Median');
set(medians, 'Color', colors(2, :), 'LineWidth', 1.5);  % 中位线颜色设为红色

xlabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$e_{\mbox{\boldmath{$R$}}}$ [-]', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'LineWidth', 1);
grid on;