%% 误差带图（均值±标准差）
clc,clear,close all;
load ../results/expA_main.mat

red = [255, 47, 23]/255;
blue = [33, 77, 169]/255;
green = [39, 132, 60]/255;
black = [0, 0, 0]/255;
alg_labels = {'LM','ELM','Ours'};
alg_colors = [red; blue; black]; % 去掉green

figure('Name','位置误差带');
tiledlayout(1, 4, 'Padding','compact', 'TileSpacing','compact');
set(gcf, 'Units', 'centimeters', 'Position', 2*[0, 0, 17.8, 3], 'color', 'w');
for mode = 1:num_cond
    nexttile();
    hold on;
    h = gobjects(numel(alg_labels),1);
    idx = 1;
    for i = [1,2,4] % 只显示LM, ELM, Ours
        mu = mean(squeeze(err_pos(:,mode,i,:)),2);
        sigma = std(squeeze(err_pos(:,mode,i,:)),0,2);
        fill([offset_list, fliplr(offset_list)], ...
             [mu+sigma; flipud(mu-sigma)]', ...
             alg_colors(idx,:), 'FaceAlpha', 0.2, 'EdgeColor','none');
        h(idx) = plot(offset_list, mu, 'Color', alg_colors(idx,:), 'LineWidth',1.5);
        idx = idx + 1;
    end
    xlim([0, 0.05]); xticks(0:0.025:0.05);
    % ylim([0, 0.05]); yticks([0, 0.05]);
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    xlabel('Uncertainty strength of $\mbox{\boldmath{$p$}}_0$ [m]', 'Interpreter','latex', 'fontsize', 12);
    ylabel('$e_{\mbox{\boldmath{$p$}}}$ [m]', 'Interpreter','latex', 'fontsize', 12);
    title(labels{mode}, 'Interpreter','latex', 'fontsize', 14);
    if mode == 4
        legend(h, alg_labels, 'Location','northeastoutside', 'Interpreter','latex', 'fontsize', 12, 'box', 'off');
    end
    grid on;
    % box on;
    hold off;
end

figure('Name','旋转误差带');
tiledlayout(1, 4, 'Padding','compact', 'TileSpacing','compact');
set(gcf, 'Units', 'centimeters', 'Position', 2*[0, 0, 17.8, 3], 'color', 'w');
for mode = 1:num_cond
    nexttile();
    hold on;
    h = gobjects(numel(alg_labels),1);
    idx = 1;
    for i = [1,2,4] % 只显示LM, ELM, Ours
        mu = mean(squeeze(err_rot(:,mode,i,:)),2);
        sigma = std(squeeze(err_rot(:,mode,i,:)),0,2);
        fill([offset_list, fliplr(offset_list)], ...
             [mu+sigma; flipud(mu-sigma)]', ...
             alg_colors(idx,:), 'FaceAlpha', 0.2, 'EdgeColor','none');
        h(idx) = plot(offset_list, mu, 'Color', alg_colors(idx,:), 'LineWidth',1.5);
        idx = idx + 1;
    end
    xlim([0, 0.05]); xticks(0:0.025:0.05);
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
    xlabel('Uncertainty strength of $\mbox{\boldmath{$p$}}_0$ [m]', 'Interpreter','latex', 'fontsize', 12);
    ylabel('$e_{R}$ (Frobenius)', 'Interpreter','latex', 'fontsize', 12);
    title(labels{mode}, 'Interpreter','latex', 'fontsize', 14);
    if mode == 4
        legend(h, alg_labels, 'Location','northeastoutside', 'Interpreter','latex', 'fontsize', 12, 'box', 'off');
    end
    grid on;
    hold off;
end