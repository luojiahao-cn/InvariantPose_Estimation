%% 误差带图（均值±标准差）
clc,clear,close all;
load ../results/expA_main.mat

red = [255, 47, 23]/255;
blue = [33, 77, 169]/255;
green = [39, 132, 60]/255;
alg_colors = [red; blue; green; 0,0,0];
labels = {'Free', 'Disturbance', 'Noise', 'Disturbance and noise'};

figure('Name','位置误差带');
tiledlayout(1, 4, 'Padding','compact', 'TileSpacing','compact');
set(gcf, 'Units', 'centimeters', 'Position', 2*[0, 0, 17.8, 3], 'color', 'w');
for mode = 1:num_cond
    nexttile();
    hold on;
    h = gobjects(num_alg,1);
    for i = 1:num_alg
        mu = mean(squeeze(err_pos(:,mode,i,:)),2);
        sigma = std(squeeze(err_pos(:,mode,i,:)),0,2);
        fill([offset_list, fliplr(offset_list)], ...
             [mu+sigma; flipud(mu-sigma)]', ...
             alg_colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor','none');
        h(i) = plot(offset_list, mu, 'Color', alg_colors(i,:), 'LineWidth',1.5);
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

figure('Name','旋转误差带','Color','white');
tiledlayout(1, 4, 'Padding','compact', 'TileSpacing','compact');
for mode = 1:num_cond
    nexttile();
    hold on;
    h = gobjects(num_alg,1);
    for i = 1:num_alg
        mu = mean(squeeze(err_rot(:,mode,i,:)),2);
        sigma = std(squeeze(err_rot(:,mode,i,:)),0,2);
        fill([offset_list, fliplr(offset_list)], ...
             [mu+sigma; flipud(mu-sigma)]', ...
             alg_colors(i,:), 'FaceAlpha',0.15, 'EdgeColor','none');
        h(i) = plot(offset_list, mu, 'Color', alg_colors(i,:), 'LineWidth',2);
    end
    xlabel('p初值偏离度 (m)');
    ylabel('旋转误差 (Frobenius)');
    title(labels{mode});
    if mode == 1
        legend(h, alg_labels, 'Location','northwest');
    end
    grid on;
    hold off;
end