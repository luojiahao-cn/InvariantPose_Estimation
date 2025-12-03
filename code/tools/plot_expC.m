function plot_expC()
% plot_expC.m
% 可视化 expC_main.mat 的定位误差结果
% 绘制4种算法在4类工况下的位置误差和旋转误差箱线图
% 结果图片自动保存到 figures/ 目录

close all;
S = load('../results/expC_main.mat');

alg_labels = string(S.alg_labels);
labels = string(S.labels);
num_cond = S.num_cond;
err_pos = S.err_pos;
err_rot = S.err_rot;

figdir = '../figures/';
if ~exist(figdir, 'dir'), mkdir(figdir); end

% 位置误差箱线图
figure('Name','位置误差','Color','w');
for mode = 1:num_cond
    subplot(2,2,mode);
    boxplot(squeeze(err_pos(mode,:,:))', alg_labels);
    title(labels(mode));
    ylabel('位置误差 (m)');
    grid on;
end
sgtitle('多磁铁实验 - 位置误差分布');
saveas(gcf, fullfile(figdir, 'expC_pos_error_box.png'));

% 旋转误差箱线图
figure('Name','旋转误差','Color','w');
for mode = 1:num_cond
    subplot(2,2,mode);
    boxplot(squeeze(err_rot(mode,:,:))', alg_labels);
    title(labels(mode));
    ylabel('旋转误差 (Frobenius)');
    grid on;
end
sgtitle('多磁铁实验 - 旋转误差分布');
saveas(gcf, fullfile(figdir, 'expC_rot_error_box.png'));

fprintf('已保存箱线图到 %s\n', figdir);
end
