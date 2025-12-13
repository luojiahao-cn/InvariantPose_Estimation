close all; clear;

fontsize = 16;

% ===================== Figure & Layout =====================
figure;
set(gcf, 'color', 'w', 'units', 'centimeters', 'position', [10, 10, 19.8, 12]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ===================== Common Plot Params =====================
group_spacing = 3.0;     % 大组之间间距（Config 之间）
within_offset = 0.35;    % 同一大组内 7 个箱线展开宽度
box_width     = 0.22;

% ===================== Top: gamma -> pos_error =====================
load('../../results/batch_W.mat', 'batch_results');

num_dlist = length(batch_results);
num_W     = length(batch_results{1});

pos_error = zeros(num_dlist, num_W, 300);
for d = 1:num_dlist
    num_W_i = length(batch_results{d});
    for w = 1:num_W_i
        pos_error(d, w, :) = [batch_results{d}(w).results.ours_pos_error]';
    end
end

[num_dlist, num_W, N] = size(pos_error);
colors = lines(num_dlist);

ax1 = nexttile; hold(ax1, 'on');

for d = 1:num_dlist
    for w = 1:num_W
        data = squeeze(pos_error(d, w, :)); % 300x1
        x = d * group_spacing + (w - (num_W+1)/2) * within_offset;

        boxchart(ax1, x * ones(size(data)), data, ...
            'BoxWidth', box_width, ...
            'BoxFaceColor', colors(d,:), ...
            'MarkerStyle', 'none');
    end
end

xlabel(ax1, '$\log_{10}(\gamma)$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel(ax1, '$e_{\mbox{\boldmath{$p$}}}$ [m]', 'Interpreter', 'latex', 'FontSize', fontsize);
set(ax1, 'FontSize', fontsize, 'LineWidth', 1);
grid(ax1, 'on');
box on;
% 28 ticks for gamma
gamma_log10 = [-5 -4 -3 -2 -1 0 1];

all_x = zeros(num_dlist*num_W, 1);
all_labels = strings(num_dlist*num_W, 1);
idx = 1;
for d = 1:num_dlist
    for w = 1:num_W
        all_x(idx) = d * group_spacing + (w - (num_W+1)/2) * within_offset;
        all_labels(idx) = sprintf('%d', gamma_log10(w));
        idx = idx + 1;
    end
end

xticks(ax1, all_x);
xticklabels(ax1, all_labels);
set(ax1, 'TickLabelInterpreter', 'latex');
xlim(ax1, [min(all_x)-within_offset, max(all_x)+within_offset]);

% Legend handles (create once)
h = gobjects(num_dlist, 1);
for d = 1:num_dlist
    h(d) = plot(ax1, nan, nan, 's', ...
        'MarkerFaceColor', colors(d,:), ...
        'MarkerEdgeColor', colors(d,:));
end

% ===================== Bottom: mu -> rot_error =====================
ax2 = nexttile; hold(ax2, 'on');

load('../../results/batch_mu.mat', 'batch_results');

num_dlist = length(batch_results);
num_W     = length(batch_results{1});

rot_error = zeros(num_dlist, num_W, 300);
for d = 1:num_dlist
    num_W_i = length(batch_results{d});
    for w = 1:num_W_i
        rot_error(d, w, :) = [batch_results{d}(w).results.ours_rot_error]';
    end
end

[num_dlist, num_W, N] = size(rot_error);
colors = lines(num_dlist); % 保持与上图一致

for d = 1:num_dlist
    for w = 1:num_W
        data = squeeze(rot_error(d, w, :)); % 300x1
        x = d * group_spacing + (w - (num_W+1)/2) * within_offset;

        boxchart(ax2, x * ones(size(data)), data, ...
            'BoxWidth', box_width, ...
            'BoxFaceColor', colors(d,:), ...
            'MarkerStyle', 'none');
    end
end

xlabel(ax2, '$\log_{10}(\mu)$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel(ax2, '$e_{\mbox{\boldmath{$R$}}}$ [rad]', 'Interpreter', 'latex', 'FontSize', fontsize);
set(ax2, 'FontSize', fontsize, 'LineWidth', 1);
grid(ax2, 'on');
box on;
% 28 ticks for mu
mu_log10 = [-3 -2 -1 0 1 2 3];

all_x = zeros(num_dlist*num_W, 1);
all_labels = strings(num_dlist*num_W, 1);
idx = 1;
for d = 1:num_dlist
    for w = 1:num_W
        all_x(idx) = d * group_spacing + (w - (num_W+1)/2) * within_offset;
        all_labels(idx) = sprintf('%d', mu_log10(w));
        idx = idx + 1;
    end
end

xticks(ax2, all_x);
xticklabels(ax2, all_labels);
set(ax2, 'TickLabelInterpreter', 'latex');
xlim(ax2, [min(all_x)-within_offset, max(all_x)+within_offset]);
ylim(ax2, [0, 1]);

% ========================================================================
%                     Shared Legend: one per Config (1x4)
% ========================================================================
% 用“代理句柄”确保 legend 每个 Config 只出现一次，并且颜色正确
h = gobjects(num_dlist, 1);
for d = 1:num_dlist
    h(d) = plot(ax1, nan, nan, 's', ...
        'MarkerFaceColor', colors(d,:), ...
        'MarkerEdgeColor', colors(d,:), ...
        'MarkerSize', 8);
end

lgd = legend(h, {'Redundant configuration 1', 'Redundant configuration 2', 'Minimal configuration 1', 'Minimal configuration 2'}, ...
    'Interpreter', 'latex', 'FontSize', fontsize-2);
lgd.NumColumns = 2;
lgd.Orientation = 'horizontal';
lgd.Layout.Tile = 'north';