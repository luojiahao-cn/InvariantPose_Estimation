function plot_batch_results_2d(batch_results, ideal_path, params)
% PLOT_BATCH_RESULTS_2D 绘制批量实验结果 (2D 平面投影，仅绘制当前段)
% 输入：
%   batch_results - 由 run_batch_experiments 返回的批量结果结构体数组
%   ideal_path    - 3xN 矩阵，包含 NaN 间隔的理想坐标
%   params        - 实验参数结构体
close all;
if nargin < 3, params = struct(); end

%% 1. 数据提取与预处理
num_points = length(batch_results);
num_methods = 4;
method_names = {'ours', 'lm', 'elm', 'fischer'};

all_coords = [batch_results.test_point];
all_p_true = [all_coords.p_true]; 

% 自动分段逻辑
bin_indices = ones(1, num_points);
num_segments = 1;

try
    data_path = fullfile(fileparts(mfilename('fullpath')), '..', 'exp', 'path_ind.mat');
    if exist(data_path, 'file')
        data_meta = load(data_path, 'magnetic_char_lengths');
        if isfield(data_meta, 'magnetic_char_lengths')
            lengths = data_meta.magnetic_char_lengths;
            if sum(lengths) == num_points
                num_segments = length(lengths);
                bin_indices = zeros(1, num_points);
                curr_ptr = 1;
                for k = 1:num_segments
                    L = lengths(k);
                    bin_indices(curr_ptr:curr_ptr+L-1) = k;
                    curr_ptr = curr_ptr + L;
                end
            end
        end
    end
catch
end

% 计算全局误差用于颜色映射 (Colormap)
pos_error_matrix = zeros(num_points, num_methods);
rot_error_matrix = zeros(num_points, num_methods);
for i = 1:num_points
    for j = 1:num_methods
        method = method_names{j};
        pos_error_matrix(i, j) = batch_results(i).summary.(method).pos_mean;
        rot_error_matrix(i, j) = 2*asin(batch_results(i).summary.(method).rot_mean/(2*sqrt(2)));
    end
end

% 归一化综合误差到[0,1]
delta = 2;
rot_err_clamped = rot_error_matrix;
rot_err_clamped(rot_err_clamped > delta) = delta;
error_sum_matrix = pos_error_matrix + rot_err_clamped;
err_max = max(error_sum_matrix(:));
err_min = min(error_sum_matrix(:));
if err_max > err_min
    error_norm = (error_sum_matrix - err_min) / (err_max - err_min);
else
    error_norm = ones(size(error_sum_matrix)) * 0.5;
end

%% 2. 绘制分段图
if num_segments == 1
    render_one_figure('Overview', batch_results, ideal_path, error_norm, ...
                     pos_error_matrix, rot_error_matrix, bin_indices, [], params);
else
    mag_labels = 'MAGNETIC';
    for k = 1:num_segments
        if num_segments == 8
            tag = mag_labels(k);
        else
            tag = sprintf('Seg%d', k);
        end
        render_one_figure(tag, batch_results, ideal_path, error_norm, ...
                         pos_error_matrix, rot_error_matrix, bin_indices, k, params);
    end
end

end

%% ========================== 核心渲染子函数 ========================== %%
function render_one_figure(title_tag, batch_results, ideal_path, error_norm, ...
                          pos_err_mat, rot_err_mat, bin_indices, highlight_idx, params)
    
    num_points = length(batch_results);
    num_methods = 4;
    method_labels = {'Ours', '(1)', '(4)', 'Fischer'};
    all_coords = [batch_results.test_point];
    all_p_true = [all_coords.p_true];
    fontsize = 16;
    n_colors = 256;
    custom_colormap = slanCM('gor', n_colors);

    % 获取传感器个数
    n_sensors = 0;
    if isfield(params, 'sensor') && isfield(params.sensor, 'd_list')
        n_sensors = size(params.sensor.d_list, 2);
    elseif isfield(params, 'd_list')
        n_sensors = size(params.d_list, 2);
    end

    fig_name = sprintf('Results 2D: %s', title_tag);
    % 2D 布局通常较宽
    fig_handle = figure('Color', 'white', 'units', 'centimeters', 'Position', [2, 2, 18.2, 7], 'Name', fig_name);
    tl = tiledlayout(2, num_methods, 'TileSpacing', 'compact', 'Padding', 'tight');
    
    % --- 绘制 2D 轨迹 (XY 平面) ---
    for j = 1:num_methods
        nexttile(j);
        hold on; grid on; box on;
        
        % 1. 绘制理想路径 (仅绘制当前段)
        if ~isempty(ideal_path) && ~isempty(highlight_idx)
            nan_locs = find(isnan(ideal_path(1, :)));
            start_i = 1;
            for s = 1:(length(nan_locs) + 1)
                if s <= length(nan_locs), end_i = nan_locs(s) - 1; else, end_i = size(ideal_path, 2); end
                if s == highlight_idx && end_i >= start_i
                    seg = ideal_path(:, start_i:end_i);
                    plot(seg(2, :), seg(3, :), '-', 'Color', [0 0.4470 0.7410 0.5], 'LineWidth', 2.5);
                end
                start_i = end_i + 2;
            end
        end
        
        % 2. 筛选数据
        if ~isempty(highlight_idx)
            mask = (bin_indices == highlight_idx);
        else
            mask = true(1, num_points);
        end
        
        % 提取点颜色
        colors = zeros(sum(mask), 3);
        subset_err_norm = error_norm(mask, j);
        for i = 1:length(subset_err_norm)
            c_idx = max(1, min(n_colors, round(subset_err_norm(i) * (n_colors - 1)) + 1));
            colors(i, :) = custom_colormap(c_idx, :);
        end
        
        % % 3. 叠加估计点云 (Jitter Cloud)
        % all_p_est = zeros(3, num_points);
        method_keys_map = {'ours', 'lm', 'elm', 'fischer'};
        % curr_method = method_keys_map{j};
        % for i = find(mask)
        %     all_p_est(:, i) = batch_results(i).results(1).(['p_' curr_method]);
        % end
        % % 使用彩色半透明散点表示抖动
        % scatter(all_p_est(2, mask)', all_p_est(3, mask)', 25, colors, 'filled', ...
        %         'MarkerFaceAlpha', 0.4);

        % 4. 绘制真实位置 (加大点径并添加黑边，确保可见)
        scatter(all_p_true(2, mask)', all_p_true(3, mask)', 25, colors, 'filled', ...
                'LineWidth', 0.8);
        
        % 2D 视图设置
        xlabel('$y$ [m]', 'FontSize', fontsize-2, 'Interpreter', 'latex');
        if j == 1
            ylabel('$z$ [m]', 'FontSize', fontsize-2, 'Interpreter', 'latex');
            yticks([-0.1, 0, 0.1]);
            text(-0.96, -0.33, ['$x=$', num2str(seg(1, 1),1), '[m]'], 'fontsize', fontsize, 'Interpreter', 'latex','units','normalized');
        else
            yticklabels([]); % 保持紧凑，只显示首列刻度
        end
        ylim([-0.12, 0.12]);
        axis equal; grid minor;
        set(gca, 'FontSize', fontsize-2, 'TickLabelInterpreter', 'latex', 'TickDir', 'out');
        title(method_labels{j}, 'FontSize', fontsize, 'Interpreter', 'latex');
        
        % 在顶部中心标注当前 X 深度 (告诉读者这是 3D 结构的切面)
        % if ~isempty(x_val_str)
        %     text(0.5, 0.92, x_val_str, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        %          'FontSize', fontsize-4, 'Interpreter', 'latex', 'FontWeight', 'bold');
        % end
    end
    
    % --- 绘制误差箱线图 ---
    for j = 1:num_methods
        subtl = tiledlayout(tl, 2, 1);
        subtl.Layout.Tile = num_methods + j;
        
        % 筛选当前段数据
        if ~isempty(highlight_idx)
            mask = (bin_indices == highlight_idx);
        else
            mask = true(1, num_points);
        end
        
        % R 误差
        nexttile(subtl, 1);
        r_errs = rot_err_mat(mask, j);
        gc_r = categorical(repmat("Rotation", length(r_errs), 1));
        boxchart(gc_r, r_errs, 'BoxFaceColor', [0.2 0.2 0.2], 'Orientation', 'horizontal');
        if j == 1, yticklabels('$e_{\mathbf{R}}$ [rad]'); else, yticklabels(''); end
        set(gca, 'FontSize', fontsize-2, 'TickDir', 'out', 'TickLabelInterpreter', 'latex'); box on;
        
        % p 误差
        nexttile(subtl, 2);
        p_errs = 1e3 * pos_err_mat(mask, j); % mm
        gc_p = categorical(repmat("Position", length(p_errs), 1));
        boxchart(gc_p, p_errs, 'BoxFaceColor', [0.2 0.2 0.2], 'Orientation', 'horizontal');
        xlim([0, 10]);
        if j == 1, yticklabels('$e_{\mathbf{p}}$ [mm]'); else, yticklabels(''); end
        set(gca, 'FontSize', fontsize-2, 'TickDir', 'out', 'TickLabelInterpreter', 'latex'); box on;
    end

    % --- 导出图像 (已按要求注释) ---
    if ~exist('results', 'dir'), mkdir('results'); end
    filename = sprintf('results/PoseEst_2D_%s_Ns%d.png', title_tag, n_sensors);
    try
        export_fig(fig_handle, filename, '-nocrop', '-r600', '-transparent');
        fprintf('[render_one_figure] 2D 版图像已导出: %s\n', filename);
    catch ME
        saveas(fig_handle, filename);
    end
end
