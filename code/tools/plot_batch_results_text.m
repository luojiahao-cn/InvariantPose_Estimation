function plot_batch_results_text(batch_results, ideal_path, params)
% PLOT_BATCH_RESULTS_TEXT 绘制批量实验结果 (误差信息以文本形式置于标题中)
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
                fprintf('[plot_batch_results_text] 成功按元数据分组: %s\n', mat2str(lengths));
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

% 归一化综合误差到[0,1] 用于颜色显示
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
    fontsize = 14;
    n_colors = 256;
    custom_colormap = slanCM('gor', n_colors);

    % 获取传感器个数
    n_sensors = 0;
    if isfield(params, 'sensor') && isfield(params.sensor, 'd_list')
        n_sensors = size(params.sensor.d_list, 2);
    elseif isfield(params, 'd_list')
        n_sensors = size(params.d_list, 2);
    end

    fig_name = sprintf('Results: %s', title_tag);
    % 减少高度，因为去掉了 boxcharts
    fig_handle = figure('Color', 'white', 'units', 'centimeters', 'Position', [2, 2, 28, 7.5], 'Name', fig_name);
    tl = tiledlayout(1, num_methods, 'TileSpacing', 'compact', 'Padding', 'tight');
    
    % --- 绘制 3D 散点图 ---
    for j = 1:num_methods
        nexttile;
        hold on; grid on; box on;
        
        % 1. 绘制空间参考
        [gx, gy, gz] = meshgrid(linspace(-0.1, 0.1, 10));
        scatter3(gx(:), gy(:), gz(:), 2, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.05);

        % 2. 绘制理想路径
        if ~isempty(ideal_path)
            nan_locs = find(isnan(ideal_path(1, :)));
            start_i = 1;
            for s = 1:(length(nan_locs) + 1)
                if s <= length(nan_locs), end_i = nan_locs(s) - 1; else, end_i = size(ideal_path, 2); end
                if end_i >= start_i
                    seg = ideal_path(:, start_i:end_i);
                    is_main = (~isempty(highlight_idx) && s == highlight_idx);
                    if is_main
                        l_color = [0 0.4470 0.7410 0.6];
                        l_width = 1.8;
                    else
                        l_color = [0.8 0.8 0.8 0.1];
                        l_width = 1.0;
                    end
                    plot3(seg(1, :), seg(2, :), seg(3, :), '-', 'Color', l_color, 'LineWidth', l_width);
                end
                start_i = end_i + 2;
            end
        end
        
        % 3. 筛选与准备数据
        if ~isempty(highlight_idx)
            mask = (bin_indices == highlight_idx);
        else
            mask = true(1, num_points);
        end
        
        % 计算当前段的统计信息 (ep 取 mm, er 取 rad)
        curr_p_errs_mm = 1e3 * pos_err_mat(mask, j);
        curr_r_errs_rad = rot_err_mat(mask, j);
        
        mu_p = mean(curr_p_errs_mm);
        std_p = std(curr_p_errs_mm);
        mu_r = mean(curr_r_errs_rad);
        std_r = std(curr_r_errs_rad);

        colors = zeros(num_points, 3);
        alphas = ones(num_points, 1);
        for i = 1:num_points
            if ~isempty(highlight_idx) && bin_indices(i) ~= highlight_idx
                colors(i, :) = [0.85, 0.85, 0.85];
                alphas(i) = 0.1;
            else
                c_idx = max(1, min(n_colors, round(error_norm(i, j) * (n_colors - 1)) + 1));
                colors(i, :) = custom_colormap(c_idx, :);
                alphas(i) = 1.0;
            end
        end
        
        % 4. 叠加估计 jitter
        all_p_est = zeros(3, num_points);
        method_keys_map = {'ours', 'lm', 'elm', 'fischer'};
        curr_method = method_keys_map{j};
        for i = 1:num_points
            all_p_est(:, i) = batch_results(i).results(1).(['p_' curr_method]);
        end
        scatter3(all_p_est(1, :), all_p_est(2, :), all_p_est(3, :), 8, colors, 'filled', ...
                 'MarkerFaceAlpha', 'flat', 'AlphaData', alphas * 0.05, 'MarkerEdgeAlpha', 0);

        % 5. 真实位置
        scatter3(all_p_true(1, :), all_p_true(2, :), all_p_true(3, :), 40, colors, 'filled', ...
                 'MarkerFaceAlpha', 'flat', 'MarkerEdgeAlpha', 'flat', 'AlphaData', alphas);
        
        % 视图设置
        xlim([-0.12, 0.12]); ylim([-0.12, 0.12]); zlim([-0.12, 0.12]);
        xticks(-0.1:0.1:0.1); yticks(-0.1:0.1:0.1); zticks(-0.1:0.1:0.1);
        xlabel('$x$ [m]', 'Interpreter', 'latex'); ylabel('$y$ [m]', 'Interpreter', 'latex'); zlabel('$z$ [m]', 'Interpreter', 'latex');
        axis square; view(122, 13);
        set(gca, 'FontSize', fontsize-2, 'TickLabelInterpreter', 'latex');
        
        % --- 关键修改：将 ep 和 er 统计信息放入题目 ---
        % 使用 LaTeX 格式: ep: mu ± std mm, er: mu ± std rad
        title_str = sprintf('%s\\newline\\fontsize{10}{12}\\selectfont $e_p: %.2f \\pm %.2f$ mm, $e_R: %.2f \\pm %.2f$ rad', ...
                            method_labels{j}, mu_p, std_p, mu_r, std_r);
        title(title_str, 'FontSize', fontsize, 'Interpreter', 'latex');
    end

    % --- 导出图像 ---
    if ~exist('results', 'dir'), mkdir('results'); end
    filename = sprintf('results/PoseEst_Text_%s_Ns%d.png', title_tag, n_sensors);
    try
        export_fig(fig_handle, filename, '-nocrop', '-r600', '-transparent');
        fprintf('[render_one_figure] 文本标注版图像已导出: %s\n', filename);
    catch
        saveas(fig_handle, filename);
    end
end
