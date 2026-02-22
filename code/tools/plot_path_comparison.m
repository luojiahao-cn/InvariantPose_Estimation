function plot_path_comparison(batch_results, params)
% PLOT_PATH_COMPARISON 绘制路径估计结果对比，重点对比 Ours (SSL) 和其他经典算法
% 输入：
%   batch_results - 批量实验结果 (如果为空则尝试加载 results/path_cubic_nested.mat)
%   params        - 实验参数

%% 1. 数据加载与预处理
if nargin < 1 || isempty(batch_results)
    % data_path = './results/path_cubic_nested_2.mat';
    data_path = '../results/point_test_cone.mat';
    if exist(data_path, 'file')
        fprintf('[plot_path_comparison] 正在从 %s 加载数据...\n', data_path);
        loaded = load(data_path);
        batch_results = loaded.batch_results;
        params = loaded.params;
    else
        error('[plot_path_comparison] 未提供数据且未找到 %s', data_path);
    end
end

if nargin < 2 || isempty(params)
    params = struct();
end

num_points = length(batch_results);
method_tags = {'ours', 'lm', 'elm', 'fischer'};
method_labels = {'Ours (SSL)', '(3)', '(6)', 'Fischer'};
num_methods = length(method_tags);

% 提取真实轨迹
all_p_true = zeros(3, num_points);
for i = 1:num_points
    all_p_true(:, i) = batch_results(i).test_point.p_true;
end

% 提取位置误差和主轴指向误差 (1 - |r'*r_hat|)
pos_err_mat = zeros(num_points, num_methods);
rot_err_mat = zeros(num_points, num_methods);

for i = 1:num_points
    s = batch_results(i).summary;
    for j = 1:num_methods
        tag = method_tags{j};
        
        % 位置误差
        if isfield(s, tag)
            pos_err_mat(i, j) = s.(tag).pos_mean;
        else
            pos_err_mat(i, j) = NaN;
        end
        
        % 指向误差 (r-axis error)
        if strcmp(tag, 'ours')
             if isfield(s.ours, 'direct_r_mean_SSL')
                rot_err_mat(i, j) = s.ours.direct_r_mean_SSL;
             else
                rot_err_mat(i, j) = s.ours.r_mean;
             end
        else
            if isfield(s, tag) && isfield(s.(tag), 'r_mean')
                rot_err_mat(i, j) = s.(tag).r_mean;
            else
                rot_err_mat(i, j) = NaN;
            end
        end
    end
end

%% 2. 误差融合与颜色映射 (Early Plotting Strategy)
% 我们将位置误差 (m) 和指向误差 (无量纲 0-1) 结合
% 通常位置误差在 0-0.05m 左右，指向误差在 0-0.2 左右
% 这里采用简单的加权融合或分别截断归一化
max_p_err = 0.1; % 2cm 截断
max_r_err = 0.1;  % 0.1 截断 (约 25度)
% pos_err_mat(pos_err_mat > max_p_err) = max_p_err;
% rot_err_mat(rot_err_mat > max_r_err) = max_r_err;

% 采用固定量程归一化 [0, 0.2]，确保不同实验图之间的颜色具有可比性
% 0.2 对应 Colorbar 的最大值
error_combined = pos_err_mat + rot_err_mat; % 简单加权融合 Riemannian metric

error_combined = (error_combined - min(error_combined(:))) / (max(error_combined(:)) - min(error_combined(:)) + eps);

% Colormap: 绿 -> 黄 -> 红
n_colors = 256;
try
    % 如果有 slanCM 扩展包
    custom_colormap = slanCM('gor', n_colors);
catch
    % 自定义绿-红渐变
    c1 = [0, 1, 0]; % 绿
    c2 = [1, 1, 0]; % 黄
    c3 = [1, 0, 0]; % 红
    custom_colormap = [interp1([0, 0.5, 1], [c1; c2; c3], linspace(0, 1, n_colors))];
end

%% 3. 绘图渲染
fig_handle = figure('Color', 'white', 'Name', 'Path Estimation Comparison', ...
                   'Units', 'centimeters', 'Position', [2, 2, 33, 12.03]);
tiledlayout(5, num_methods, 'Padding', 'compact', 'TileSpacing', 'compact');

fontsize = 16;

for j = 1:num_methods
    % --- 绘制 3D 轨迹图 (占用前 3 行) ---
    nexttile(j, [3, 1]);
    hold on; grid on; box on;
    
    % 绘制串联线 (背景线)
    plot3(all_p_true(1, :), all_p_true(2, :), all_p_true(3, :), '-', ...
          'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    
    % 绘制轨迹点，根据误差着色
    for i = 1:num_points
        c_idx = max(1, min(n_colors, round(error_combined(i, j) * (n_colors - 1)) + 1));
        scatter3(all_p_true(1, i), all_p_true(2, i), all_p_true(3, i), ...
                 30, custom_colormap(c_idx, :), 'filled', 'MarkerFaceAlpha', 1);
    end
    
    % 设置视角
    view(-14, 7);
    axis equal;

    xlim([-0.1, 0.1]); ylim([-0.1, 0.1]); zlim([-0.1, 0.1]);

    
    % 标签与标题
    if j == 1
        xlabel('$x$ [m]', 'Interpreter', 'latex');
        ylabel('$y$ [m]', 'Interpreter', 'latex');
        zlabel('$z$ [m]', 'Interpreter', 'latex');
    end
    
    title(method_labels{j}, 'FontSize', fontsize + 2, 'Interpreter', 'latex');
    set(gca, 'FontSize', fontsize, 'TickLabelInterpreter', 'latex', 'TickDir', 'out');
end

% --- 绘制误差箱线图 (ep 和 er) ---
for j = 1:num_methods
    % Position Error Boxchart (Row 4)
    nexttile(3 * num_methods + j);
    ep_data = pos_err_mat(:, j) * 1e3;
    boxchart(categorical(repmat("ep", length(ep_data), 1)), ep_data, ...
             'BoxFaceColor', [0.2 0.2 0.2], 'Orientation', 'horizontal');
    grid on; box on;
    if j == 1, yticklabels('$e_{\mbox{\boldmath{$p$}}}$ [mm]'); else, yticklabels(''); end
    set(gca, 'FontSize', fontsize, 'TickLabelInterpreter', 'latex', 'TickDir', 'out');

    xlim([0, 1]);

    % Rotation Error Boxchart (Row 5)
    nexttile(4 * num_methods + j);
    er_data = rot_err_mat(:, j);
    boxchart(categorical(repmat("er", length(er_data), 1)), er_data, ...
             'BoxFaceColor', [0.4 0.4 0.4], 'Orientation', 'horizontal');
    grid on; box on;
    if j == 1, yticklabels('$e_{\mbox{\boldmath{$r$}}}$ [rad]'); else, yticklabels(''); end
    set(gca, 'FontSize', fontsize, 'TickLabelInterpreter', 'latex', 'TickDir', 'out');

    xlim([0, 0.1]);
end

% 全局 Colorbar - 置于最下方横向显示
cb = colorbar;
cb.Layout.Tile = 'south';
cb.Orientation = 'horizontal';
% cb.Label.String = '';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fontsize;
cb.Ticks = [0, 0.5, 1];
cb.TickLabels = {'min', '$e_{\mbox{\boldmath{$p$}}} + e_{\mbox{\boldmath{$r$}}}$', 'max'};
cb.TickLabelInterpreter = 'latex';
colormap(custom_colormap);

% 保存结果
if ~exist('results', 'dir'), mkdir('results'); end

% 获取传感器数量用于命名
num_sensors = 0;
if isfield(params, 'num_sensors')
    num_sensors = params.num_sensors;
elseif isfield(params, 'sensor') && isfield(params.sensor, 'd_list')
    num_sensors = size(params.sensor.d_list, 2);
end

if num_sensors > 0
    filename = sprintf('results/Path_Comparison_N%d.png', num_sensors);
else
    filename = 'results/Path_Comparison_Comparison.png';
end

saveas(fig_handle, filename);
fprintf('[plot_path_comparison] 图像已保存至: %s\n', filename);

end
