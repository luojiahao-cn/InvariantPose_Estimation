function plot_principal_axis_results_cone(batch_results, params)
% PLOT_PRINCIPAL_AXIS_RESULTS_CONE 绘制主轴估计结果对比，重点对比 Ours (SSL) 和其他经典算法
% 输入：
%   batch_results - 批量实验结果 (如果为空则尝试加载 results/path_cubic_nested.mat)
%   params        - 实验参数

%% 1. 数据加载与预处理
% if nargin < 1 || isempty(batch_results)
%     % data_path = './results/path_cubic_nested_2.mat';
%     data_path = '../results/point_test_cone.mat';
%     if exist(data_path, 'file')
%         fprintf('[plot_principal_axis_results_cone] 正在从 %s 加载数据...\n', data_path);
%         loaded = load(data_path);
%         batch_results = loaded.batch_results;
%         params = loaded.params;
%     else
%         error('[plot_principal_axis_results_cone] 未提供数据且未找到 %s', data_path);
%     end
% end

% if nargin < 2 || isempty(params)
%     params = struct();
% end

num_points = length(batch_results);
method_tags = {'ours', 'rgd', 'sca', 'lm', 'elm'};
method_labels = {'Ours (SSL)', 'RGD', 'SCA', '(3)', '(6)'};
num_methods = length(method_tags);

% 提取真实轨迹与理论主轴方向 (针对 Cone 轨迹，理论主轴 r 应与位置 p 方向一致)
all_p_true = zeros(3, num_points);
all_r_theory = zeros(3, num_points);
for i = 1:num_points
    p_true = batch_results(i).test_point.p_true;
    all_p_true(:, i) = p_true;
    % 理论方向：由坐标原点指向磁铁位置的单位向量，仅用于绘图参考
    all_r_theory(:, i) = p_true / norm(p_true);
end

% 提取位置误差和指向误差 (保持与 plot_path_comparison.m 一致的误差源)
pos_err_mat = zeros(num_points, num_methods);
rot_err_mat = zeros(num_points, num_methods);

for i = 1:num_points
    s = batch_results(i).summary;
    for j = 1:num_methods
        tag = method_tags{j};
        
        % 位置误差：ours, rgd, sca 使用相同的 p_ours
        if strcmp(tag, 'rgd') || strcmp(tag, 'sca')
            if isfield(s, 'ours'), pos_err_mat(i, j) = s.ours.pos_mean; end
        elseif isfield(s, tag)
            pos_err_mat(i, j) = s.(tag).pos_mean;
        else
            pos_err_mat(i, j) = NaN;
        end
        
        % 指向误差 (r-axis error) - 采用与 plot_path_comparison.m 一致的提取逻辑
        if strcmp(tag, 'ours')
             if isfield(s.ours, 'direct_r_mean_SSL')
                rot_err_mat(i, j) = s.ours.direct_r_mean_SSL;
             else
                rot_err_mat(i, j) = s.ours.r_mean;
             end
          elseif strcmp(tag, 'rgd')
                 if isfield(s.ours, 'direct_r_mean_RGD')
                     rot_err_mat(i, j) = s.ours.direct_r_mean_RGD;
                 elseif isfield(s.ours, 'direct_r_mean_RO')
                     rot_err_mat(i, j) = s.ours.direct_r_mean_RO;
                 end
        elseif strcmp(tag, 'sca')
             if isfield(s.ours, 'direct_r_mean_MM')
                rot_err_mat(i, j) = s.ours.direct_r_mean_MM;
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
% 将截断后的误差归一化到 [0, 1] 范围
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
fig_handle = figure('Color', 'white', 'Name', 'Principal Axis Estimation Comparison', ...
                   'Units', 'centimeters', 'Position', [2, 2, 38, 12.03]);
tiledlayout(5, num_methods, 'Padding', 'compact', 'TileSpacing', 'compact');

fontsize = 16;

for j = 1:num_methods
    % --- 绘制 3D 轨迹图 (占用前 3 行) ---
    nexttile(j, [3, 1]);
    hold on; grid on; box on;
    
    % 绘制串联线 (背景线)
    plot3(all_p_true(1, :), all_p_true(2, :), all_p_true(3, :), '-', ...
          'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    
    % 绘制轨迹点
    arrow_len = 0.02; % 箭头长度 2cm

    for i = 1:num_points
        c_idx = max(1, min(n_colors, round(error_combined(i, j) * (n_colors - 1)) + 1));
        arrow_color = custom_colormap(c_idx, :);
        
        % 绘制散点 (黑色，缩小)
        scatter3(all_p_true(1, i), all_p_true(2, i), all_p_true(3, i), ...
                 12, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
        
        % 提取真值 r_true (如果有) 或 径向理论向量 (all_r_theory)
        r_plot = all_r_theory(:, i); 
        if isfield(batch_results(i).test_point, 'r_true')
            r_plot = batch_results(i).test_point.r_true;
        end
        
        % 修正：如果箭头反向了，在此处取反 (根据用户反馈全部反向)
        r_plot = -r_plot;
        
        % 绘制带颜色的箭头 (方向固定为真值，颜色反映估算误差)
        if norm(r_plot) > 1e-6
            r_plot = r_plot / norm(r_plot) * arrow_len;
            quiver3(all_p_true(1, i), all_p_true(2, i), all_p_true(3, i), ...
                    r_plot(1), r_plot(2), r_plot(3), ...
                    0, 'Color', arrow_color, 'LineWidth', 1.5, 'MaxHeadSize', 0.8, 'AutoScale', 'off'); 
        end
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
    
        % xlim([0, 0.1]);
        axis auto
    
        % Rotation Error Boxchart (Row 5)
        nexttile(4 * num_methods + j);
        er_data = rot_err_mat(:, j);
        boxchart(categorical(repmat("er", length(er_data), 1)), er_data, ...
                 'BoxFaceColor', [0.4 0.4 0.4], 'Orientation', 'horizontal');
        grid on; box on;
        if j == 1, yticklabels('$e_{\mbox{\boldmath{$r$}}}$ [rad]'); else, yticklabels(''); end
        set(gca, 'FontSize', fontsize, 'TickLabelInterpreter', 'latex', 'TickDir', 'out');
        
        % xlim([0, 0.1]);
        axis auto
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
    filename = sprintf('results/Path_Cone_N%d.png', num_sensors);
else
    filename = 'results/Path_Cone.png';
end

saveas(fig_handle, filename);
fprintf('[plot_path_Cone] 图像已保存至: %s\n', filename);

end
