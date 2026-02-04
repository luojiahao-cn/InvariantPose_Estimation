function plot_principal_axis_results(batch_results, params)
% PLOT_PRINCIPAL_AXIS_RESULTS 专门用于主轴估计性能对比的绘图函数
% 输入：
%   batch_results - 由 run_batch_experiments 返回的批量结果结构体数组
%   params        - 实验参数结构体
close all;
if nargin < 2, params = struct(); end

%% 1. 数据提取与预处理
num_points = length(batch_results);
num_methods = 5;
method_labels = {'SCA', 'RO', 'SDP', 'SDP RED', 'SPEC'};

all_coords = [batch_results.test_point];
all_p_true = [all_coords.p_true]; 
all_theta_true = [all_coords.theta_true];

% 虽然位置几乎不变，但我们仍提取 p 的平均值作为绘图中心
p_avg = mean(all_p_true, 2);

% 计算主轴 (principal axis) 指向误差对比
rot_error_matrix = zeros(num_points, num_methods);
for i = 1:num_points
    if isfield(batch_results(i).summary, 'ours')
        s = batch_results(i).summary.ours;
        rot_error_matrix(i, 1) = s.direct_r_mean_MM;           % SCA
        rot_error_matrix(i, 2) = s.direct_r_mean_RO;           % RO
        rot_error_matrix(i, 3) = s.direct_r_mean_SDP;          % SDP
        rot_error_matrix(i, 4) = s.direct_r_mean_SDP_reduced;  % Ours (SDP RED)
        rot_error_matrix(i, 5) = s.direct_r_mean_spec;         % SPEC
    end
end

% 归一化误差到[0,1]用于颜色映射 (越绿误差越小，越红越大)
delta_err = 0.05; 
err_clamped = rot_error_matrix;
err_clamped(err_clamped > delta_err) = delta_err;
err_max = max(err_clamped(:));
err_min = min(err_clamped(:));
if err_max > err_min
    error_norm = (err_clamped - err_min) / (err_max - err_min);
else
    error_norm = ones(size(err_clamped)) * 0.5;
end

%% 2. 绘图执行
render_one_figure('Principal Axis Estimation', batch_results, p_avg, all_theta_true, error_norm, ...
                 rot_error_matrix, params);

end

%% ========================== 核心渲染子函数 ========================== %%
function render_one_figure(title_tag, batch_results, p_avg, all_theta_true, error_norm, ...
                          rot_err_mat, params)
    
    num_points = length(batch_results);
    num_methods = 5;
    method_labels = {'SCA', 'RO', 'SDP', 'SDP RED', 'SPEC'};
    fontsize = 16;
    n_colors = 256;
    % slanCM('gor') 是从绿到红的渐变，若无此函数则回退到 jet 翻转
    try
        custom_colormap = slanCM('gor', n_colors);
    catch
        custom_colormap = flipud(jet(n_colors)); 
    end

    % 获取参数中的磁矩默认指向
    m_hat = [1; 0; 0];
    if isfield(params, 'magnet') && isfield(params.magnet, 'm_hat')
        temp_m = params.magnet.m_hat;
        m_hat = temp_m(:, 1);
    end

    n_sensors = 0;
    if isfield(params, 'sensor') && isfield(params.sensor, 'd_list')
        n_sensors = size(params.sensor.d_list, 2);
    end

    fig_name = sprintf('Results: %s', title_tag);
    % 调整画布比例以容纳 5 个图
    fig_handle = figure('Color', 'white', 'units', 'centimeters', 'Position', [2, 2, 28, 8], 'Name', fig_name);
    tl = tiledlayout(3, num_methods, 'TileSpacing', 'loose', 'Padding', 'tight');
    
    % --- 绘制 3D 矢量分布图 ---
    for j = 1:num_methods
        nexttile(j, [2, 1]);
        hold on; grid on; box on;
        
        % 1. 绘制极简空间参考
        [gx, gy, gz] = meshgrid(linspace(-0.05, 0.05, 5));
        scatter3(gx(:)+p_avg(1), gy(:)+p_avg(2), gz(:)+p_avg(3), 1, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.05);

        % 2. 绘制带颜色的方向矢量 (抽样显示)
        step = max(1, floor(num_points / 40)); 
        subset_idx = 1:step:num_points;
        
        for idx = subset_idx
             % 轴角转旋转矩阵
             th = all_theta_true(:, idx);
             angle = norm(th);
             if angle < 1e-10
                 R = eye(3);
             else
                 paxis = th / angle;
                 K = [0 -paxis(3) paxis(2); paxis(3) 0 -paxis(1); -paxis(2) paxis(1) 0];
                 R = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;
             end
             r_vec = R * m_hat;
             
             % 计算该点对应的颜色索引
             c_idx = max(1, min(n_colors, round(error_norm(idx, j) * (n_colors - 1)) + 1));
             c = custom_colormap(c_idx, :);
             
             % 绘制箭头
             quiver3(p_avg(1), p_avg(2), p_avg(3), ...
                     r_vec(1)*0.04, r_vec(2)*0.04, r_vec(3)*0.04, ...
                     'Color', [c, 0.8], 'LineWidth', 1.5, 'AutoScale', 'off', 'MaxHeadSize', 0.5);
        end
        
        % 视图设置
        axis equal; axis square; view(122, 13);
        xlabel('$x$ [m]', 'Interpreter', 'latex');
        ylabel('$y$ [m]', 'Interpreter', 'latex');
        zlabel('$z$ [m]', 'Interpreter', 'latex');
        
        % 限制范围在中心点附近
        range = 0.05;
        xlim([p_avg(1)-range, p_avg(1)+range]);
        ylim([p_avg(2)-range, p_avg(2)+range]);
        zlim([p_avg(3)-range, p_avg(3)+range]);
        
        set(gca, 'FontSize', fontsize-2, 'TickLabelInterpreter', 'latex', 'TickDir', 'out', 'LineWidth', 1);
        title(method_labels{j}, 'FontSize', fontsize, 'Interpreter', 'latex');
    end
    
    % --- 绘制误差箱线图 (仅保留 e_r) ---
    for j = 1:num_methods
        subtl = tiledlayout(tl, 1, 1);
        subtl.Layout.Tile = 10 + j;
        
        nexttile(subtl);
        r_errs = rot_err_mat(:, j);
        gc_r = categorical(repmat("Rotation", length(r_errs), 1));
        boxchart(gc_r, r_errs, 'BoxFaceColor', [0.2 0.2 0.2], 'Orientation', 'horizontal');
        
        if j == 1, yticklabels('$e_{\mathbf{r}}$'); else, yticklabels(''); end
        set(gca, 'FontSize', fontsize-2, 'TickDir', 'out', 'TickLabelInterpreter', 'latex', 'LineWidth', 1); 
        grid on; box on;
    end

    % --- 导出图像 ---
    if ~exist('results', 'dir'), mkdir('results'); end
    filename = sprintf('results/PrincipalAxis_Comparison_Ns%d.png', n_sensors);
    
    try
        export_fig(fig_handle, filename, '-nocrop', '-r600', '-transparent');
        fprintf('[plot_principal_axis_results] 图像已导出: %s\n', filename);
    catch
        saveas(fig_handle, filename);
    end
end
