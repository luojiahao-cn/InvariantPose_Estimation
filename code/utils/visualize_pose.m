function visualize_pose(m_pos, m_hat, m_norm, p_true, R_true, results, d_list)
    figure;
    hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('磁铁与传感器阵列位姿可视化');

    % 绘制磁铁位置
    h1 = scatter3(m_pos(1,:), m_pos(2,:), m_pos(3,:), 80, 'r', 'filled');
    h2 = gobjects(1, size(m_pos,2));
    for i = 1:size(m_pos,2)
        % 绘制磁化方向
        h2(i) = quiver3(m_pos(1,i), m_pos(2,i), m_pos(3,i), ...
            m_hat(1,i), m_hat(2,i), m_hat(3,i), ...
            0.03 * m_norm(i)/max(m_norm), 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    end

    % 绘制真实传感器阵列参考点
    h3 = scatter3(p_true(1), p_true(2), p_true(3), 100, 'b', 'filled');

    % 计算真实传感器阵列矩形顶点
    rect_local = [d_list(:,2), d_list(:,3), d_list(:,5), d_list(:,4)];
    rect_true = p_true + R_true * rect_local;
    h4 = fill3(rect_true(1,:), rect_true(2,:), rect_true(3,:), [0.7 0.9 0.7], 'FaceAlpha',0.5, 'EdgeColor','g', 'LineWidth',2);

    % 绘制真实阵列方向向量
    dir_vec_true = R_true(:,1) * 0.01;
    h5 = quiver3(p_true(1), p_true(2), p_true(3), dir_vec_true(1), dir_vec_true(2), dir_vec_true(3), ...
        1, 'b', 'LineWidth', 2, 'MaxHeadSize', 2);

    % 支持的算法及其显示属性
    algs = {
        'lm',    'p_lm',    'R_lm',    'm',           [0.9 0.7 0.9], 'LM估计';
        'elm',   'p_elm',   'R_elm',   [0.2 0.6 1],   [0.2 0.6 1],   'ELM估计';
        'Rlm',   'p_Rlm',   'R_Rlm',   [0.2 0.8 0.2], [0.2 0.8 0.2], 'RLm估计';
        % 新算法可在此添加
        'ours',  'p_ours',  'R_ours',  [1 0.6 0.2],   [1 0.6 0.2],   'OURS估计';
    };

    % 过滤异常点（距离真实参考点超过阈值）
    threshold = 0.5;
    valid_idx = true(1, length(results));
    for k = 1:length(results)
        dists = [];
        for a = 1:size(algs,1)
            p_field = algs{a,2};
            if isfield(results(k), p_field)
                dists(end+1) = norm(results(k).(p_field) - p_true);
            end
        end
        if ~isempty(dists) && max(dists) > threshold
            valid_idx(k) = false;
        end
    end
    results = results(valid_idx);

    % 批量绘制所有实验结果
    h_alg = cell(1, size(algs,1));
    for k = 1:length(results)
        for a = 1:size(algs,1)
            p_field = algs{a,2};
            R_field = algs{a,3};
            scatter_color = algs{a,4};
            rect_color = algs{a,5};
            if isfield(results(k), p_field) && isfield(results(k), R_field)
                p_val = results(k).(p_field);
                h_alg{a}(end+1) = scatter3(p_val(1), p_val(2), p_val(3), 60, scatter_color, 'filled');
                rect = p_val + results(k).(R_field) * rect_local;
                fill3(rect(1,:), rect(2,:), rect(3,:), rect_color, 'FaceAlpha',0.3, 'EdgeColor',rect_color, 'LineWidth',1);
                dir_vec = results(k).(R_field)(:,1) * 0.01;
                quiver3(p_val(1), p_val(2), p_val(3), dir_vec(1), dir_vec(2), dir_vec(3), ...
                    'LineWidth', 1, 'Color', scatter_color, 'MaxHeadSize', 1);
            end
        end
    end

    % 图例句柄和标签
    legend_handles = [h1, h2(1), h3, h4, h5];
    legend_labels = {'磁铁位置', '磁化方向', '真实参考点', '真实阵列', '真实方向'};
    for a = 1:size(algs,1)
        if ~isempty(h_alg{a})
            legend_handles(end+1) = h_alg{a}(1);
            legend_labels{end+1} = algs{a,6};
        end
    end
    legend(legend_handles, legend_labels, 'Location', 'bestoutside');

    % 设置视角以真实传感器阵列参考点为中心
    camtarget(p_true'); % 视线焦点
    campos(p_true' + [0.1 0.1 0.1]); % 视点位置，可根据实际场景调整距离
    camva(6); % 视场角，可根据实际场景调整
    axis vis3d; % 保持三维比例
    view(3);    % 强制3D视图

    hold off;
end