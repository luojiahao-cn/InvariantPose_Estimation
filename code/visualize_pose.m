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

    % 过滤异常点（距离真实参考点超过阈值）
    threshold = 0.5; % 距离阈值（单位：米）
    valid_idx = true(1, length(results));
    for k = 1:length(results)
        dist_lm = inf;
        dist_elm = inf;
        dist_prop = inf;
        if isfield(results(k), 'p_lm')
            dist_lm = norm(results(k).p_lm - p_true);
        end
        if isfield(results(k), 'p_elm')
            dist_elm = norm(results(k).p_elm - p_true);
        end
        if isfield(results(k), 'p_prop')
            dist_prop = norm(results(k).p_prop - p_true);
        end

        if max([dist_lm, dist_elm, dist_prop]) > threshold
            valid_idx(k) = false;
        end
    end
    results = results(valid_idx);

    % 批量绘制所有实验结果
    h_lm = [];
    h_elm = [];
    h_prop = [];
    for k = 1:length(results)
        % LM
        if isfield(results(k), 'p_lm')
            h_lm(end+1) = scatter3(results(k).p_lm(1), results(k).p_lm(2), results(k).p_lm(3), 60, 'm', 'filled');
            rect_lm = results(k).p_lm + results(k).R_lm * rect_local;
            fill3(rect_lm(1,:), rect_lm(2,:), rect_lm(3,:), [0.9 0.7 0.9], 'FaceAlpha',0.3, 'EdgeColor','m', 'LineWidth',1);
            dir_vec_lm = results(k).R_lm(:,1) * 0.01;
            quiver3(results(k).p_lm(1), results(k).p_lm(2), results(k).p_lm(3), dir_vec_lm(1), dir_vec_lm(2), dir_vec_lm(3), ...
                'LineWidth', 1, 'Color', 'm', 'MaxHeadSize', 1);
        end
        % ELM
        if isfield(results(k), 'p_elm')
            h_elm(end+1) = scatter3(results(k).p_elm(1), results(k).p_elm(2), results(k).p_elm(3), 60, [0.2 0.6 1], 'filled');
            rect_elm = results(k).p_elm + results(k).R_elm * rect_local;
            fill3(rect_elm(1,:), rect_elm(2,:), rect_elm(3,:), [0.2 0.6 1], 'FaceAlpha',0.3, 'EdgeColor',[0.2 0.6 1], 'LineWidth',1);
            dir_vec_elm = results(k).R_elm(:,1) * 0.01;
            quiver3(results(k).p_elm(1), results(k).p_elm(2), results(k).p_elm(3), dir_vec_elm(1), dir_vec_elm(2), dir_vec_elm(3), ...
                'LineWidth', 1, 'Color', [0.2 0.6 1], 'MaxHeadSize', 1);
        end
        % OURS
        if isfield(results(k), 'p_prop')
            h_prop(end+1) = scatter3(results(k).p_prop(1), results(k).p_prop(2), results(k).p_prop(3), 60, [1 0.6 0.2], 'filled');
            rect_prop = results(k).p_prop + results(k).R_prop * rect_local;
            fill3(rect_prop(1,:), rect_prop(2,:), rect_prop(3,:), [1 0.6 0.2], 'FaceAlpha',0.3, 'EdgeColor',[1 0.6 0.2], 'LineWidth',1);
            dir_vec_prop = results(k).R_prop(:,1) * 0.01;
            quiver3(results(k).p_prop(1), results(k).p_prop(2), results(k).p_prop(3), dir_vec_prop(1), dir_vec_prop(2), dir_vec_prop(3), ...
                'LineWidth', 1, 'Color', [1 0.6 0.2], 'MaxHeadSize', 1);
        end
    end

    % 分别设置图例，确保颜色和内容对应
    legend([h1, h2(1), h3, h4, h5, h_lm(1), h_elm(1), h_prop(1)], {
        '磁铁位置',... 
        '磁化方向',... 
        '真实参考点',... 
        '真实阵列',... 
        '真实方向',... 
        'LM估计',...
        'ELM估计',... 
        'OURS估计'
    }, 'Location', 'bestoutside');

    % 设置视角以真实传感器阵列参考点为中心
    camtarget(p_true'); % 视线焦点
    campos(p_true' + [0.1 0.1 0.1]); % 视点位置，可根据实际场景调整距离
    camva(6); % 视场角，可根据实际场景调整
    axis vis3d; % 保持三维比例
    view(3);    % 强制3D视图

    hold off;
end