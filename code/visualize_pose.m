function visualize_pose(m_pos, m_hat, m_norm, p_true, R_true, p_est, R_est, d_list)
    figure;
    hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('磁铁与传感器阵列位姿可视化');

    % 绘制磁铁位置
    scatter3(m_pos(1,:), m_pos(2,:), m_pos(3,:), 80, 'r', 'filled');
    for i = 1:size(m_pos,2)
        % 绘制磁化方向
        quiver3(m_pos(1,i), m_pos(2,i), m_pos(3,i), ...
            m_hat(1,i), m_hat(2,i), m_hat(3,i), ...
            0.03 * m_norm(i)/max(m_norm), 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    end

    % 绘制真实传感器阵列参考点
    scatter3(p_true(1), p_true(2), p_true(3), 100, 'b', 'filled');
    % 绘制估计传感器阵列参考点
    scatter3(p_est(1), p_est(2), p_est(3), 100, 'm', 'filled');

    % 计算真实传感器阵列矩形顶点
    rect_local = [d_list(:,2), d_list(:,3), d_list(:,5), d_list(:,4)];
    rect_true = p_true + R_true * rect_local;
    fill3(rect_true(1,:), rect_true(2,:), rect_true(3,:), [0.7 0.9 0.7], 'FaceAlpha',0.5, 'EdgeColor','g', 'LineWidth',2);

    % 计算估计传感器阵列矩形顶点
    rect_est = p_est + R_est * rect_local;
    fill3(rect_est(1,:), rect_est(2,:), rect_est(3,:), [0.9 0.7 0.9], 'FaceAlpha',0.5, 'EdgeColor','m', 'LineWidth',2);

    % 绘制真实阵列方向向量
    dir_vec_true = R_true(:,1) * 0.01;
    quiver3(p_true(1), p_true(2), p_true(3), dir_vec_true(1), dir_vec_true(2), dir_vec_true(3), ...
        1, 'b', 'LineWidth', 2, 'MaxHeadSize', 2);

    % 绘制估计阵列方向向量
    dir_vec_est = R_est(:,1) * 0.01;
    quiver3(p_est(1), p_est(2), p_est(3), dir_vec_est(1), dir_vec_est(2), dir_vec_est(3), ...
        1, 'm', 'LineWidth', 2, 'MaxHeadSize', 2);

    legend({'磁铁','磁化方向','真实参考点','估计参考点','真实阵列','估计阵列','真实方向','估计方向'}, 'Location', 'best');
    view([45 30]);
    hold off;
end