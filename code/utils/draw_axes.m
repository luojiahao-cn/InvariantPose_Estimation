function draw_axes(T, label, scale)
    % draw_axes 画出给定齐次变换矩阵的坐标轴，并可以调整箭头大小
    % 输入:
    % T - 4x4 齐次变换矩阵
    % scale - 箭头长度的缩放因子

    % 检查输入矩阵尺寸
    if size(T, 1) ~= 4 || size(T, 2) ~= 4
        error('输入矩阵必须是4x4的');
    end

    % 检查label参数
    if nargin < 2
        label = ''; % 如果没有指定，不显示标签
    end

    % 检查scale参数
    if nargin < 3
        scale = 0.2; % 如果没有指定，默认为1
    end

    % 坐标轴的原点
    origin = T(1:3, 4);

    % 坐标轴方向
    x_axis = T(1:3, 1);
    y_axis = T(1:3, 2);
    z_axis = T(1:3, 3);

    % 规范化向量并应用缩放因子
    x_axis = x_axis / norm(x_axis) * scale;
    y_axis = y_axis / norm(y_axis) * scale;
    z_axis = z_axis / norm(z_axis) * scale;

    % 绘制x轴
    quiver3(origin(1), origin(2), origin(3), x_axis(1), x_axis(2), x_axis(3), 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5);

    % 如果指定了标签，则添加
    if ~isempty(label)
        text(origin(1), origin(2), origin(3), ['  \{', label, '\}'], 'FontSize', 24, 'Interpreter', 'latex');
    end

    % 绘制y轴
    quiver3(origin(1), origin(2), origin(3), y_axis(1), y_axis(2), y_axis(3), 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % 绘制z轴
    quiver3(origin(1), origin(2), origin(3), z_axis(1), z_axis(2), z_axis(3), 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% 
%     % 设置图形属性
%     axis equal;
%     grid on;
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     view(3); % 三维视图
end
