function [p_est, R_est, stats] = estimate_pose_ours(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, params)
% PROPOSED_METHOD_POSE_ESTIMATION 使用所提方法估计传感器姿态（位置和方向）
% 该方法使用公式(11)
% 输入参数：
%   b_total  - 3×N传感器测量矩阵（局部坐标系），单位：T
%   d_list   - 3×N位移矩阵，传感器在参考坐标系中的偏移，单位：m
%   m_pos    - 3×K磁铁位置矩阵，单位：m
%   m_hat    - 3×K磁化方向单位向量（归一化）
%   m_norm   - 1×K磁矩幅值向量，单位：A·m²
%   p_init   - 3×1初始位置估计
%
% 输出参数：
%   p_est_22 - 第二阶段优化后的位置估计
%   R_est    - 估计的旋转矩阵
%   stats    - 包含中间结果和统计信息的结构体
R_init = MatrixExp3(VecToso3(theta_init));
num_sensors = size(b_total, 2);
%% 构建磁场差矩阵和位移矩阵
[X_opt, ~, D_delta, B_delta] = lc_grad_tensor_estimator(b_total, d_list);

%% Stage #1: Estimate for position \hat{p}
[Q, ~, ~] = qr(d_list'); % Note: must return 3 outputs even we only need Q
r = rank(d_list); % 构型判据
Q_bar = Q(:, r+1:end);
b_bar = b_total * Q_bar; % 计算bBar

fun22 = @(p) obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt, params.W);
% 粗搜索+精搜索
% p_est = grid_search(p_init, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt, lb_p, ub_p, params.W);
p_est = p_init;
% options.FunctionTolerance = 1e-8;
[p_est, ~, ~, ~, output] = lsqnonlin(fun22, p_est, lb_p, ub_p, options);
% output.message
%% Stage #2: Estimate for rotation \hat{R}
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B_matrix = b_p * ones(1, num_sensors);
B_bar = B_matrix * Q_bar;
[R_init_est1, R_init_est2] = estimateR(b_bar, B_bar, A_p, X_opt, D_delta, B_delta);

beta = 1e2;
% 应该要以权重项为基准
R_SSL = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init, params.mu, beta, params.R_true); % using R_init

R_est = R_SSL.R;
stats.X_opt = X_opt;         % 估计的梯度矩阵
stats.b_bar = b_bar;         % 测量磁场投影矩阵
stats.B_bar = B_bar;         % 预测磁场投影矩阵
stats.A_p = A_p;             % 位置估计处的梯度张量
stats.R_est_init1 = R_init_est1;
stats.R_est_init2 = R_init_est2;
stats.R_iter_history = R_SSL.R_iter_history; % 每次迭代的R
stats.delta_history = R_SSL.delta_history;   % 每次迭代的delta
end

%% ----------------------------Functions-------------------------------  %%
%% Estimate p
function res = obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X, W)
    [b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);
    B_bar = b_p * ones(1, num_sensors) * Q_bar;
    % term1 = norm(B_bar, 'fro') - norm(b_bar, 'fro');
    term1 = sort(eig(B_bar * B_bar'), 'descend')  - sort(eig(b_bar * b_bar'), 'descend'); % 矩阵特征值匹配
    term2 = sort(eig(A_p), 'descend') - sort(eig(X), 'descend'); % 特征值匹配
    res = [term1; term2]; % 只需要匹配最大的两个特征值，因为迹0
    % res = [term1; term2; term5];
    % res = W * res;
end

%% 网格粗搜索
function p_est = grid_search(p_init, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt, lb_p, ub_p, W)
    % 解析 options
    epsilon       = 0.1;
    grid_step     = 0.002;
    num_iter      = 3;
    shrink_factor = 1.4;

    % 粗到细网格搜索估计位置
    cost_prev = inf;

    p_est = p_init;
    cost_history = nan(num_iter, 1);
    p_history    = nan(3, num_iter);

    for it = 1:num_iter

        % 当前搜索立方体范围
        x_min = max(lb_p(1), p_est(1) - epsilon/2);
        x_max = min(ub_p(1), p_est(1) + epsilon/2);
        y_min = max(lb_p(2), p_est(2) - epsilon/2);
        y_max = min(ub_p(2), p_est(2) + epsilon/2);
        z_min = max(lb_p(3), p_est(3) - epsilon/2);
        z_max = min(ub_p(3), p_est(3) + epsilon/2);

        x_grid = x_min:grid_step:x_max;
        y_grid = y_min:grid_step:y_max;
        z_grid = z_min:grid_step:z_max;

        [Xg, Yg, Zg] = ndgrid(x_grid, y_grid, z_grid);
        P_grid = [Xg(:), Yg(:), Zg(:)];

        for k = 1:size(P_grid, 1)
            p = P_grid(k, :)';

            cost = norm(obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt, W));
            if cost < cost_prev
                cost_prev = cost;
                p_est     = p;
            end
        end

        % shading interp;
        cost_history(it) = cost_prev;
        p_history(:,it)  = p_est;

        % 缩小搜索范围和步长
        epsilon   = epsilon   / shrink_factor;
        grid_step = grid_step / shrink_factor;
    end

end