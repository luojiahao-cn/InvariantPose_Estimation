function [p_est, R_est, stats] = estimate_pose_fischer( ...
    b_total, d_list, m_pos, m_hat, m_norm, ...
    theta_init, p_init, options, lb_p, ub_p, params)

    %#ok<*INUSD>  % 先不使用 theta_init 和 params，保留接口
    num_sensors = size(b_total, 2);

    if num_sensors > 6
        d_list = d_list(:, [1, 3, 7]);
        b_total = b_total(:, [1, 3, 7]);
    elseif (num_sensors <= 6) && (num_sensors > 3)
        d_list = d_list(:, [1, 2, 4]);
        b_total = b_total(:, [1, 2, 4]);
    else
        d_list = d_list(:, 1:3);
        b_total = b_total(:, 1:3);
    end
    num_sensors = 3;

    %% 利用三颗传感器估计局部梯度张量 X_opt（6×5 线性系统）
    % s1 设为参考点
    d = d_list(:, 2) - d_list(:, 1);   % s2 - s1
    e = d_list(:, 3) - d_list(:, 1);   % s3 - s1

    Db12 = b_total(:, 2) - b_total(:, 1);  % b2 - b1
    Db13 = b_total(:, 3) - b_total(:, 1);  % b3 - b1

    dx = d(1); dy = d(2); dz = d(3);
    ex = e(1); ey = e(2); ez = e(3);

    % 未知向量 x = [Gxx, Gxy, Gxz, Gyy, Gyz]'
    % Gzz = -Gxx - Gyy 由无散条件给出
    C = [ dx    dy    dz     0     0  ;
          0     dx    0      dy    dz ;
         -dz    0     dx    -dz    dy ;
          ex    ey    ez     0     0  ;
          0     ex    0      ey    ez ;
         -ez    0     ex    -ez    ey ];

    h = [ Db12(1);
          Db12(2);
          Db12(3);
          Db13(1);
          Db13(2);
          Db13(3) ];

    rankC = rank(C);
    if rankC < 5
        warning('estimate_pose_fischer:rankDeficient', ...
            'Gradient system matrix C is rank deficient. Sensor geometry may be ill conditioned.');
    end

    x = C \ h;   % 最小二乘解（满秩时等价于 pinv(C)*h）

    Gxx = x(1);
    Gxy = x(2);
    Gxz = x(3);
    Gyy = x(4);
    Gyz = x(5);
    Gzz = -Gxx - Gyy;

    X_opt = [ Gxx Gxy Gxz ;
              Gxy Gyy Gyz ;
              Gxz Gyz Gzz ];

    %% 构建磁场差矩阵和位移矩阵
    % X_opt_ours = lc_grad_tensor_estimator(b_total, d_list);
    %% 传感器系特征分解，得到 R_sv
    [V_s, D_s] = eig(X_opt);
    [lambda_s, idx_s] = sort(diag(D_s), 'ascend');
    V_s = V_s(:, idx_s);

    % 保证是右手正交基
    if det(V_s) < 0
        V_s(:,1) = -V_s(:,1);
    end
    Rsv = V_s;             % eigenvector matrix in {s}
    lambda_s = lambda_s(:);

    %% 构造测量不变量（模长和特征值）
    b_s = mean(b_total, 2);           % 用所有传感器平均场近似局部 b
    norm_b_meas = norm(b_s);
    norm_G_meas = norm(X_opt, 'fro');

    invariants_meas = [norm_b_meas; norm_G_meas; lambda_s];

    %% 解析 options
    epsilon       = 0.1; % 0.3
    grid_step     = 0.02; % 0.02
    num_iter      = 4;
    shrink_factor = 1.4;

    %% 粗到细网格搜索估计位置
    p_est = p_init;
    cost_prev = inf;

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

            % 模型场和梯度（世界坐标系）
            [b_w, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);

            % 模型不变量
            norm_b_mod = norm(b_w);
            norm_G_mod = norm(A_p, 'fro');
            lambda_mod = sort(eig(A_p), 'ascend');
            lambda_mod = lambda_mod(:);

            invariants_mod = [norm_b_mod; norm_G_mod; lambda_mod];

            % 代价函数
            cost = norm(invariants_mod - invariants_meas);

            if cost < cost_prev
                cost_prev = cost;
                p_est     = p;
            end
        end

        cost_history(it) = cost_prev;
        p_history(:,it)  = p_est;

        % 缩小搜索范围和步长
        epsilon   = epsilon   / shrink_factor;
        grid_step = grid_step / shrink_factor;
    end

    %% 姿态估计：在 p_est 处再算一次梯度，特征向量配准
    [b_p, A_p_final] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);

    % % % 对d_list'进行QR分解
    % [Q, ~] = qr(d_list');
    % r = rank(d_list); % 构型判据
    % Q_bar = Q(:, r+1:end);
    % b_bar = b_total * Q_bar; % 计算bBar

    % B_matrix = b_p * ones(1, num_sensors);
    % B_bar = B_matrix * Q_bar;

    % [R_est, ~] = estimateR(b_bar, B_bar, A_p_final, X_opt);

    [V_w, D_w] = eig(A_p_final);
    [lambda_w, idx_w] = sort(diag(D_w), 'ascend');
    V_w = V_w(:, idx_w);

    if det(V_w) < 0
        V_w(:,1) = -V_w(:,1);
    end

    Rwv = V_w;        % eigenvector matrix in {w}
    lambda_w = lambda_w(:);

    % 传感器系到世界系的旋转
    R_est = Rwv * Rsv';

    %% 输出一些调试信息
    stats = struct();
    stats.X_opt        = X_opt;
    stats.cost_history = cost_history;
    stats.p_history    = p_history;
end
