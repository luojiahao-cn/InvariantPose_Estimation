clc,clear,close all
addpath('./utils')
addpath('./Functions')
addpath('./tools')

%% ========== 参数配置 ==========
params = get_experiment_params();

% 固定源信息
pC = params.magnet.m_pos;   % 3xM
M  = size(pC,2);

% 设置随机种子
rng(params.experiment.random_seed);

% 传感器构型（去中心化）
d_list = [
    [1e-3; 0; 0], ...
    [0; 0; 0], ...
    [0; 0; 1e-3]
];
d_list = d_list - mean(d_list, 2);
params.sensor.d_list = d_list;

%% ======== 生成测试点位 =========
num_test_points = 300;
test_points = generate_test_points(params, num_test_points, 'random');

%% ======== 预计算 Q_bar（与构型相关，和测试点无关）========
% Q_bar 需要满足 d_list' * Q_bar = 0
Q_bar = null(d_list'); % N x (N-r)  注意N=num_sensors，此处需要先知道N
% 但 N 只有在 generate_magnetic_data 后才知道
% 所以这里先留空，进入循环后第一次测量时再算也行
Q_bar_ready = false;

%% ======== 执行实验 ==========
err_norm.R = zeros(1, num_test_points);
err_norm.p = zeros(1, num_test_points);

% 保存基线驱动（bg）工作点：三个磁铁都有任务
bg_m_norm = params.magnet.m_norm;   % 1xM 或 Mx1，按你结构

m_norm_variation = 5; % 小增量幅值

for test_point_idx = 1:num_test_points

    %% --- 设置当前真值位姿 ---
    params.ground_truth.p_true = test_points(:, test_point_idx).p_true;
    params.ground_truth.theta_true = test_points(:, test_point_idx).theta_true;

    p_true = params.ground_truth.p_true;
    R_true = MatrixExp3(VecToso3(params.ground_truth.theta_true));

    %% --- 基线测量：同一位姿下，所有源都在 bg_m_norm ---
    params.magnet.m_norm = bg_m_norm;
    [b_bg, ~, ~] = generate_magnetic_data(params);   % 3xN
    num_sensors = size(b_bg, 2);

    % 只需算一次 Q_bar 与 g（与 num_sensors 相关，但对所有 test point 相同）
    if ~Q_bar_ready
        Q_bar = null(d_list');  % 这里 d_list' 是 N x 3，其中 N=num_sensors
        if size(Q_bar,1) ~= num_sensors
            error('Q_bar row size mismatch: check d_list dimension equals num_sensors.');
        end
        g = Q_bar' * ones(num_sensors, 1);
        if norm(g) < 1e-10
            error('g = Q_bar''*1 nearly zero; b_mag ill-conditioned for this array.');
        end
        Q_bar_ready = true;
    end

    %% --- 轮流小增益差分，得到每个源的 rho_i ---
    rho = zeros(3, M);

    for i = 1:M
        % 仅第 i 个源加增量
        params.magnet.m_norm = bg_m_norm;
        params.magnet.m_norm(i) = params.magnet.m_norm(i) + m_norm_variation;

        [b_plus, ~, ~] = generate_magnetic_data(params);   % 3xN
        dB = b_plus - b_bg;                                % 差分 -> 仅 Δu_i 响应（理想线性下）

        % ======== 计算差分中心场估计 Δb_i,s ========
        % dB*Q_bar ≈ (Δb_s) * (1^T Q_bar) -> 最小二乘恢复 Δb_s
        b_mag = dB * Q_bar * g / (norm(g)^2);              % 3x1

        % ======== 计算差分梯度张量 ΔX_i ========
        X_opt = lc_grad_tensor_estimator(dB, params.sensor.d_list); % 3x3

        % 可选：条件数保护，避免数值爆炸
        if cond(X_opt) > 1e10
            warning('test %d, source %d: X_opt ill-conditioned (cond=%.2e)', ...
                test_point_idx, i, cond(X_opt));
        end

        % ======== 欧拉齐次关系：rho_i = -3 * (ΔX_i)^{-1} * (Δb_i) ========
        rho(:,i) = -3 * (X_opt \ b_mag);
    end

    %% ======== 刚体配准求 (R, p) ========
    % 模型：rho_i = R^T (p - pC_i)
    % 令 s_i = -rho_i，则 pC_i = p + R*s_i
    S = -rho;   % 3xM  sensor points
    W = pC;     % 3xM  world points

    % Kabsch（更稳）
    s_bar = mean(S, 2);
    w_bar = mean(W, 2);
    S0 = S - s_bar;
    W0 = W - w_bar;

    H = S0 * W0';
    [U,~,V] = svd(H);
    D = diag([1, 1, det(V*U')]);
    R_est = V * D * U';
    p_est = w_bar - R_est * s_bar;

    %% ======== 误差统计 ========
    % 旋转误差：用 geodesic 距离更稳：acos((trace(R_true^T R_est)-1)/2)
    R_rel = R_true' * R_est;
    angle_err = acos(max(-1,min(1,(trace(R_rel)-1)/2)));
    err_norm.R(test_point_idx) = angle_err;

    err_norm.p(test_point_idx) = norm(p_true - p_est);

end

figure;
plot(err_norm.p)
xlabel('测试点索引');
ylabel('位置误差 norm(p_{true}-p_{est}) (m)');

figure;
plot(err_norm.R)
xlabel('测试点索引');
ylabel('旋转误差 angle(R_{true},R_{est}) (rad)');
