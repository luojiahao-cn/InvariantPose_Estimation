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
X_opt = lc_grad_tensor_estimator(b_total, d_list);

%% 阶段检查 对应公式2
% [~, A_p_true] = calcFieldAndGradient(params.p_true, m_pos, m_hat, m_norm);
% R_true = params.R_true;
% X_true = R_true'*A_p_true*R_true;
% norm(R_true'*b_p*ones(1,num_sensors)+R_true'*A_p*R_true*d_list - b_total, 'fro')
% norm(R_true'*A_p*R_true*D_matrix - B_matrix, 'fro')

%% Stage #1: Estimate for position \hat{p}
% 对d_list'进行QR分解
[Q, ~] = qr(d_list');
r = rank(d_list); % 构型判据
Q_bar = Q(:, r+1:end);
b_bar = b_total * Q_bar; % 计算bBar
% 优化第一阶段位置
fun22 = @(p) obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt);
p_est = lsqnonlin(fun22, p_init, lb_p, ub_p, options);
%% Stage #2: Estimate for rotation \hat{R}
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B_matrix = b_p * ones(1, num_sensors);
B_bar = B_matrix * Q_bar;
[R_init_est1, R_init_est2] = estimateR(b_bar, B_bar, A_p, X_opt);

mu = 1;
beta = 1e2;
R_PPI = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init_est1, mu, beta, params.R_true); % using R_init
stats.X_opt = X_opt;         % 估计的梯度矩阵
stats.R_iter_history = R_PPI.R_iter_history; % 每次迭代的R
stats.delta_history = R_PPI.delta_history;   % 每次迭代的delta
R_est.PPI = R_PPI.R; % 返回最终R
R_est.R_init_est1 = R_init_est1;
R_est.R_init_est2 = R_init_est2;
end

%% ----------------------------Functions-------------------------------  %%
%% Estimate p
function res = obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X)
    [b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);
    term1 = norm(b_p * ones(1, num_sensors) * Q_bar, 'fro') - norm(b_bar, 'fro');
    term2 = trace(A_p*A_p) - trace(X*X);
    term3 = det(A_p) - det(X);
    res = [term1; term2; term3];
end