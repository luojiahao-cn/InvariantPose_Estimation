function results = run_beta_experiment(exp_idx, num_experiments, params, ...
    b_total, d_list, p_true)
% RUN_BETA_EXPERIMENT 执行beta参数实验
% 输入：
%   exp_idx        - 当前实验索引
%   num_experiments - 总实验次数
%   params         - 实验参数结构体
%   b_total        - 3×N传感器测量矩阵（局部坐标系）
%   d_list         - 3×N传感器位移矩阵
%   p_true         - 真实位置
%   R_true         - 真实旋转矩阵
%   lb_p           - 位置下界
%   ub_p           - 位置上界
% 输出：
%   result - 包含所有算法结果和误差的结构体

fprintf('\n===== 实验 %d/%d =====\n', exp_idx, num_experiments);

% 提取参数
m_pos = params.magnet.m_pos;
m_hat = params.magnet.m_hat;
m_norm = params.magnet.m_norm;
theta_true = params.ground_truth.theta_true;
p_uncertainty = params.uncertainty.p_uncertainty;
r_uncertainty = params.uncertainty.r_uncertainty;
mu = params.optimization.mu;

% ==== 初始值猜测 ====
init_error = -1 + 2 * rand(3,1); % [-1, 1]
p_init = p_true + p_uncertainty * init_error;
theta_init = theta_true + r_uncertainty * init_error;
% ==== 算法调用 ====
R_true = MatrixExp3(VecToso3(theta_true));
R_init = MatrixExp3(VecToso3(theta_init));
%% 构建磁场差矩阵和位移矩阵
num_sensors = size(b_total, 2);
X_opt = lc_grad_tensor_estimator(b_total, d_list); % 估计局部梯度张量

% 对d_list'进行QR分解
[Q, ~] = qr(d_list');
r = rank(d_list); % 构型判据
Q_bar = Q(:, r+1:end);
b_bar = b_total * Q_bar; % 计算bBar

[b_p, A_p] = calcFieldAndGradient(p_init, m_pos, m_hat, m_norm); % 基于初始位置估计计算磁场和梯度张量
B_matrix = b_p * ones(1, num_sensors);
B_bar = B_matrix * Q_bar;

[R_init_est1, R_init_est2] = estimateR(b_bar, B_bar, A_p, X_opt); % 一次初始化结果
% 把R_init_est1和R_init_est2，以及eR = norm(R_true - R_init_est1, 'fro')和eR = norm(R_true - R_init_est2, 'fro')存入R_beta_history

R_beta_history = {};
beta_vec = 0:0.1:5;
for beta = beta_vec
    R_beta = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init, mu, beta, R_true); % using R_init
    R_beta.eR = norm(R_true - R_beta.R, 'fro');
    R_beta_history{end+1} = R_beta;
    R_beta
end

% results由三部分组成：
% 1. R_beta_history
% 2. R_init_est1 and R_init_est2和R_true的误差等关系
% 3. R_true

R_init_est = struct('R_init_est1', R_init_est1, 'R_init_est2', R_init_est2, 'eR_init_est1', norm(R_true - R_init_est1, 'fro'), 'eR_init_est2', norm(R_true - R_init_est2, 'fro'));

results.num_beta = numel(beta_vec);
results.beta_vec = beta_vec;
results.R_beta_history = R_beta_history;
results.R_init_est = R_init_est;
results.R_true = R_true;

% % LM
% [p_lm, R_lm, stats_lm] = estimate_pose_lm( ...
%     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );

% % ELM
% [p_elm, R_elm, stats_elm] = estimate_pose_elm( ...
%     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );

% [p_Rlm, R_Rlm, stats_Rlm] = estimate_R_lm( ...
%     b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_true, options );
end

