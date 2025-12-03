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
pairs = nchoosek(1:num_sensors, 2);
D_matrix = zeros(3, size(pairs,1));
B_matrix = zeros(3, size(pairs,1));

for idx = 1:size(pairs,1)
    i = pairs(idx, 1);
    j = pairs(idx, 2);
    D_matrix(:, idx) = d_list(:, j) - d_list(:, i);
    B_matrix(:, idx) = b_total(:, j) - b_total(:, i);
end
%% 阶段检查 对应公式2
% [~, A_p_true] = calcFieldAndGradient(params.p_true, m_pos, m_hat, m_norm);
% R_true = params.R_true;
% X_true = R_true'*A_p_true*R_true;
% norm(R_true'*b_p*ones(1,num_sensors)+R_true'*A_p*R_true*d_list - b_total, 'fro')
% norm(R_true'*A_p*R_true*D_matrix - B_matrix, 'fro')
%% 估计局部梯度张量
% 构建选择矩阵S
S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
    0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
    0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

% 构建完整约束矩阵C
C_matrix = kron(D_matrix', eye(3)) * S;
h_vector = B_matrix(:);

% 求解最小二乘问题
x_opt = pinv(C_matrix) * h_vector;
X_opt = reshape(S * x_opt, 3, 3);  % 估计梯度（传感器坐标系）
%% Stage #1: Estimate for position \hat{p}
% 对d_list'进行QR分解
[Q, ~] = qr(d_list');
r = rank(d_list); % 构型判据
Q_bar = Q(:, r+1:end);
b_bar = b_total * Q_bar; % 计算bBar
% 优化第一阶段位置
fun22 = @(p) obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt);
% Jp = numJacobian(fun22, p_true)
% rank(Jp)
% svd(Jp)
% c = cond(Jp)
p_est = lsqnonlin(fun22, p_init, lb_p, ub_p, options);
%% Stage #2: Estimate for rotation \hat{R}
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B_matrix = b_p * ones(1, num_sensors);
B_bar = B_matrix * Q_bar;
[R_init_est1, R_init_est2] = estimateR(b_bar, B_bar, A_p, X_opt);

L = norm(A_p)*norm(X_opt)/(norm(B_bar)*norm(b_bar));
% L = 1000;
% L = norm(A_p)*norm(X_opt);
% mu = L; %L; % regularization parameter
% mu = 1e3;
% beta = 1e-3;
% 使用外部beta变量

R_est = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init, params.mu, params.beta); % using R_init
% R_est = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init_est2, params.mu, params.beta); % using R_init_est
stats.X_opt = X_opt;         % 估计的梯度矩阵
stats.R_iter_history = R_est.R_iter_history; % 每次迭代的R
stats.delta_history = R_est.delta_history;   % 每次迭代的delta
R_est = R_est.R; % 返回最终R
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

function J = numJacobian(fun, x)
f0 = fun(x);
n  = numel(x);
m  = numel(f0);
J  = zeros(m, n);
h  = 1e-6;  % 步长，可调

for i = 1:n
    xh      = x;
    xh(i)   = xh(i) + h;
    J(:, i) = (fun(xh) - f0) / h;  % 前向差分
end
end
