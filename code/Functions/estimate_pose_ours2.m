function [p_est, R_est, stats] = estimate_pose_ours(b_total, d_list, m_pos, m_hat, m_norm, p_init, R_init, R_true, p_true, options, lb_p, ub_p, varargin)
% PROPOSED_METHOD_POSE_ESTIMATION 使用所提方法估计传感器姿态（位置和方向）
%
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

num_sensors = size(b_total, 2);

% --------- 支持可选beta参数 ---------
if ~isempty(varargin) && isnumeric(varargin{1})
    beta = varargin{1};
else
    beta = 1e-3;
end

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
[~, A_p_true] = calcFieldAndGradient(p_true, m_pos, m_hat, m_norm);
X_true = R_true'*A_p_true*R_true;
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
lb = lb_p(:);
ub = ub_p(:);
fun22 = @(p) obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt);
% Jp = numJacobian(fun22, p_true)
% rank(Jp)
% svd(Jp)
% c = cond(Jp)
p_est = lsqnonlin(fun22, p_init, lb, ub, options);
%% Stage #2: Estimate for rotation \hat{R}
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B_p = b_p * ones(1, num_sensors);
B_bar = B_p * Q_bar;
[~, R_init_est] = estiamteR(b_bar, B_bar, A_p, X_opt);
% R_est = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init_est, beta);
R_est = estimateR_iter(A_p, B_p, b_total, d_list, R_init);
stats.X_opt = X_opt;         % 估计的梯度矩阵
stats.R_iter_history = R_est.R_iter_history; % 每次迭代的R
stats.delta_history = R_est.delta_history;   % 每次迭代的delta
R_est = R_est.R; % 返回最终R
end

%% ----------------------------Functions-------------------------------  %%
%% Calc magntic field and gradient tensor
function [b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm)
b_p = zeros(3, 1);
A_p = zeros(3, 3);
for i = 1:size(m_pos,2)
    r = p - m_pos(:, i);
    [B, gradB] = dipole_b_and_gradb(r, m_hat(:, i), m_norm(i));
    b_p = b_p + B;
    A_p = A_p + gradB;
end
end
%% Estimate p
function res = obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X)
[b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);
term1 = norm(b_p * ones(1, num_sensors) * Q_bar, 'fro') - norm(b_bar, 'fro');
term2 = trace(A_p*A_p) - trace(X*X);
term3 = det(A_p) - det(X);
res = [term1; term2; term3];
end
%% Estimate R
function R_struct = estimateR_iter(A_p, B_p, b_list, d_list, R_init)
%% Iterative method:
A_bar = A_p;
B_bar = d_list * b_list' + b_list * d_list';

C_bar = A_p * A_p;
D_bar = d_list * d_list';

E_bar = 2 * (B_p * b_list' - A_p * B_p * d_list');

A_t = A_bar + min(eig(A_bar)) * eye(3);
B_t = B_bar + min(eig(B_bar)) * eye(3);

% L = 1000;
% L = norm(A_p)*norm(X_opt);
% mu = L; %L; % regularization parameter
% mu = 1e6;
% 使用外部beta变量
% M = @(R) mu * B_bar * b_bar' + 2 * At * R * Xt + beta * R;
M = @(R) 2 * A_t * R * B_t - 2 * C_bar * R * D_bar + E_bar;
% f = @(R) trace(R' * A_p * R * X_opt + mu * R' * B_bar * b_bar' + beta * R);
% fbar = @(R) trace(R' * At * R * Xt + mu * R' * B_bar * b_bar' + beta * R);
% skw = @(R) 0.5 * (R - R');
% g = @(R) 2 * skw(R' * (2 * A_p * R * X_opt + mu * B_bar * b_bar' + beta * eye(3)));
% gbar = @(R) 2 * skw(R' * (2 * At * R * Xt + mu * B_bar * b_bar' + beta * eye(3)));
%% Initial condition
delta = 1e6; 
k = 0; kmax = 1000;
R = R_init;
R_iter_history = {};
delta_history = [];
while k < kmax && delta > 1e-5
    [U, ~, V] = svd(M(R));
    R_opt = U*diag([1,1,det(U*V')])*V';
    delta = norm(R_opt - R, 'fro');
    R_iter_history{end+1} = R_opt;
    delta_history(end+1) = delta;
    R = R_opt;
    k = k + 1;
end
if k == kmax
    warning('迭代未收敛，可能需要调整参数');
end
R_struct.R = R;
R_struct.R_iter_history = R_iter_history;
R_struct.delta_history = delta_history;
end

function [R_opt1, R_opt2] = estiamteR(b_bar, B_bar, A_p, X)
% estimate R using R^T A R = X
[PA, LA] = eig(A_p);
[PX, LX] = eig(X);
[lam, IndA] = sort(diag(LA), 'ascend');
PA = PA(:,IndA);
[sig, IndX] = sort(diag(LX), 'ascend');
PX = PX(:,IndX);

% ---- 规范为 SO(3) 的特征向量基 ----
if det(PA) < 0, PA(:,1) = -PA(:,1); end
if det(PX) < 0, PX(:,1) = -PX(:,1); end

% ---- 枚举 6 个置换，做最佳匹配 ----
perms3 = perms(1:3);
bestCost = inf; bestP = eye(3);
for k = 1:size(perms3,1)
    p = perms3(k,:);
    P = eye(3); P = P(:, p);                 % 置换矩阵
    cost = sum((lam(p) - sig).^2);           % 等价于最大化 sum(sig .* lam(p))
    if cost < bestCost
        bestCost = cost; bestP = P;
    end
end

% ---- 保证 Y ∈ SO(3) ----
Y = bestP;
if det(Y) < 0
    Y = Y * diag([1 1 -1]);                  % 用符号阵把行列式调成 +1
end

% ---- 回代得到最优解 R ----
R_opt1 = PA * Y * PX';
% R_opt1 = PA * diag(sign(diag(PX'*PA)))*PX';

% estimate R using SVD
M = B_bar * b_bar';
[U, ~, V] = svd(M);
R_opt2 = V*diag([1,1,det(V*U')])*U';
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
