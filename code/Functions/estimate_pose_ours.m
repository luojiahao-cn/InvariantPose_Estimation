function [p_est, R_est, stats] = estimate_pose_ours(b_total, d_list, m_pos, m_hat, m_norm, p_init, R_true, p_true, options, lb_p, ub_p)
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
p_est = lsqnonlin(fun22, p_init, lb, ub, options);
%% Stage #2: Estimate for rotation \hat{R}
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B_bar = b_p * ones(1, num_sensors) * Q_bar;
R_est = estimateR_iter(b_bar, B_bar, A_p, X_opt);
R_true
%% 保存中间结果
stats.X_opt = X_opt;         % 估计的梯度矩阵
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
function R = estimateR_iter(b_bar, B_bar, A_p, X)
    %% Iterative method:
    LA = min(eig(A_p));
    LX = min(eig(X));
    At = A_p - LA * eye(3);
    Xt = X - LX * eye(3);
    L = 4*norm(At)*norm(Xt);
    mu = 0; %1/L; % regularization parameter
    beta = 0; % regularization parameter for R
    M = @(R) mu * B_bar * b_bar' + 2 * At * R * Xt + beta * R;
    f = @(R) trace(R' * A_p * R * X + mu * R' * B_bar * b_bar' + beta * R);
    fbar = @(R) trace(R' * At * R * Xt + mu * R' * B_bar * b_bar' + beta * R);
    skw = @(R) 0.5 * (R - R');
    %% Initial condition
    R = estiamteR(b_bar, B_bar, A_p, X);
    f(R) 
    fbar(R) + LA * trace(X) + LX * trace(A_p) - 3 * LA * LX

    skw(R'*(A_p*R*X))
    delta = 1e6; k = 0; kmax = 20;
    while k < kmax && delta > 1e-6
        [U, Sigma, V] = svd(M(R));
        R_opt = U*diag([1,1,det(U*V')])*V';
        delta = norm(R_opt - R, 'fro');
        R = R_opt;
        k = k + 1;
    end
end

function [R_opt1, R_opt2] = estiamteR(b_bar, B_bar, A_p, X)
    % estimate R using R^T A R = X
    [PA, LA] = eig(A_p);
    [PX, LX] = eigs(X);
    [~, IndA] = sort(diag(LA), 'ascend');
    PA = PA(:,IndA);
    [~, IndX] = sort(diag(LX), 'ascend');
    PX = PX(:,IndX);
    Y = PA*diag(sign(diag(PX'*PA)))*PX';
    [U, ~, V] = svd(Y);
    R_opt1 = U*diag([1,1,det(U*V')])*V';
    % S = diag(sign(diag(PX'*PA)));
    % if det(S)*det(PA)*det(PX) < 0
    %     d = diag(PX'*PA);
    %     [~, ind] = min(abs(d));
    %     S(ind, ind) = -S(ind, ind); % 确保行列式为正
    % end
    % R_opt1 = PA*S*PX';
    % estimate R using SVD
    M = B_bar * b_bar';
    [U, ~, V] = svd(M);
    R_opt2 = V*diag([1,1,det(V*U')])*U';
end