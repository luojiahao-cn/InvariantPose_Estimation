function [p_est, R_est3, stats] = estimate_pose_ours(b_total, d_list, m_pos, m_hat, m_norm, p_init, R_true, p_true, options)
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
%% 步骤1: 构建磁场差矩阵和位移矩阵
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
[~, A_p] = calcFieldAndGradient(p_true, m_pos, m_hat, m_norm);
X_true = R_true'*A_p*R_true;
% norm(R_true'*b_p*ones(1,num_sensors)+R_true'*A_p*R_true*d_list - b_total, 'fro')
% norm(R_true'*A_p*R_true*D_matrix - B_matrix, 'fro')
%% 步骤2: 构建梯度约束方程，计算Q_bar
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

% 对d_list'进行QR分解
[Q, ~] = qr(d_list');
r = rank(d_list); 
Q_bar = Q(:, r+1:end);
%% 步骤3: 位置估计第一阶段（公式20）
% 优化第一阶段位置
fun22 = @(p) obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_total, Q_bar, X_opt);
p_est = lsqnonlin(fun22, p_init, [], [], options);
%% 步骤4: 方向估计（公式24-25）
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B = b_total * Q_bar;
C = b_p * ones(1, num_sensors) * Q_bar;
[R_est, R_est2] = estiamteR(p_est, m_pos, m_hat, m_norm, num_sensors, b_total, Q_bar, X_opt);
% norm(R_est-R_true,'fro')
% norm(R_est-R_est2,'fro')
mu = norm(B,'fro')*norm(C,'fro')/(norm(X_opt, 'fro')*norm(A_p, 'fro'));
R_est3 = estimateR_iter(p_est, m_pos, m_hat, m_norm, num_sensors, b_total, Q_bar, X_opt, mu);
% R_est4 = softProcrustesGrad(R_est,B,C,A_p,X_opt,mu);
% R_est5 = softProcrustesLM(R_est,B,C,A_p,X_opt,norm(B*C','fro')/norm(A_p*X_opt,'fro'));
% clc
% norm(R_true - R_est, 'fro')
% norm(R_true - R_est2, 'fro')
% norm(R_true - R_est3, 'fro')
% norm(R_true - R_est4, 'fro')
% norm(R_true - R_est5, 'fro')
%% 保存中间结果
stats.X_opt = X_opt;         % 估计的梯度矩阵
end

%% ----------------------------Functions-------------------------------  %% 
%% calc magntic field and gradient tensor
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
%% estimate p
function res = obj_fun22(p, m_pos, m_hat, m_norm, num_sensors, b_total, Q_bar, X)
    [b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);
    % 公式(22)的三个残差项
    term1 = norm(b_p * ones(1, num_sensors) * Q_bar, 'fro') - norm(b_total * Q_bar, 'fro');
    term2 = trace(A_p*A_p) - trace(X*X);
    term3 = det(A_p) - det(X);
    res = [term1; term2; term3];
end
%% Estimate R
function R_opt = estimateR_iter(p_est, m_pos, m_hat, m_norm, num_sensors, b_total, Q_bar, X, mu)
    %% Initial condition
    [b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
    [U, LA] = eig(A_p);
    [V, LX] = eigs(X);
    [~, Ind] = sort(diag(LA), 'ascend');
    U = U(:,Ind);
    [~, IndX] = sort(diag(LX), 'ascend');
    V = V(:,IndX);
    R = U*diag(sign(diag(V'*U)))*V';
    if det(R) < 0
        R = -R;
    end
    %% Iterative method:
    skew = @(X) 0.5*(X-X');
    % tr(BC^T R + muARXR^T) = tr((BC^T + muXR^TA)R) = tr(MR)
    B = b_total * Q_bar * (b_p*ones(1,num_sensors)*Q_bar)';
    f = @(R) -trace((B+mu*X*R'*A_p)*R);
    for i = 1:20
        R_opt = R;
        % Procrustes
        M = B + mu * X * R' * A_p;
        [U, ~, V] = svd(M);
        R = V*diag([1,1,det(V*U')])*U'; % R_{k+1}^{-}
        Omega = MatrixLog3(R*R_opt');
        % Armijo backtracking method
        alpha = 1;
        cArmijo = 0.1;
        gradf = 2*skew(R'*(B'+2*mu*A_p*R*X));
        while(f(R_opt*MatrixExp3(alpha*Omega)) >= (f(R_opt) + cArmijo*alpha*1/2*trace(gradf'*Omega)))
            tau = 0.5;
            alpha = tau * alpha;
            if alpha <= 1e-4
                break
            end
        end
        R = R_opt*MatrixExp3(alpha*Omega);

        if norm(R - R_opt, 'fro') < 0.01
            R_opt = R;
            disp('converged')
            break
        end
    end
end

function [R_opt1, R_opt2] = estiamteR(p_est, m_pos, m_hat, m_norm, num_sensors, b_total, Q_bar, X)
    [b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
    % estimate R using R^T A R = X
    [U, LA] = eig(A_p);
    [V, LX] = eigs(X);
    [~, Ind] = sort(diag(LA), 'ascend');
    U = U(:,Ind);
    [~, IndX] = sort(diag(LX), 'ascend');
    V = V(:,IndX);
    R_opt1 = U*diag(sign(diag(V'*U)))*V';
    if det(R_opt1) < 0
        R_opt1 = -R_opt1;
    end
    % estimate R using SVD
    M = b_total * Q_bar * (b_p*ones(1,num_sensors)*Q_bar)';
    [U, ~, V] = svd(M);
    R_opt2 = V*diag([1,1,det(V*U')])*U';
    
    % norm(b_total * Q_bar - R_opt1' * (b_p*ones(1,num_sensors)*Q_bar))
    % norm(b_total * Q_bar - R_opt2' * (b_p*ones(1,num_sensors)*Q_bar))

end

function R = softProcrustesGrad(R0,B,C,A,X,mu,maxIter,tol)
% softProcrustesGrad   Riemannian gradient descent on SO(3)
%
%   R = softProcrustesGrad(B,C,A,X,mu)
%   R = softProcrustesGrad(B,C,A,X,mu,maxIter,tol)
%
% INPUT
%   B,C,A,X : 3×3 real matrices
%   mu      : positive scalar, penalty weight
%   maxIter : (optional) maximum iterations, default 200
%   tol     : (optional) gradient‑norm tolerance, default 1e‑8
%
% OUTPUT
%   R       : 3×3 rotation matrix (det=+1)

    if nargin < 7,  maxIter = 200;  end
    if nargin < 8,  tol     = 1e-8; end

    % ---------- helper lambdas ----------
    hat  = @(v) [   0   -v(3)  v(2);
                  v(3)    0   -v(1);
                 -v(2)  v(1)    0  ];        % so(3) "hat"
    vee  = @(S) [S(3,2); S(1,3); S(2,1)];     % inverse‑hat
    skew = @(M) 0.5*(M-M.');                  % skew‑symmetric part
    cost = @(R) norm(B-R.'*C,'fro')^2 + ...
                mu*norm(X-R.'*A*R,'fro')^2;   % objective J(R)
    % ---------- 0. initialiation -----------
    R = R0;
    % ---------- 1. main iteration ----------
    cArmijo = 0.1;       beta = 0.5;      % line‑search params
    for k = 1:maxIter
        % 1.1 Riemannian gradient (body‑frame) --------------------
        G  = skew(R.'*C*B.') + ...
              2*mu*skew(X*R.'*A - A*R*X);     % Eq.(2) in说明
        g  = vee(G);                          % 3×1 vector form
        ng = norm(g);

        if ng < tol,   break;  end            % convergence test

        % 1.2 search direction (steepest descent) ----------------
        dir = -g;                             % δ_k  (3×1)

        % 1.3 Armijo back‑tracking line search -------------------
        alpha = 1.0;
        Jcur  = cost(R);
        while true
            Rnew = R * expm(hat(alpha*dir));  % re‑project via exp
            if cost(Rnew) <= Jcur - cArmijo*alpha*ng^2
                break;                        % sufficient decrease
            end
            alpha = alpha*beta;
            if alpha < 1e-12
                warning('Step size underflow — stop.'); 
                R = Rnew;  return;
            end
        end

        % 1.4 accept step ---------------------------------------
        R = Rnew;
    end
end

function R = softProcrustesLM(R0,B,C,A,X,mu,maxIter,tol,lambda0)
% softProcrustesLM   Levenberg–Marquardt solver on SO(3) for
%                    J(R)=‖B-RᵀC‖_F² + μ‖X-RᵀAR‖_F².
%
%   R = softProcrustesLM(B,C,A,X,mu)
%   R = softProcrustesLM(B,C,A,X,mu,maxIter,tol,lambda0)
%
% INPUT
%   B,C,A,X   3×3 real matrices (B,C 任意；A,X 对称更典型)
%   mu        惩罚权重 μ>0
%   maxIter   (可选) 最大迭代步, 默认 100
%   tol       (可选) ‖δ‖ 收敛阈值, 默认 1e‑10
%   lambda0   (可选) 初始阻尼 λ₀,  默认 1e‑3
%
% OUTPUT
%   R         3×3 旋转矩阵 (det≈+1)

    if nargin<7, maxIter = 100; end
    if nargin<8, tol     = 1e-10; end
    if nargin<9, lambda0 = 1e-3; end

    %% -------- helper lambdas --------
    hat  = @(v) [  0   -v(3)  v(2);
                  v(3)   0   -v(1);
                 -v(2)  v(1)   0 ];
    vee  = @(S) [S(3,2); S(1,3); S(2,1)];
    skew = @(M) 0.5*(M-M.');
    % 目标函数
    cost = @(R) norm(B-R.'*C,'fro')^2 + ...
                mu*norm(X-R.'*A*R,'fro')^2;

    %% -------- 0. 初始R：一次SVD-Procrustes --------
    R = R0;

    %% -------- 1. LM 主循环 --------
    lambda = lambda0;
    nuInc  = 2;     % 增大λ的倍率
    nuDec  = 0.5;   % 减小λ的倍率

    for k = 1:maxIter
        % ---- 1.1 残差 r ∈ R^18 ---------------------
        r1 = (B - R.'*C);                 % 9×1
        r2 = (X - R.'*A*R);               % 9×1
        r1 = r1(:);
        r2 = r2(:);
        r  = [r1; sqrt(mu)*r2];               % 18×1

        % ---- 1.2 雅可比 J ∈ R^{18×3} --------------
        % 对误差项  E = B - RᵀC  有  dE/dδ = -(CR)×
        J1 = -crossMat(R'*C);                  % 9×3
        % 对误差项  F = X - RᵀAR 有  dF/dδ = -(AR - RXRᵀA)×
        J2 = -crossMat(A*R - R*X);            % 9×3
        J  = [J1; sqrt(mu)*J2];               % 18×3

        % ---- 1.3 LM 增量：解 (JᵀJ+λI)δ = -Jᵀr -------
        H     = J.'*J + lambda*eye(3);
        delta = -H \ (J.'*r);                 % 3×1
        dNorm = norm(delta);

        % ---- 1.4 收敛检测 ---------------------------
        if dNorm < tol,  break;  end

        % ---- 1.5 执行更新并评估 ----------------------
        Rnew = R * expm(hat(delta));
        if cost(Rnew) < cost(R)          % 成功，接受步长
            R      = Rnew;
            lambda = max(lambda*nuDec, 1e-12);
        else                              % 失败，加大阻尼
            lambda = lambda*nuInc;
        end
    end
end

%% -------- 附：生成交叉乘矩阵 (v×·) 的辅助函数 --------
function J = crossMat(M)
% 将 3×3 矩阵 M 映射到 9×3 Jacobian:  (Mδ) vec = (M ×) δ
    J = [  0     -M(3,1)  M(2,1);
           0     -M(3,2)  M(2,2);
           0     -M(3,3)  M(2,3);
          M(3,1)  0     -M(1,1);
          M(3,2)  0     -M(1,2);
          M(3,3)  0     -M(1,3);
         -M(2,1)  M(1,1)  0;
         -M(2,2)  M(1,2)  0;
         -M(2,3)  M(1,3)  0 ];
end
