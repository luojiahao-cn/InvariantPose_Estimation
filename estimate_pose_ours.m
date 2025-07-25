function [p_est_22, R_est, stats] = estimate_pose_ours(b_total, d_list, m_pos, m_hat, m_norm, p_init)
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
% b_p = zeros(3, 1);
% A_p = zeros(3, 3);
% for i = 1:size(m_pos,2)
%     r = p_init - m_pos(:, i);
%     [B, gradB] = dipole_b_and_gradb(r, m_hat(:, i), m_norm(i));
%     b_p = b_p + B;
%     A_p = A_p + gradB;
% end
% norm(R_true'*b_p*ones(1,num_sensors)+R_true'*A_p*R_true*d_list - b_total, 'fro')
% norm(R_true'*A_p*R_true*D_matrix - B_matrix, 'fro')
% X_true = R_true'*A_p*R_true;
%% 步骤2: 构建梯度约束方程
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

%% 步骤3: 位置估计第一阶段（公式20）
% 对d_list'进行QR分解
[Q, ~] = qr(d_list');
r = rank(d_list); 
Q_bar = Q(:, r+1:end);
bQ_bar = b_total * Q_bar;

% 优化第一阶段位置
fun20 = @(p) obj_fun20(p, m_pos, m_hat, m_norm, bQ_bar, Q_bar);
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off', ...
    'MaxIterations', 200, ...
    'FunctionTolerance', 1e-8, ...
    'StepTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 3000);
p_est = lsqnonlin(fun20, p_init, [], [], options);

%% 步骤4: 位置估计第二阶段（公式22）
fun22 = @(p) diag([1,0.5,0.5])*obj_fun22(p, m_pos, m_hat, m_norm, bQ_bar, Q_bar, X_opt);
p_est_22 = lsqnonlin(fun22, p_init, [], [], options);

%% 步骤5: 方向估计（公式24-25）
% 计算在p_est_22处的全局坐标系梯度矩阵
A_p_opt = zeros(3,3);
% p_est_22 = p_init;
for magnet_idx = 1:size(m_pos,2)
    r = p_est_22 - m_pos(:, magnet_idx);
    [~, gradB] = dipole_b_and_gradb(r, m_hat(:, magnet_idx), m_norm(magnet_idx));
    A_p_opt = A_p_opt + gradB;
end

% 特征分解对齐
[U_A, Lambda_A] = eig(A_p_opt);
[~, idx_A] = sort(diag(Lambda_A), 'descend');
U_A = U_A(:,idx_A);

[V_X, Lambda_X] = eig(X_opt);
[~, idx_X] = sort(diag(Lambda_X), 'descend');
V_X = V_X(:,idx_X);

% 计算旋转矩阵估计
R_sgn = [1,1,1;
    1,1,-1;
    1,-1,1;
    -1,1,1];

for i = 1:4
    R_opt(:,:,i) = U_A * diag(R_sgn(i,:)) * V_X';
    if det(R_opt(:,:,i)) < 1
        R_opt(:,:,i) = - R_opt(:,:,i);
    end
    theta(i) = norm(so3ToVec(MatrixLog3(R_opt(:,:,i))));
end
[~, ind] = min(theta);
R_est = R_opt(:,:,ind);

% 特征分解对齐
% [U_A, Lambda_A] = eig(A_p_opt);
% [V_X, Lambda_X] = eig(X_opt);

% [~, idx_A] = sort(diag(Lambda_A), 'descend');
% [~, idx_X] = sort(diag(Lambda_X), 'descend');

% U_A = U_A(:, idx_A);
% V_X = V_X(:, idx_X);

% 计算旋转矩阵估计
% R_est = U_A * V_X';

% % 确保SO(3)
% if det(R_est) < 0
%     U_A(:, 3) = -U_A(:, 3);
%     R_est = U_A * V_X';
% end

%% 保存中间结果
stats.X_opt = X_opt;         % 估计的梯度矩阵
stats.A_p_opt = A_p_opt;     % 参考点的全局梯度矩阵
stats.p_est = p_est;          % 第一阶段位置估计
stats.BQ_bar = bQ_bar;        % BQ_bar矩阵
stats.Q_bar = Q_bar;          % 零空间矩阵

end

%% 目标函数（公式20）
function f = obj_fun20(p, m_pos, m_hat, m_norm, bQ_bar, Q_bar)
    b_p = zeros(3, 1);
    for i = 1:size(m_pos,2)
        r = p - m_pos(:, i);
        b_i = dipole_b_and_gradb(r, m_hat(:, i), m_norm(i));
        b_p = b_p + b_i(:, 1);
    end
    
    bQ_norm = norm(b_p * ones(1, size(Q_bar,1)) * Q_bar, 'fro');
    bQ_bar_norm = norm(bQ_bar, 'fro');
    
    f = (bQ_norm - bQ_bar_norm)^2;
end

%% 目标函数（公式22）
function res = obj_fun22(p, m_pos, m_hat, m_norm, bQ_bar, Q_bar, X_opt)
    b_p = zeros(3, 1);
    A_p = zeros(3, 3);
    for i = 1:size(m_pos,2)
        r = p - m_pos(:, i);
        [B, gradB] = dipole_b_and_gradb(r, m_hat(:, i), m_norm(i));
        b_p = b_p + B;
        A_p = A_p + gradB;
    end
    
    % 公式(22)的三个残差项
    term1 = norm(b_p * ones(1, size(Q_bar,1)) * Q_bar, 'fro') - norm(bQ_bar, 'fro');
    term2 = trace(A_p*A_p) - trace(X_opt*X_opt);
    term3 = det(A_p) - det(X_opt);
    
    res = [term1; term2; term3];
end

