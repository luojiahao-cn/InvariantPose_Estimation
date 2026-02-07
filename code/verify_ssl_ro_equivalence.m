%% 验证 SSL 与 RO 的等价性 (SSL-RO Equivalence Verification)
% 目标：验证当 SSL 的 beta 取动态值 beta_prox = 1/eta + r' * gradE 时，
%      其单步更新结果与黎曼梯度下降 (RO) 完全一致。

clc; clear; close all;
addpath('./Functions')
addpath('./utils')

% 1. 模拟一个简单的球面优化问题数据
% 目标函数: f(r) = r' A r - 2 g' r, s.t. ||r||=1
A = [2, 0.1, 0.5; 0.1, 3, -0.2; 0.5, -0.2, 1.5];
g_vec = [1; 2; 0.5];
r = [1; 0; 0]; % 初始点
r = r / norm(r);

% 定义欧几里得梯度计算函数
calc_gradE = @(x) 2*A*x - 2*g_vec;

% 2. 参数设置
eta = 0.5; % 设定一个步长

% --- RO (Riemannian Optimization) 更新步骤 ---
gradE = calc_gradE(r);
% 投影到切平面得到黎曼梯度
P = eye(3) - (r * r');
gradR = P * gradE;
% 沿负梯度移动并归一化 (Retraction)
r_ro_next = (r - eta * gradR) / norm(r - eta * gradR);

% --- SSL (Successive Semidefinite Lifting) 更新步骤 ---
% 使用推导的动态 beta: beta_prox = 1/eta + r' * gradE
% 注意：由于 SSL 默认公式是 m = -gradE + beta_prox * r
% 我们需要对应调整缩放尺度：
% RO 的分子是: (1 + eta * r' * gradE) * r - eta * gradE
% 将其除以 eta 得到更新方向为: (1/eta + r' * gradE) * r - gradE
beta_prox_dynamic = 1/eta + r' * gradE;
m_ssl = -gradE + beta_prox_dynamic * r;
r_ssl_next = m_ssl / norm(m_ssl);

% 3. 结果对比
fprintf('--- 单步更新验证 (eta = %.2f) ---\n', eta);
fprintf('RO  更新后的向量: [%.6f, %.6f, %.6f]\n', r_ro_next);
fprintf('SSL 更新后的向量: [%.6f, %.6f, %.6f]\n', r_ssl_next);
fprintf('两者之间的二范数距离: %.2e\n', norm(r_ro_next - r_ssl_next));

% 4. 多步迭代验证
fprintf('\n--- 10步迭代验证 ---\n');
r_ro = r;
r_ssl = r;
for i = 1:10
    % RO step
    gE_ro = calc_gradE(r_ro);
    gR = (eye(3) - r_ro*r_ro') * gE_ro;
    r_ro = (r_ro - eta * gR) / norm(r_ro - eta * gR);
    
    % SSL step with dynamic beta
    gE_ssl = calc_gradE(r_ssl);
    bpD = 1/eta + r_ssl' * gE_ssl;
    m = -gE_ssl + bpD * r_ssl;
    r_ssl = m / norm(m);
    
    dist = norm(r_ro - r_ssl);
    fprintf('第 %2d 步距离: %.2e\n', i, dist);
end

if dist < 1e-12
    fprintf('\n验证成功：SSL 在动态 beta 设置下等价于黎曼梯度下降 (RO)！\n');
else
    fprintf('\n验证失败，请检查推导细节。\n');
end
