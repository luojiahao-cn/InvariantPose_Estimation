clc,clear,close all
addpath('./utils')
addpath('./Functions')
addpath('./tools')

%% ========== 参数配置 ==========
params = get_experiment_params();
R_true = MatrixExp3(VecToso3(params.ground_truth.theta_true));
p_true = params.ground_truth.p_true;
pC = params.magnet.m_pos;   % 3xM
M  = size(pC,2);
% 设置随机种子
rng(params.experiment.random_seed);
d_list = [
    [1e-3; 0; 0], ...
    [0; 0; 0], ...
    [0; 0; 1e-3]
];
d_list = d_list - mean(d_list, 2); % 去中心化
params.sensor.d_list = d_list;
%% ======== 执行实验 ==========
% 当前测试生成磁场数据
b_total = generate_magnetic_data(params);
num_sensors = size(b_total, 2);
% 设定时分增加幅值：
m_norm_variation = 5; % 幅值变化范围的比例

for i = 1:3
    params.magnet.m_norm(i) = params.magnet.m_norm(i) + m_norm_variation;

    [b_total_varied, ~, gradb_total] = generate_magnetic_data(params);

    b_difference = b_total_varied - b_total;

    % ======== 计算b_mag ========
    r = rank(d_list');
    Q_bar = null(d_list');   % 保证 d_list' * Q_bar = 0

    g = Q_bar' * ones(num_sensors, 1);
    b_mag = b_difference * Q_bar * g / norm(g)^2; % 计算b_mag
    %% ========= 计算梯度场 ==========
    X_opt = lc_grad_tensor_estimator(b_difference, params.sensor.d_list);

    rho(:,i) = -3 * (X_opt \ b_mag);

    params.magnet.m_norm(i) = params.magnet.m_norm(i) - m_norm_variation;
end

%% ======== 计算(R, p) =========

% 这里“源”的数量等于你循环里 i=1:3 的数量
% pC: 3xM  (世界系下每个驱动源的位置)
% rho: 3xM (由单源测得的 rho_i = R^T (p - pC_i) )
assert(size(rho,2) == M, 'rho 列数必须与源数量一致');
assert(size(pC,1) == 3 && size(rho,1) == 3, 'pC/rho 必须是 3xM');

% 传感器系点：s_i = -rho_i，使得  pC_i = p + R*s_i
S = -rho;     % 3xM  (sensor frame “points”)
W = pC;       % 3xM  (world frame points)

% --- 用“互差向量”构造 Wahba / Kabsch ---
pairs = nchoosek(1:M, 2);
K = size(pairs,1);

Ds = zeros(3, K);  % sensor-frame differences
Dw = zeros(3, K);  % world-frame differences

for k = 1:K
    i = pairs(k,1);
    j = pairs(k,2);
    Ds(:,k) = S(:,i) - S(:,j);   % d_s = s_i - s_j
    Dw(:,k) = W(:,i) - W(:,j);   % d_w = w_i - w_j
end

% 可选：给每个差向量加权（更鲁棒，尤其某些源远/噪声大时）
% 这里默认等权。若想加权，可用 w(k)=norm(Dw(:,k)) 或其他置信度指标
w = ones(1,K);

% 形成协方差矩阵 H = sum_k w_k * d_s_k * d_w_k^T
H = zeros(3,3);
for k = 1:K
    H = H + w(k) * (Ds(:,k) * Dw(:,k).');
end

% SVD 求解旋转：Dw ≈ R * Ds  =>  H = UΣV^T, R = V D U^T
[U,~,V] = svd(H);
D = diag([1, 1, det(V*U')]);     % 保证 det(R)=+1
R_est = V * D * U';

% --- 用 R 估计 p：p = mean_i (w_i - R*s_i) ---
p_candidates = zeros(3, M);
for i = 1:M
    p_candidates(:,i) = W(:,i) - R_est * S(:,i);
end
p_est = mean(p_candidates, 2);

% --- 自检残差 ---
res = zeros(1,M);
for i = 1:M
    res(i) = norm(W(:,i) - (p_est + R_est*S(:,i)));
end

fprintf('Estimated R:\n'); disp(R_est);
fprintf('Estimated p:\n'); disp(p_est);
fprintf('Per-source residual norms:\n'); disp(res);
fprintf('Mean residual: %.3e\n', mean(res));

%% 对比R_true, p_true
norm(R_true - R_est, 'fro')
norm(p_true - p_est)

% 你也可以把结果写回 params 或输出
% params.est.R = R_est;
% params.est.p = p_est;
