function result = run_single_experiment(exp_idx, num_experiments, params, ...
    b_total, d_list, p_true, R_true, lb_p, ub_p)
% RUN_SINGLE_EXPERIMENT 执行单次实验
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
options = params.optimization.options;
mu = params.optimization.mu;
beta = params.optimization.beta;
W = params.optimization.W;
% ==== 初始值猜测 ====

% 初始化p_init为在以p_true为球心，半径为p_uncertainty，p_init(3)>=0的球内随机采样
% 球面内均匀采样，生成半径为 <= p_uncertainty 的向量
u = randn(3,1);       % 随机方向
u(3) = abs(u(3));
u = u / norm(u);      % 单位向量
r = (rand())^(1/3) * p_uncertainty;  % 半径分布
p_init = p_true + r * u;

u_theta = randn(3,1);
u_theta = u_theta / norm(u_theta);
r_theta = (rand())^(1/3) * r_uncertainty;
theta_init = theta_true + r_theta * u_theta;

fprintf('theta_init = [%.4f, %.4f, %.4f]\n', theta_init(1), theta_init(2), theta_init(3));
fprintf('p_init = [%.4f, %.4f, %.4f]\n', p_init(1), p_init(2), p_init(3));

% ==== 算法调用 ====
% LM
%% 计算方法的执行时间
tic;
[p_lm, R_lm, stats_lm] = estimate_pose_lm( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
time_lm = toc;

tic;
% ELM
[p_elm, R_elm, stats_elm] = estimate_pose_elm( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );
time_elm = toc;

% 所提算法
tic;
alg_params = struct('R_true', R_true, 'p_true', p_true, 'mu', mu, 'beta', beta, 'W', W);
[p_ours, R_ours, stats_ours] = estimate_pose_ours( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, alg_params);
time_ours = toc;

% ==== 直接估计主轴研究 (基于ours位置) ====

% 初始值设为离散枚举的结果
opts_r = struct('u', [0;0;1]);
opts_r.r0 = MatrixExp3(VecToso3(theta_init)) * opts_r.u; % 旋转初始主轴方向
opts_r.r_true = R_true * opts_r.u; % 真实主轴方向
opts_r.alpha = [1, 0];
opts_r.beta = 1;
[prob_r, opts_r] = compute_principal_axis_prob(stats_ours.b_bar, stats_ours.B_bar, stats_ours.A_p, stats_ours.X_opt, opts_r);

r_direct_spec = estimate_principal_axis(stats_ours.b_bar, stats_ours.B_bar, stats_ours.A_p, stats_ours.X_opt, opts_r);
[r_direct_MM, info_r_MM]   = estimate_principal_axis_SCA(prob_r, opts_r);
[r_direct_SDP, info_r_SDP] = estimate_principal_axis_SDP(prob_r, opts_r);
[r_direct_RO, info_r_RO]   = estimate_principal_axis_RO(prob_r, opts_r);
% prob_r.r0 = r_direct_spec; % 将 SPEC 结果作为 SSL 的初始值
[r_direct_SSL, info_r_SSL] = estimate_principal_axis_SSL(prob_r, opts_r);


% mu = 1e-6;
% % h = @(r) r' * stats_ours.A_p * r + mu * r' * stats_ours.B_bar * stats_ours.b_bar' * opts_r.u;
% alpha = [1, 0];
% f = @(r) mu * norm(stats_ours.B_bar' * r - stats_ours.b_bar' * opts_r.u)^2;
% h = @(r) r' * stats_ours.A_p * r + f(r);
% s1 = opts_r.u' * stats_ours.X_opt * opts_r.u;
% s2 = opts_r.u' * stats_ours.X_opt * stats_ours.X_opt * opts_r.u;
% g = @(r) alpha(1) * (r' * stats_ours.A_p * r - s1)^2 + alpha(2) * (r' * stats_ours.A_p * stats_ours.A_p * r - s2)^2;
% % z = @(r) 1 - abs(opts_r.r_true' * r) ;
% f_val = [f(opts_r.r_true), f(R_ours * [0;0;1]), f(r_direct_RO), f(r_direct_spec)];
% g_val = [g(opts_r.r_true), g(R_ours * [0;0;1]), g(r_direct_RO), g(r_direct_spec)];
% f_val + g_val
% h_val = [h(opts_r.r_true), h(R_ours * [0;0;1]), h(r_direct_RO), h(r_direct_spec)]

tic;
[p_fischer, R_fischer, stats_fischer] = estimate_pose_fischer(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, alg_params);
time_fischer = toc;

% Rlm
tic;
[p_Rlm, R_Rlm, stats_Rlm] = estimate_R_lm( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );
time_Rlm = toc;

% ==== 结果存储 ====
result.p_lm    = p_lm;
result.R_lm    = R_lm;
result.p_elm   = p_elm;
result.R_elm   = R_elm;
result.p_ours  = p_ours;
result.R_ours  = R_ours;
result.r_direct_MM = r_direct_MM; % 存储直接估计的主轴
result.r_direct_SSL = r_direct_SSL;
result.r_direct_SDP = r_direct_SDP;
result.r_direct_RO  = r_direct_RO;
result.r_direct_spec = r_direct_spec;
result.p_fischer = p_fischer;
result.R_fischer = R_fischer;
result.p_Rlm   = p_Rlm;
result.R_Rlm   = R_Rlm;
result.p_init  = p_init;
result.theta_init = theta_init;
result.time_lm = time_lm;
result.time_elm = time_elm;
result.time_ours = time_ours;
result.time_fischer = time_fischer;
result.time_Rlm = time_Rlm;

% ==== 误差分析 ====
% LM
r_true = opts_r.r_true;
result.lm_pos_error = norm(p_lm - p_true);
result.lm_rot_error = norm(R_lm - R_true, 'fro');
r_hat_lm = R_lm * opts_r.u;
r_hat_lm = r_hat_lm / norm(r_hat_lm);
result.lm_r_error = 1 - abs(r_true' * r_hat_lm); % 主轴方向误差 (1-abs(inner))

% ELM
result.elm_pos_error = norm(p_elm - p_true);
result.elm_rot_error = norm(R_elm - R_true, 'fro');
r_hat_elm = R_elm * opts_r.u;
r_hat_elm = r_hat_elm / norm(r_hat_elm);
result.elm_r_error = 1 - abs(r_true' * r_hat_elm); % 主轴方向误差 (1-abs(inner))

% 所提算法
result.ours_pos_error = norm(p_ours - p_true);
result.ours_rot_error = norm(R_ours - R_true, 'fro');
r_hat_ours = R_ours * opts_r.u;
r_hat_ours = r_hat_ours / norm(r_hat_ours);
result.ours_r_error = 1 - abs(r_true' * r_hat_ours); % SSL 或内部默认估计的r误差
result.direct_r_error_MM = 1 - abs(r_true' * r_direct_MM); % 显式直接估计的r误差
result.direct_r_error_SSL = 1 - abs(r_true' * r_direct_SSL);
result.direct_r_error_SDP = 1 - abs(r_true' * r_direct_SDP);
result.direct_r_error_RO  = 1 - abs(r_true' * r_direct_RO);
result.direct_r_error_spec = 1 - abs(r_true' * r_direct_spec);

% Fischer
result.fischer_pos_error = norm(p_fischer - p_true);
result.fischer_rot_error = norm(R_fischer - R_true, 'fro');
r_hat_fischer = R_fischer * opts_r.u;
r_hat_fischer = r_hat_fischer / norm(r_hat_fischer);
result.fischer_r_error = 1 - abs(r_true' * r_hat_fischer); % 主轴方向误差 (1-abs(inner))

% Rlm
result.Rlm_pos_error = norm(p_Rlm - p_true);
result.Rlm_rot_error = norm(R_Rlm - R_true, 'fro');
r_hat_Rlm = R_Rlm * opts_r.u;
r_hat_Rlm = r_hat_Rlm / norm(r_hat_Rlm);
result.Rlm_r_error = 1 - abs(r_true' * r_hat_Rlm); 

end

