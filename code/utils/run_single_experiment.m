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

% ==== 初始值猜测 ====

% 初始化p_init为在以p_true为球心，半径为p_uncertainty，p_init(3)>=0的球内随机采样
% 球面内均匀采样，生成半径为 <= p_uncertainty 的向量
u = randn(3,1);       % 随机方向
u(3) = abs(u(3));
u = u / norm(u);      % 单位向量
r = (rand())^(1/3) * p_uncertainty;  % 半径分布
p_init = p_true + r * u;
% p_init = [0; 0; 0];

% theta_init = theta_true + r_uncertainty * init_error;
% theta_init的随机初始化，基于真实theta_true和r_uncertainty，采用旋转向量扰动
u_theta = randn(3,1);
u_theta = u_theta / norm(u_theta);
r_theta = (rand())^(1/3) * r_uncertainty;
theta_init = theta_true + r_theta * u_theta;
% theta_init = [0; 0; 0];
% 这里的theta_init为旋转向量（轴角），后续调用算法时将通过MatrixExp3(VecToso3(theta_init))转成旋转矩阵

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
alg_params = struct('R_true', R_true, 'p_true', p_true, 'mu', mu, 'beta', beta);
[p_ours, R_ours, stats_ours] = estimate_pose_ours( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, alg_params);
time_ours = toc;

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
result.lm_pos_error = norm(p_lm - p_true);
result.lm_rot_error = norm(R_lm - R_true, 'fro');

% ELM
result.elm_pos_error = norm(p_elm - p_true);
result.elm_rot_error = norm(R_elm - R_true, 'fro');

% 所提算法
result.ours_pos_error = norm(p_ours - p_true);
result.ours_rot_error = norm(R_ours - R_true, 'fro');

% Fischer
result.fischer_pos_error = norm(p_fischer - p_true);
result.fischer_rot_error = norm(R_fischer - R_true, 'fro');

% Rlm
result.Rlm_pos_error = norm(p_Rlm - p_true);
result.Rlm_rot_error = norm(R_Rlm - R_true, 'fro');

end

