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
init_error = -1 + 2 * rand(3,1); % [-1, 1]
p_init = p_true + p_uncertainty * init_error;
theta_init = theta_true + r_uncertainty * init_error;

% ==== 算法调用 ====
% LM
[p_lm, R_lm, stats_lm] = estimate_pose_lm( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );

% ELM
[p_elm, R_elm, stats_elm] = estimate_pose_elm( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p );

% 所提算法
alg_params = struct('R_true', R_true, 'p_true', p_true, 'mu', mu, 'beta', beta);
[p_ours, R_ours, stats_ours] = estimate_pose_ours( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, alg_params);

% Rlm
[p_Rlm, R_Rlm, stats_Rlm] = estimate_R_lm( ...
    b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_ours, options );

% ==== 结果存储 ====
result.p_lm    = p_lm;
result.R_lm    = R_lm;
result.p_elm   = p_elm;
result.R_elm   = R_elm;
result.p_ours  = p_ours;
result.R_ours  = R_ours;
result.p_Rlm   = p_Rlm;
result.R_Rlm   = R_Rlm;

% ==== 误差分析 ====
% LM
result.lm_pos_error = norm(p_lm - p_true);
result.lm_rot_error = norm(R_lm - R_true, 'fro');
result.lm_field_error = calculate_field_errors(p_lm, R_lm, b_total, d_list, m_pos, m_hat, m_norm);

% ELM
result.elm_pos_error = norm(p_elm - p_true);
result.elm_rot_error = norm(R_elm - R_true, 'fro');
result.elm_field_error = calculate_field_errors(p_elm, R_elm, b_total, d_list, m_pos, m_hat, m_norm);

% 所提算法
result.ours_pos_error = norm(p_ours - p_true);
result.ours_rot_error = norm(R_ours - R_true, 'fro');
result.ours_field_error = calculate_field_errors(p_ours, R_ours, b_total, d_list, m_pos, m_hat, m_norm);

% Rlm
result.Rlm_pos_error = norm(p_Rlm - p_true);
result.Rlm_rot_error = norm(R_Rlm - R_true, 'fro');
result.Rlm_field_error = calculate_field_errors(p_Rlm, R_Rlm, b_total, d_list, m_pos, m_hat, m_norm);

end

