function [p_est, R_est, stats] = estimate_pose_union(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p)
% ESTIMATE_POSE_UNION 融合LM和ELM方法的姿态估计算法
%
% 输入参数：
%   b_total    - 3×N传感器测量矩阵（局部坐标系），单位：T
%   d_list     - 3×N位移矩阵，传感器在参考坐标系中的偏移，单位：m
%   m_pos      - 3×K磁铁位置矩阵，单位：m
%   m_hat      - 3×K磁化方向单位向量（归一化）
%   m_norm     - 1×K磁矩幅值向量，单位：A·m²
%   theta_init - 3×1初始旋转向量估计 [rad]
%   p_init     - 3×1初始位置估计 [m]
%
% 输出参数：
%   p_est    - 估计的位置向量
%   R_est    - 估计的旋转矩阵
%   stats    - 包含优化统计信息的结构体

%% ===== 算法参数配置 =====
params.dynamic_weight = false;    % 启用动态权重调整
params.w_lm = 0.5;               % LM残差初始权重
params.w_elm = 0.5;              % ELM残差初始权重
params.adaptive_decay = 0.95;    % 权重自适应衰减率
params.min_weight = 0.1;         % 权重下限
params.max_weight = 0.9;         % 权重上限
params.gradient_threshold = 1e-3; % 梯度调整阈值

num_sensors = size(b_total, 2);

%% ===== 准备传感器对信息 =====
pairs = nchoosek(1:num_sensors, 2);
num_pairs = size(pairs, 1);

% 计算传感器对的磁场差和位移差
delta_b_meas = zeros(3, num_pairs);
delta_d = zeros(3, num_pairs);

for idx = 1:num_pairs
    i = pairs(idx, 1);
    j = pairs(idx, 2);
    delta_b_meas(:, idx) = b_total(:, j) - b_total(:, i);
    delta_d(:, idx) = d_list(:, j) - d_list(:, i);
end

%% ===== 优化设置 =====
x0 = [theta_init; p_init];
lb = [-inf(3,1); lb_p(:)];
ub = [inf(3,1); ub_p(:)];

% 初始化动态权重
current_w_lm = params.w_lm;
current_w_elm = params.w_elm;
iteration_count = 0;

%% ===== 运行优化 =====
[x_opt, ~, ~, exitflag, output] = lsqnonlin(...
    @(x)fused_objective(x), ...
    x0, lb, ub, options);

% 提取结果
theta_opt = x_opt(1:3);
p_est = x_opt(4:6);
R_est = MatrixExp3(VecToso3(theta_opt));

% 保存统计信息
stats.exitflag = exitflag;
stats.output = output;
stats.final_w_lm = current_w_lm;
stats.final_w_elm = current_w_elm;

%% ===== 融合目标函数 =====
function residuals = fused_objective(x)
    % 提取参数
    theta = x(1:3);
    p = x(4:6);
    R = MatrixExp3(VecToso3(theta));
    
    % ===== 1. LM残差部分（单个传感器） =====
    lm_residuals = zeros(3*num_sensors, 1);
    for i = 1:num_sensors
        % 传感器全局位置
        p_sensor_i = p + R * d_list(:, i);
        
        % 计算磁场（全局坐标系）
        B_global = zeros(3,1);
        for j = 1:size(m_pos, 2)
            r = p_sensor_i - m_pos(:, j);
            [B, ~] = dipole_b_and_gradb(r, m_hat(:, j), m_norm(j));
            B_global = B_global + B;
        end
        
        % 转换到传感器坐标系
        b_pred = R' * B_global;
        b_meas = b_total(:, i);
        
        % 残差（应用权重）
        lm_residuals(3*(i-1)+1:3*i) = sqrt(current_w_lm) * (b_pred - b_meas);
    end
    
    % ===== 2. ELM残差部分（传感器对） =====
    elm_residuals = zeros(3*num_pairs, 1);
    for idx = 1:num_pairs
        i = pairs(idx, 1);
        j = pairs(idx, 2);
        
        % 计算传感器位置
        p_i = p + R * d_list(:, i);
        p_j = p + R * d_list(:, j);
        
        % 计算传感器i处的磁场
        B_i_global = zeros(3,1);
        for k = 1:size(m_pos, 2)
            r_i = p_i - m_pos(:, k);
            [B, ~] = dipole_b_and_gradb(r_i, m_hat(:, k), m_norm(k));
            B_i_global = B_i_global + B;
        end
        b_i_pred = R' * B_i_global;
        
        % 计算传感器j处的磁场
        B_j_global = zeros(3,1);
        for k = 1:size(m_pos, 2)
            r_j = p_j - m_pos(:, k);
            [B, ~] = dipole_b_and_gradb(r_j, m_hat(:, k), m_norm(k));
            B_j_global = B_j_global + B;
        end
        b_j_pred = R' * B_j_global;
        
        % 预测磁场差
        delta_b_pred = b_j_pred - b_i_pred;
        
        % 残差（应用权重）
        elm_residuals(3*(idx-1)+1:3*idx) = sqrt(current_w_elm) * (delta_b_pred - delta_b_meas(:, idx));
    end
    
    % ===== 3. 融合残差 =====
    residuals = [lm_residuals; elm_residuals];
end

%% ===== 动态权重调整函数 =====
function stop = optim_output(x, optimValues, state)
    stop = false;
    
    if strcmp(state, 'iter')
        iteration_count = iteration_count + 1;
        
        % 每3次迭代调整一次权重
        if mod(iteration_count, 3) == 0 && params.dynamic_weight
            % 根据当前梯度信息调整权重
            current_gradient = norm(optimValues.gradient);
            
            if current_gradient > params.gradient_threshold
                % 梯度较大时增加LM权重（强调精确建模）
                current_w_lm = min(params.max_weight, current_w_lm * 1.1);
                current_w_elm = max(params.min_weight, current_w_elm * 0.9);
            else
                % 梯度较小时增加ELM权重（强调稳定性）
                current_w_elm = min(params.max_weight, current_w_elm * 1.1);
                current_w_lm = max(params.min_weight, current_w_lm * 0.9);
            end
            
            % 归一化权重
            total_w = current_w_lm + current_w_elm;
            current_w_lm = current_w_lm / total_w;
            current_w_elm = current_w_elm / total_w;
            
            fprintf('[UNION] 迭代 %d: 权重调整 w_lm=%.3f, w_elm=%.3f, 梯度=%.4e\n', ...
                iteration_count, current_w_lm, current_w_elm, current_gradient);
        end
    end
end

end

