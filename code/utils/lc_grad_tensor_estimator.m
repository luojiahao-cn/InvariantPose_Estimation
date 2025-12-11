function [X_opt, x_opt] = lc_grad_tensor_estimator(b_total, d_list)
    num_sensors = size(b_total, 2);
    pairs = nchoosek(1:num_sensors, 2);
    D_matrix = zeros(3, size(pairs,1));
    B_matrix = zeros(3, size(pairs,1));

    for idx = 1:size(pairs,1)
        i = pairs(idx, 1);
        j = pairs(idx, 2);
        D_matrix(:, idx) = d_list(:, j) - d_list(:, i);
        B_matrix(:, idx) = b_total(:, j) - b_total(:, i);
    end
    %% 估计局部梯度张量
    % 构建选择矩阵S
    S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
        0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
        0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

    % 构建完整约束矩阵C
    C_matrix = kron(D_matrix', eye(3)) * S;

    if rank(C_matrix) < 5
        warning('lc_grad_tensor_estimator:rankDeficient', ...
            'Gradient system matrix C is rank deficient. Sensor geometry may be ill conditioned.');
    end

    h_vector = B_matrix(:);

    % 求解最小二乘问题
    x_opt = pinv(C_matrix) * h_vector;
    X_opt = reshape(S * x_opt, 3, 3);  % 估计梯度（传感器坐标系）

    %% 阶段检查 对应公式2
    % [~, A_p_true] = calcFieldAndGradient(params.p_true, m_pos, m_hat, m_norm);
    % R_true = params.R_true;
    % X_true = R_true'*A_p_true*R_true;
    % norm(R_true'*b_p*ones(1,num_sensors)+R_true'*A_p*R_true*d_list - b_total, 'fro')
    % norm(R_true'*A_p*R_true*D_matrix - B_matrix, 'fro')

end