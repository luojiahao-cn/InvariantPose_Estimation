function [X_opt, x_opt, D_delta, B_delta] = lc_grad_tensor_estimator(b_total, d_list)
% LC_GRAD_TENSOR_ESTIMATOR Estimate local gradient tensor from sensor measurements
% Uses pairwise differences between sensors to construct a linear system

    num_sensors = size(b_total, 2);
    pairs = nchoosek(1:num_sensors, 2);
    D_delta = zeros(3, size(pairs,1));
    B_delta = zeros(3, size(pairs,1));

    for idx = 1:size(pairs,1)
        i = pairs(idx, 1);
        j = pairs(idx, 2);
        D_delta(:, idx) = d_list(:, j) - d_list(:, i);
        B_delta(:, idx) = b_total(:, j) - b_total(:, i);
    end

    % Selection matrix S for gradient tensor components
    S = [1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0;
        0,1,0,0,0; 0,0,0,1,0; 0,0,0,0,1;
        0,0,1,0,0; 0,0,0,0,1; -1,0,0,-1,0];

    C_matrix = kron(D_delta', eye(3)) * S;

    if rank(C_matrix) < 5
        warning('Gradient system matrix C is rank deficient.');
    end

    h_vector = B_delta(:);
    x_opt = pinv(C_matrix) * h_vector;
    X_opt = reshape(S * x_opt, 3, 3);
end
