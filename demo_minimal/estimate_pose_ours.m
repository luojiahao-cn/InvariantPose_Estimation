function [p_est, R_est, stats] = estimate_pose_ours(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, params)
% PROPOSED_METHOD_POSE_ESTIMATION Two-stage method:
% Stage 1: Position estimation via eigenvalue matching
% Stage 2: Rotation estimation via iterative Procrustes

R_init = MatrixExp3(VecToso3(theta_init));
num_sensors = size(b_total, 2);

[X_opt, ~, D_delta, B_delta] = lc_grad_tensor_estimator(b_total, d_list);

%% Stage #1: Estimate position using eigenvalue matching
[Q, ~, ~] = qr(d_list');
r = rank(d_list);
Q_bar = Q(:, r+1:end);
b_bar = b_total * Q_bar;

% Grid search for initial position
best_p = p_init;
best_cost = inf;

% Coarse grid search
grid_size = 0.02;
for px = lb_p(1):grid_size:ub_p(1)
    for py = lb_p(2):grid_size:ub_p(2)
        for pz = lb_p(3):grid_size:ub_p(3)
            p = [px; py; pz];
            cost = obj_cost(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt);
            if cost < best_cost
                best_cost = cost;
                best_p = p;
            end
        end
    end
end

% Refine with fminunc
p_est = best_p;
options_fminunc = optimoptions('fminunc', 'Display', 'off', 'MaxIterations', 100, 'StepTolerance', 1e-8);
p_est = fminunc(@(p) obj_cost(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt), p_est, options_fminunc);

% Additional refinement with bounds using fmincon
p_est = fmincon(@(p) obj_cost(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X_opt), p_est, [], [], [], [], lb_p, ub_p);

%% Stage #2: Estimate rotation
[b_p, A_p] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);
B_matrix = b_p * ones(1, num_sensors);
B_bar = B_matrix * Q_bar;
[R_init_est1, R_init_est2] = estimateR(b_bar, B_bar, A_p, X_opt, D_delta, B_delta);

mu = params.optimization.mu;
beta = params.optimization.beta;
R_true = params.R_true;
R_PPI = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init, mu, beta, R_true);

R_est = R_PPI.R;
stats.X_opt = X_opt;
stats.R_est_init1 = R_init_est1;
stats.R_est_init2 = R_init_est2;
end

function cost = obj_cost(p, m_pos, m_hat, m_norm, num_sensors, b_bar, Q_bar, X)
    [b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);
    B_bar = b_p * ones(1, num_sensors) * Q_bar;

    % Eigenvalue matching cost (use real part to handle numerical issues)
    eig_B = real(eig(B_bar * B_bar'));
    eig_b = real(eig(b_bar * b_bar'));
    eig_A = real(eig(A_p));
    eig_X = real(eig(X));

    term1 = sum((sort(eig_B, 'descend') - sort(eig_b, 'descend')).^2);
    term5 = sum((sort(eig_A, 'descend') - sort(eig_X, 'descend')).^2);

    cost = term1 + term5;
end
