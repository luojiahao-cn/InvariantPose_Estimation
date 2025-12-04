%% Estimate R
function R_struct = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init, mu, beta)
    %% Iterative method:
    LA = min(eig(A_p)) - 1e-3;
    LX = min(eig(X_opt)) - 1e-3;
    Abar = A_p - LA * eye(3);
    Xbar = X_opt - LX * eye(3);
    L = 4 * norm(Abar) * norm(Xbar);
    beta = beta * L;

    M = @(R) 2 * Abar * R * Xbar + mu * B_bar * b_bar';

    Mbar = @(R) M(R) + beta * R;
    % f = @(R) trace(R' * A_p * R * X_opt + mu * R' * B_bar * b_bar' + beta * R);
    % fbar = @(R) trace(R' * At * R * Xt + mu * R' * B_bar * b_bar' + beta * R);
    % skw = @(R) 0.5 * (R - R');
    % g = @(R) 2 * skw(R' * (2 * A_p * R * X_opt + mu * B_bar * b_bar' + beta * eye(3)));
    % gbar = @(R) 2 * skw(R' * (2 * At * R * Xt + mu * B_bar * b_bar' + beta * eye(3)));
    %% Initial condition
    delta = 1e6; k = 0; kmax = 1000;
    R = R_init;
    R_iter_history = {};
    delta_history = [];
    while k < kmax && delta > 1e-5
        [U, ~, V] = svd(Mbar(R));
        R_opt = U*diag([1,1,det(U*V')])*V';
        delta = norm(R_opt - R, 'fro');
        R_iter_history{end+1} = R_opt;
        delta_history(end+1) = delta;
        R = R_opt;
        k = k + 1;
    end
    if k == kmax
        warning('迭代未收敛，可能需要调整参数');
    end
    R_struct.R = R;
    R_struct.R_iter_history = R_iter_history;
    R_struct.delta_history = delta_history;
    end