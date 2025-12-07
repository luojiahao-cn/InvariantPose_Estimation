%% Estimate R
function R_struct = estimateR_iter(b_bar, B_bar, A_p, X_opt, R_init, mu, beta, R_true)
    %% Iterative method:
    LA = min(eig(A_p)) - eps;
    LX = min(eig(X_opt)) - eps;
    Abar = A_p - LA * eye(3);
    Xbar = X_opt - LX * eye(3);

    L = 4 * norm(Abar) * norm(Xbar);
    % beta = beta * L;
    mu = mu * L;
    beta = beta * L / mu;

    M = @(R) 2 * Abar * R * Xbar + mu * B_bar * b_bar';

    Mbar = @(R) M(R) + beta * R;
    % f = @(R) trace(R' * A_p * R * X_opt + mu * R' * B_bar * b_bar' + beta * R);
    fbar = @(R) trace(R' * Abar * R * Xbar + mu * R' * B_bar * b_bar' + beta * R);
    % skw = @(R) 0.5 * (R - R');
    % g = @(R) 2 * skw(R' * (2 * A_p * R * X_opt + mu * B_bar * b_bar' + beta * eye(3)));
    % gbar = @(R) 2 * skw(R' * (2 * At * R * Xt + mu * B_bar * b_bar' + beta * eye(3)));
    %% Initial condition
    delta = 1e6; k = 0; kmax = 1e4;
    R = R_init;
    R_iter_history{1} = R_init;
    delta_history = [];
    eR_history = norm(R_true - R_init, 'fro');
    f_history = fbar(R_init);
    deltabd = 0;
    while k < kmax && delta > deltabd
        [U, ~, V] = svd(Mbar(R));
        R_opt = U*diag([1,1,det(U*V')])*V';
        delta = norm(R_opt - R, 'fro');
        R = R_opt;
        k = k + 1;

        R_iter_history{end+1} = R_opt;
        eR_history(end+1) = norm(R_true - R_opt, 'fro');
        delta_history(end+1) = delta;
        f_history(end+1) = fbar(R_opt);
    end
    if k == kmax
        warning('迭代未收敛，可能需要调整参数');
    end

    R_struct.R = R;
    R_struct.R_iter_history = R_iter_history;
    R_struct.eR_history = eR_history;
    R_struct.L = L;
    R_struct.k = k;
    R_struct.kmax = kmax;
    R_struct.deltabd = 1e-5;
    R_struct.beta = beta;
    R_struct.delta_history = delta_history;
    R_struct.f_history = f_history;
    end