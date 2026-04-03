function [p_est, R_est, stats] = estimate_pose_fischer( ...
    b_total, d_list, m_pos, m_hat, m_norm, ...
    theta_init, p_init, options, lb_p, ub_p, params)
% ESTIMATE_POSE_FISCHER Grid search with eigenvalue invariants
% Uses coarse-to-fine grid search for position, eigen decomposition for rotation

num_sensors = size(b_total, 2);

% Select 3 sensors for gradient estimation
if num_sensors > 6
    d_list = d_list(:, [1, 3, 7]);
    b_total = b_total(:, [1, 3, 7]);
elseif (num_sensors <= 6) && (num_sensors > 3)
    d_list = d_list(:, [1, 2, 4]);
    b_total = b_total(:, [1, 2, 4]);
else
    d_list = d_list(:, 1:3);
    b_total = b_total(:, 1:3);
end
num_sensors = 3;

% Estimate local gradient tensor using 3 sensors
d = d_list(:, 2) - d_list(:, 1);
e = d_list(:, 3) - d_list(:, 1);

Db12 = b_total(:, 2) - b_total(:, 1);
Db13 = b_total(:, 3) - b_total(:, 1);

dx = d(1); dy = d(2); dz = d(3);
ex = e(1); ey = e(2); ez = e(3);

% Linear system for gradient components
C = [ dx    dy    dz     0     0  ;
      0     dx    0      dy    dz ;
     -dz    0     dx    -dz    dy ;
      ex    ey    ez     0     0  ;
      0     ex    0      ey    ez ;
     -ez    0     ex    -ez    ey ];

h = [ Db12(1); Db12(2); Db12(3); Db13(1); Db13(2); Db13(3) ];

x = C \ h;

Gxx = x(1); Gxy = x(2); Gxz = x(3);
Gyy = x(4); Gyz = x(5);
Gzz = -Gxx - Gyy;

X_opt = [ Gxx Gxy Gxz ;
          Gxy Gyy Gyz ;
          Gxz Gyz Gzz ];

% Eigen decomposition in sensor frame
[V_s, D_s] = eig(X_opt);
[lambda_s, idx_s] = sort(diag(D_s), 'ascend');
V_s = V_s(:, idx_s);

if det(V_s) < 0
    V_s(:,1) = -V_s(:,1);
end
Rsv = V_s;
lambda_s = lambda_s(:);

% Measurement invariants
b_s = mean(b_total, 2);
norm_b_meas = norm(b_s);
norm_G_meas = norm(X_opt, 'fro');
invariants_meas = [norm_b_meas; norm_G_meas; lambda_s];

% Grid search parameters
epsilon       = 0.1;
grid_step     = 0.02;
num_iter      = 4;
shrink_factor = 1.4;

p_est = p_init;
cost_prev = inf;

cost_history = nan(num_iter, 1);
p_history    = nan(3, num_iter);

% Coarse-to-fine grid search for position
for it = 1:num_iter
    x_min = max(lb_p(1), p_est(1) - epsilon/2);
    x_max = min(ub_p(1), p_est(1) + epsilon/2);
    y_min = max(lb_p(2), p_est(2) - epsilon/2);
    y_max = min(ub_p(2), p_est(2) + epsilon/2);
    z_min = max(lb_p(3), p_est(3) - epsilon/2);
    z_max = min(ub_p(3), p_est(3) + epsilon/2);

    x_grid = x_min:grid_step:x_max;
    y_grid = y_min:grid_step:y_max;
    z_grid = z_min:grid_step:z_max;

    [Xg, Yg, Zg] = ndgrid(x_grid, y_grid, z_grid);
    P_grid = [Xg(:), Yg(:), Zg(:)];

    for k = 1:size(P_grid, 1)
        p = P_grid(k, :)';
        [b_w, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm);
        norm_b_mod = norm(b_w);
        norm_G_mod = norm(A_p, 'fro');
        lambda_mod = sort(eig(A_p), 'ascend');
        lambda_mod = lambda_mod(:);
        invariants_mod = [norm_b_mod; norm_G_mod; lambda_mod];
        cost = norm(invariants_mod - invariants_meas);

        if cost < cost_prev
            cost_prev = cost;
            p_est = p;
        end
    end

    cost_history(it) = cost_prev;
    p_history(:,it)  = p_est;
    epsilon   = epsilon   / shrink_factor;
    grid_step = grid_step / shrink_factor;
end

% Estimate rotation via eigen decomposition
[b_p, A_p_final] = calcFieldAndGradient(p_est, m_pos, m_hat, m_norm);

[V_w, D_w] = eig(A_p_final);
[lambda_w, idx_w] = sort(diag(D_w), 'ascend');
V_w = V_w(:, idx_w);

if det(V_w) < 0
    V_w(:,1) = -V_w(:,1);
end

Rwv = V_w;
lambda_w = lambda_w(:);
R_est = Rwv * Rsv';

stats = struct();
stats.X_opt        = X_opt;
stats.cost_history = cost_history;
stats.p_history    = p_history;
end
