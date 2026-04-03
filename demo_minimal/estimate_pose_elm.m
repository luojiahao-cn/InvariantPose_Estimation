function [p_est, R_est, stats] = estimate_pose_elm(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p)
% ESTIMATE_POSE_ELM Optimizes magnetic field difference residuals
% Uses pairwise differences between sensors to eliminate sensor bias

num_sensors = size(b_total, 2);
pairs = nchoosek(1:num_sensors, 2);
num_pairs = size(pairs, 1);

delta_b_meas = zeros(3, num_pairs);
delta_d = zeros(3, num_pairs);

for idx = 1:num_pairs
    i = pairs(idx, 1);
    j = pairs(idx, 2);
    delta_b_meas(:, idx) = b_total(:, j) - b_total(:, i);
    delta_d(:, idx) = d_list(:, j) - d_list(:, i);
end

x0 = [theta_init; p_init];
lb = [-inf(3,1); lb_p(:)];
ub = [inf(3,1); ub_p(:)];

[x_opt, resnorm, residual, exitflag, output] = lsqnonlin(...
    @(x)elm_objective(x, m_pos, m_hat, m_norm, d_list, delta_d, delta_b_meas), ...
    x0, lb, ub, options);

theta_opt = x_opt(1:3);
p_est = x_opt(4:6);
R_est = MatrixExp3(VecToso3(theta_opt));

stats.resnorm = resnorm;
stats.exitflag = exitflag;
stats.output = output;

    function residuals = elm_objective(x, m_pos, m_hat, m_norm, d_list, delta_d, delta_b_meas)
        theta = x(1:3);
        p = x(4:6);
        R = MatrixExp3(VecToso3(theta));
        num_pairs = size(delta_d, 2);
        residuals = zeros(3*num_pairs, 1);

        for idx = 1:num_pairs
            i = pairs(idx, 1);
            j = pairs(idx, 2);
            p_i = p + R * d_list(:, i);
            p_j = p + R * d_list(:, j);

            B_i = zeros(3,1);
            for k = 1:size(m_pos,2)
                r_i = p_i - m_pos(:, k);
                [B_single, ~] = dipole_b_and_gradb(r_i, m_hat(:,k), m_norm(k));
                B_i = B_i + B_single;
            end

            B_j = zeros(3,1);
            for k = 1:size(m_pos,2)
                r_j = p_j - m_pos(:, k);
                [B_single, ~] = dipole_b_and_gradb(r_j, m_hat(:,k), m_norm(k));
                B_j = B_j + B_single;
            end

            b_i = R' * B_i;
            b_j = R' * B_j;
            delta_b_pred = b_j - b_i;
            residuals(3*(idx-1)+1:3*idx) = delta_b_pred - delta_b_meas(:, idx);
        end
    end
end
