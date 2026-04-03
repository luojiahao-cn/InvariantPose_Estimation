function [p_opt, R_opt, stats] = estimate_pose_lm(b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p)
% ESTIMATE_POSE_LM Levenberg-Marquardt algorithm for pose estimation
% Directly optimizes magnetic field residuals

x0 = [theta_init; p_init];
lb = [-inf(3,1); lb_p(:)];
ub = [inf(3,1); ub_p(:)];

[x_opt, resnorm, residual, exitflag, output] = lsqnonlin(...
    @(x)lm_objective(x, m_pos, m_hat, m_norm, d_list, b_total), ...
    x0, lb, ub, options);

omega_opt = x_opt(1:3);
p_opt = x_opt(4:6);
R_opt = MatrixExp3(VecToso3(omega_opt));

stats.resnorm = resnorm;
stats.residual = residual;
stats.exitflag = exitflag;
stats.output = output;

    function residuals = lm_objective(x, m_pos, m_hat, m_norm, d_list, b_list)
        omega = x(1:3);
        p = x(4:6);
        R = MatrixExp3(VecToso3(omega));
        num_sensors = size(d_list, 2);
        residuals = zeros(3*num_sensors, 1);

        for i = 1:num_sensors
            p_sensor_i = p + R * d_list(:, i);
            B_global = zeros(3,1);
            num_magnets = size(m_pos,2);
            for j = 1:num_magnets
                r = p_sensor_i - m_pos(:, j);
                [B_single, ~] = dipole_b_and_gradb(r, m_hat(:,j), m_norm(j));
                B_global = B_global + B_single;
            end
            b_pred = R' * B_global;
            b_meas = b_list(:, i);
            residuals(3*(i-1)+1:3*i) = b_pred - b_meas;
        end
    end
end
