function analyze_beta_results(batch_results)
    % 提取所有测试点的结果
    num_points = length(batch_results);
    num_trials_per_point = length(batch_results(1).results);
    % batch_results：test_point, results
    % results：beta_config, R_true, R_beta_history, R_init_est
    % R_beta_history：R, k, kmax, deltabd, beta, delta_history
    % R_init_est：R_init_est1, R_init_est2, eR_init_est1, eR_init_est2
    % R_true：R_true
    % test_point：p_true, theta_true
    % R_beta_history：R, k, kmax, deltabd, beta, delta_history
    % R_init_est：R_init_est1, R_init_est2, eR_init_est1, eR_init_est2

    %% 1. 提取所有的beta和对应的eR
    k = 1;
    beta_vec = [];
    eR_vec = [];
    iter_step = [];
    delta_history = [];
    eR_init_vec = [];
    eR_history = [];
    f_history = [];
    for j = 1:batch_results(1).results{1}.num_beta
        for i = 1:num_points
            beta_vec(i) = batch_results(i).results{k}.R_beta_history{j}.beta;
            eR_vec(i) = batch_results(i).results{k}.R_beta_history{j}.eR;
            iter_step(j) = batch_results(i).results{k}.R_beta_history{j}.k;
            delta_history{j} = batch_results(i).results{k}.R_beta_history{j}.delta_history;
            eR_init_vec(i) = batch_results(i).results{k}.R_init_est.eR_init_est;
            eR_history{j} = batch_results(i).results{k}.R_beta_history{j}.eR_history;
            f_history{j} = batch_results(i).results{k}.R_beta_history{j}.f_history;
        end
        eR_vec_mean(j) = mean(eR_vec);
        eR_vec_std(j) = std(eR_vec);
        iter_step_mean(j) = mean(iter_step);
        iter_step_std(j) = std(iter_step);
    end
end