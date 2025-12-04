function plot_beta_results(batch_results)
    % 提取所有测试点的结果
    num_points = length(batch_results);
    
    % batch_results：test_point, results
    % results：beta_config, R_true, R_beta_history, R_init_est
    % R_beta_history：R, k, kmax, deltabd, beta, delta_history
    % R_init_est：R_init_est1, R_init_est2, eR_init_est1, eR_init_est2
    % R_true：R_true
    % test_point：p_true, theta_true
    % R_beta_history：R, k, kmax, deltabd, beta, delta_history
    % R_init_est：R_init_est1, R_init_est2, eR_init_est1, eR_init_est2

    %% 1. 提取所有的beta和对应的eR
    beta_vec = [];
    eR_vec = [];
    eR_init_vec = [];
    for i = 1:num_points
        for j = 1:batch_results(i).results{1}.num_beta
            beta_vec = [beta_vec, batch_results(i).results{1}.R_beta_history{j}.beta];
            eR_vec = [eR_vec, batch_results(i).results{1}.R_beta_history{j}.eR];
            iter_step(j) = batch_results(i).results{1}.R_beta_history{j}.k;
            delta_history{j} = batch_results(i).results{1}.R_beta_history{j}.delta_history;
            eR_init_vec = [eR_init_vec, batch_results(i).results{1}.R_init_est.eR_init_est1, batch_results(i).results{1}.R_init_est.eR_init_est2];
            eR_history{j} = batch_results(i).results{1}.R_beta_history{j}.eR_history;
        end
    end
    %% 2. 绘制beta和eR的关系图
    figure
    plot(beta_vec, eR_vec, 'b-o', 'LineWidth', 2);
    for j = 1:batch_results(i).results{1}.num_beta
        plot(delta_history{j}, 'LineWidth', 2);
        hold on;
    end
    % xlabel('Iteration');
    % ylabel('delta');
    % title('delta vs iteration');
    % grid on;

    figure
    plot(beta_vec, iter_step, 'b-o', 'LineWidth', 2);

    figure
    plot(beta_vec, eR_vec, 'r-o', 'LineWidth', 2);
end