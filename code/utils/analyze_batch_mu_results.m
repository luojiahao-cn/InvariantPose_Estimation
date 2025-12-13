% batch_results
function [pos_error_mean, pos_error_std, rot_error_mean, rot_error_std] = analyze_batch_mu_results(batch_results, num_test_points, num_trials)
    pos_error_mean = zeros(num_trials, 1);
    pos_error_std = zeros(num_trials, 1);
    rot_error_mean = zeros(num_trials, 1);
    rot_error_std = zeros(num_trials, 1);
    for i = 1:num_trials % 对每个trial
        pos_error = zeros(num_test_points, 1);
        rot_error = zeros(num_test_points, 1);
        % 统计每个W下的pos_error
        for j = 1:num_test_points % 对每个test_point
            pos_error(j) = batch_results(i).results(j).ours_pos_error;
            rot_error(j) = batch_results(i).results(j).ours_rot_error;
        end
        % 统计每个W下的pos_error和rot_error
        pos_error_mean(i) = mean(pos_error);
        pos_error_std(i) = std(pos_error);
        rot_error_mean(i) = mean(rot_error);
        rot_error_std(i) = std(rot_error);
    end

    % sprintf('pos_error_mean = %.4f, pos_error_std = %.4f, rot_error_mean = %.4f, rot_error_std = %.4f', pos_error_mean, pos_error_std, rot_error_mean, rot_error_std);