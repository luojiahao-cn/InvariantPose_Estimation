function analyze_batch_results(batch_results)
% ANALYZE_BATCH_RESULTS 分析批量实验结果
% 输入：
%   batch_results - 由 run_batch_experiments 返回的批量结果结构体数组

num_points = length(batch_results);

% 收集所有测试点的统计信息
all_lm_pos_means = zeros(num_points, 1);
all_lm_rot_means = zeros(num_points, 1);
all_elm_pos_means = zeros(num_points, 1);
all_elm_rot_means = zeros(num_points, 1);
all_ours_pos_means = zeros(num_points, 1);
all_ours_rot_means = zeros(num_points, 1);
all_Rlm_pos_means = zeros(num_points, 1);
all_Rlm_rot_means = zeros(num_points, 1);

for i = 1:num_points
    all_lm_pos_means(i) = batch_results(i).summary.lm.pos_mean;
    all_lm_rot_means(i) = batch_results(i).summary.lm.rot_mean;
    all_elm_pos_means(i) = batch_results(i).summary.elm.pos_mean;
    all_elm_rot_means(i) = batch_results(i).summary.elm.rot_mean;
    all_ours_pos_means(i) = batch_results(i).summary.ours.pos_mean;
    all_ours_rot_means(i) = batch_results(i).summary.ours.rot_mean;
    all_Rlm_pos_means(i) = batch_results(i).summary.Rlm.pos_mean;
    all_Rlm_rot_means(i) = batch_results(i).summary.Rlm.rot_mean;
end

% 打印总体统计
fprintf('\n========== 批量实验结果统计 ==========\n');
fprintf('测试点总数: %d\n\n', num_points);

fprintf('位置误差统计 (m):\n');
fprintf('  LM:   均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_lm_pos_means), std(all_lm_pos_means), max(all_lm_pos_means));
fprintf('  ELM:  均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_elm_pos_means), std(all_elm_pos_means), max(all_elm_pos_means));
fprintf('  Ours: 均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_ours_pos_means), std(all_ours_pos_means), max(all_ours_pos_means));
fprintf('  Rlm:  均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_Rlm_pos_means), std(all_Rlm_pos_means), max(all_Rlm_pos_means));

fprintf('\n旋转误差统计:\n');
fprintf('  LM:   均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_lm_rot_means), std(all_lm_rot_means), max(all_lm_rot_means));
fprintf('  ELM:  均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_elm_rot_means), std(all_elm_rot_means), max(all_elm_rot_means));
fprintf('  Ours: 均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_ours_rot_means), std(all_ours_rot_means), max(all_ours_rot_means));
fprintf('  Rlm:  均值=%.6f, 标准差=%.6f, 最大值=%.6f\n', ...
    mean(all_Rlm_rot_means), std(all_Rlm_rot_means), max(all_Rlm_rot_means));

% 找出表现最好和最差的测试点
[~, best_lm_idx] = min(all_lm_pos_means);
[~, worst_lm_idx] = max(all_lm_pos_means);
[~, best_ours_idx] = min(all_ours_pos_means);
[~, worst_ours_idx] = max(all_ours_pos_means);

fprintf('\n最佳/最差测试点:\n');
p_best_lm = batch_results(best_lm_idx).test_point.p_true;
p_worst_lm = batch_results(worst_lm_idx).test_point.p_true;
p_best_ours = batch_results(best_ours_idx).test_point.p_true;
p_worst_ours = batch_results(worst_ours_idx).test_point.p_true;

fprintf('  LM最佳: 点%d, p=[%.4f,%.4f,%.4f], 误差=%.6f\n', ...
    best_lm_idx, p_best_lm(1), p_best_lm(2), p_best_lm(3), ...
    batch_results(best_lm_idx).summary.lm.pos_mean);
fprintf('  LM最差: 点%d, p=[%.4f,%.4f,%.4f], 误差=%.6f\n', ...
    worst_lm_idx, p_worst_lm(1), p_worst_lm(2), p_worst_lm(3), ...
    batch_results(worst_lm_idx).summary.lm.pos_mean);
fprintf('  Ours最佳: 点%d, p=[%.4f,%.4f,%.4f], 误差=%.6f\n', ...
    best_ours_idx, p_best_ours(1), p_best_ours(2), p_best_ours(3), ...
    batch_results(best_ours_idx).summary.ours.pos_mean);
fprintf('  Ours最差: 点%d, p=[%.4f,%.4f,%.4f], 误差=%.6f\n', ...
    worst_ours_idx, p_worst_ours(1), p_worst_ours(2), p_worst_ours(3), ...
    batch_results(worst_ours_idx).summary.ours.pos_mean);

end

