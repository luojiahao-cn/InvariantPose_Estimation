%% run_batch_demo - Batch Random Point Test Main Script
% This script generates multiple random test points, runs all algorithms on each,
% and compares the results.

clear all; close all; clc;
addpath helper

%% 1. Setup Parameters
fprintf('=== Initializing Parameters ===\n');
params = get_params();
m_pos = params.magnet.m_pos;
m_hat = params.magnet.m_hat;
m_norm = params.magnet.m_norm;
d_list = params.sensor.d_list;
lb_p = params.workspace.center - params.workspace.radius * ones(3,1);
ub_p = params.workspace.center + params.workspace.radius * ones(3,1);
options = params.optimization.options;
mu = params.optimization.mu;
beta = params.optimization.beta;
W = params.optimization.W;

%% 2. Test Configuration
num_test_points = 10;
fprintf('Generating %d random test points...\n', num_test_points);

%% 3. Pre-allocate Results Storage
results = struct('ours', {}, 'lm', {}, 'elm', {}, 'fischer', {});

%% 4. Run Batch Test
for test_idx = 1:num_test_points
    fprintf('\n--- Test Point %d/%d ---\n', test_idx, num_test_points);

    % Generate random test point
    test_points = generate_test_points(params, 1);
    p_true = test_points.p_true;
    theta_true = test_points.theta_true;

    % Set ground truth
    params.ground_truth.p_true = p_true;
    params.ground_truth.theta_true = theta_true;
    R_true = MatrixExp3(VecToso3(theta_true));
    params.R_true = R_true;

    % Generate magnetic data
    [b_total, sensor_positions] = generate_magnetic_data(params);

    % Initial estimate (identity rotation and zero position)
    theta_init = zeros(3,1);
    p_init = zeros(3,1);

    % Run Proposed Algorithm
    try
        [p_est_ours, R_est_ours, ~] = estimate_pose_ours(...
            b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, params);
        error_ours.pos(test_idx) = norm(p_est_ours - p_true);
        error_ours.rot(test_idx) = norm(R_est_ours - R_true, 'fro') / sqrt(2);
    catch ME
        error_ours.pos(test_idx) = NaN;
        error_ours.rot(test_idx) = NaN;
        fprintf('Ours error: %s\n', ME.message);
    end

    % Run LM Algorithm
    try
        [p_est_lm, R_est_lm, ~] = estimate_pose_lm(...
            b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p);
        error_lm.pos(test_idx) = norm(p_est_lm - p_true);
        error_lm.rot(test_idx) = norm(R_est_lm - R_true, 'fro') / sqrt(2);
    catch ME
        error_lm.pos(test_idx) = NaN;
        error_lm.rot(test_idx) = NaN;
        fprintf('LM error: %s\n', ME.message);
    end

    % Run ELM Algorithm
    try
        [p_est_elm, R_est_elm, ~] = estimate_pose_elm(...
            b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p);
        error_elm.pos(test_idx) = norm(p_est_elm - p_true);
        error_elm.rot(test_idx) = norm(R_est_elm - R_true, 'fro') / sqrt(2);
    catch ME
        error_elm.pos(test_idx) = NaN;
        error_elm.rot(test_idx) = NaN;
        fprintf('ELM error: %s\n', ME.message);
    end

    % Run Fischer Algorithm
    try
        [p_est_fischer, R_est_fischer, ~] = estimate_pose_fischer(...
            b_total, d_list, m_pos, m_hat, m_norm, theta_init, p_init, options, lb_p, ub_p, params);
        error_fischer.pos(test_idx) = norm(p_est_fischer - p_true);
        error_fischer.rot(test_idx) = norm(R_est_fischer - R_true, 'fro') / sqrt(2);
    catch ME
        error_fischer.pos(test_idx) = NaN;
        error_fischer.rot(test_idx) = NaN;
        fprintf('Fischer error: %s\n', ME.message);
    end

    % Print current test point errors
    if ~isnan(error_ours.pos(test_idx))
        fprintf('Ours:     pos=%.4f mm, rot=%.4f rad\n', error_ours.pos(test_idx)*1e3, error_ours.rot(test_idx));
    else
        fprintf('Ours:     pos=NaN, rot=NaN\n');
    end
    if ~isnan(error_lm.pos(test_idx))
        fprintf('LM:       pos=%.4f mm, rot=%.4f rad\n', error_lm.pos(test_idx)*1e3, error_lm.rot(test_idx));
    else
        fprintf('LM:       pos=NaN, rot=NaN\n');
    end
    if ~isnan(error_elm.pos(test_idx))
        fprintf('ELM:      pos=%.4f mm, rot=%.4f rad\n', error_elm.pos(test_idx)*1e3, error_elm.rot(test_idx));
    else
        fprintf('ELM:      pos=NaN, rot=NaN\n');
    end
    if ~isnan(error_fischer.pos(test_idx))
        fprintf('Fischer:  pos=%.4f mm, rot=%.4f rad\n', error_fischer.pos(test_idx)*1e3, error_fischer.rot(test_idx));
    else
        fprintf('Fischer:  pos=NaN, rot=NaN\n');
    end
end

%% 5. Statistics
fprintf('\n========== Statistics ==========\n');
fprintf('Algorithm | Position Error (mm)       | Rotation Error (rad)\n');
fprintf('          | Mean   Std    Max         | Mean   Std    Max\n');
fprintf('----------------------------------------------------------------\n');

print_stats('Ours', error_ours);
print_stats('LM', error_lm);
print_stats('ELM', error_elm);
print_stats('Fischer', error_fischer);

fprintf('\nTest Complete!\n');

%% Helper Function
function print_stats(name, error_struct)
    pos = error_struct.pos;
    rot = error_struct.rot;
    pos_mean = nanmean(pos) * 1e3;
    pos_std = nanstd(pos) * 1e3;
    pos_max = nanmax(pos) * 1e3;
    rot_mean = nanmean(rot);
    rot_std = nanstd(rot);
    rot_max = nanmax(rot);
    fprintf('%-8s | %6.2f  %6.2f  %6.2f  | %6.4f  %6.4f  %6.4f\n', name, pos_mean, pos_std, pos_max, rot_mean, rot_std, rot_max);
end
