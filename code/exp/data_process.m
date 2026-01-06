clc,clear,close all
%% Environment setup
this_file = mfilename('fullpath');
exp_dir = fileparts(this_file);
repo_root = fileparts(exp_dir);
addpath(fullfile(repo_root, 'utils'));
addpath(fullfile(repo_root, 'Functions'));
addpath(fullfile(repo_root, 'tools'));
%%
load('scan_records.mat', 'scan_records');
data = scan_records;

R.R1 = [[1;0;0], [0;0;-1], [0;1;0]];
R.R2 = [[1;0;0], [0;0;1], [0;-1;0]];
R.R3 = [[-1;0;0], [0;0;1], [0;1;0]];
R.R4 = [[-1;0;0], [0;0;-1], [0;-1;0]];

de1 = 1.0e-3; de2 = 2.0e-3; de3 = 2.1e-3;
d_list = [[-de2/2; de3/2; -de1],...
    [-de2/2; de3/2; 0],...
    [-de2/2; de3/2; de1],...
    [-de2/2; -de3/2; -de1],...
    [-de2/2; -de3/2; 0],...
    [-de2/2; -de3/2; de1],...
    [de2/2; de3/2; -de1],...
    [de2/2; de3/2; 0],...
    [de2/2; de3/2; de1],...
    [de2/2; -de3/2; -de1],...
    [de2/2; -de3/2; 0],...
    [de2/2; -de3/2; de1]];

% figure
% plot3(0, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
% hold on;
% draw_axes([0;0;0], eye(3));
% draw_axes(d_list(:,1), R.R1);
% % draw_axes(d_list(:,2), R.R1);
% % draw_axes(d_list(:,3), R.R1);
% draw_axes(d_list(:,4), R.R2);
% % draw_axes(d_list(:,5), R.R2);
% % draw_axes(d_list(:,6), R.R2);
% draw_axes(d_list(:,7), R.R3);
% % draw_axes(d_list(:,8), R.R3);
% % draw_axes(d_list(:,9), R.R3);
% draw_axes(d_list(:,10), R.R4);
% % draw_axes(d_list(:,11), R.R4);
% % draw_axes(d_list(:,12), R.R4);
% scatter3(d_list(1,:), d_list(2,:), d_list(3,:), 'filled');
%%

[b_total, p_true, R_true, m_hat, m_pos] = process_matrix(data, R);

load('scan_records_bg.mat', 'scan_records');
data = scan_records;
b_total_bg = process_matrix(data, R);

b_msr = b_total - b_total_bg;
b_msr = b_msr * 1e-4; % 从Gs转换为特斯拉
%%
params = get_experiment_params();
params.sensor.d_list = d_list;
params.magnet.m_hat = m_hat;
params.magnet.m_pos = m_pos;

%% 建立一个基于err的下降搜索来估磁矩大小
m_norm_init = 0 * ones(1, 2);
lb = -1000;
ub = 1000;

objective = @(m) calculate_total_err(m, p_true, R_true, b_msr, params);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
% [m_norm_opt, min_err] = fmincon(objective, m_norm_init, [], [], [], [], lb, ub, [], options);

% params.magnet.m_norm = m_norm_opt;
% fprintf('Optimized m_norm: [%.2f, %.2f], Total Error: %.4f\n', m_norm_opt(1), m_norm_opt(2), min_err);

% Optimized m_norm: [75.49, 132.31], Total Error: 1.6521
% params.magnet.m_norm = [75.49, 132.31];
%%
objective(-89.29)

% figure;
% plot3(0, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
% hold on;
% % quiver3(pos_mat(1,1), pos_mat(2,1), pos_mat(3,1), orn_mat(1,1), orn_mat(2,1), orn_mat(3,1), 'k');
% scatter3(pos_mat(1,:), pos_mat(2,:), pos_mat(3,:));
% R = quat2rotm([orn_mat(4,1), orn_mat(1:3,1)']);

% quiver3(pos_mat(1,1), pos_mat(2,1), pos_mat(3,1), R(1,3), R(2,3), R(3,3), 'k');
% quiver3(m_pos(1,1), m_pos(2,1), m_pos(3,1), m_hat(1,1), m_hat(2,1), m_hat(3,1), 'r', 'LineWidth', 2);
% quiver3(m_pos(1,2), m_pos(2,2), m_pos(3,2), m_hat(1,2), m_hat(2,2), m_hat(3,2), 'b', 'LineWidth', 2);
% axis equal
%%
function [b_total, p_true, R_true, m_hat, m_pos] = process_matrix(data, R)

    % Rsa1 = eul2rotm([0.029868, 0.0392182, -2.4342], 'XYZ');
    % Tsa1 = RpToTrans(Rsa1,[1.42402, 0.400196, 0.0528185]');

    % Rsa2 = eul2rotm([-0.0121089, 0.00823125, 2.6721], 'XYZ');
    % Tsa2 = RpToTrans(Rsa2,[1.33131, -0.650143, 0.0313949]');  

    N = numel(data);
    p_true = nan(3, N);
    b_total = nan(3, 12, N);
    R_true = nan(3, 3, N);
    mag1_pos_mat = nan(3, N);
    mag2_pos_mat = nan(3, N);
    mag1_orn_mat = nan(4, N);
    mag2_orn_mat = nan(4, N);

    for i = 1:N
        p_true(:, i) = data(i).poses.diana7.position;
        R_true(:, :, i) = quat2rotm([data(i).poses.diana7.orientation(4), data(i).poses.diana7.orientation(1), data(i).poses.diana7.orientation(2), data(i).poses.diana7.orientation(3)]);
        mag1_pos_mat(:, i) = data(i).poses.arm1.position;
        mag2_pos_mat(:, i) = data(i).poses.arm2.position;
        mag1_orn_mat(:, i) = data(i).poses.arm1.orientation;
        mag2_orn_mat(:, i) = data(i).poses.arm2.orientation;
        b_total(:, :, i) = data(i).hall_matrix;
        b_total(:, 1:3, i) = R.R1 * b_total(:, 1:3, i);
        b_total(:, 4:6, i) = R.R2 * b_total(:, 4:6, i);
        b_total(:, 7:9, i) = R.R3 * b_total(:, 7:9, i);
        b_total(:, 10:12, i) = R.R4 * b_total(:, 10:12, i);
    end

    m_pos = [mean(mag1_pos_mat, 2), mean(mag2_pos_mat, 2)];
    mag1_orn = mean(mag1_orn_mat, 2);
    mag2_orn = mean(mag2_orn_mat, 2);

    Rmag1 = quat2rotm([mag1_orn(4), mag1_orn(1), mag1_orn(2), mag1_orn(3)]);
    Rmag2 = quat2rotm([mag2_orn(4), mag2_orn(1), mag2_orn(2), mag2_orn(3)]);

    m_hat = [Rmag1*[0;0;1], Rmag2*[0;0;1]];

end

function total_err = calculate_total_err(m_norm, p_true, R_true, b_msr, params)
    params.magnet.m_norm = m_norm * ones(1,2);
    N = size(p_true, 2);
    errs = zeros(1, N);
    for i = 1:N
        params.ground_truth.p_true = p_true(:, i);
        params.ground_truth.theta_true = so3ToVec(MatrixLog3((R_true(:, :, i))'));
        b_total = generate_magnetic_data(params);
        errs(i) = sum(vecnorm(b_msr(:,:,i) - b_total)) / size(b_msr, 2);
    end
    total_err = sum(errs);
end

function draw_axes(p, R)
    quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), 1e-3, 'r');
    quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), 1e-3, 'g');
    quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), 1e-3, 'b');
end