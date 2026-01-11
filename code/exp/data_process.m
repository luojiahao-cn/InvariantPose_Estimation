clc,clear,close all
%% Environment setup
this_file = mfilename('fullpath');
exp_dir = fileparts(this_file);
repo_root = fileparts(exp_dir);
addpath(fullfile(repo_root, 'utils'));
addpath(fullfile(repo_root, 'Functions'));
addpath(fullfile(repo_root, 'tools'));
%%
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
d_list = d_list - mean(d_list, 2);
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
load('scan_records.mat', 'scan_records');
[b_total, p_true, R_true, m_hat, m_pos, test_points, radius] = process_matrix(scan_records, R);

load('scan_records_bg.mat', 'scan_records');
b_total_bg = process_matrix(scan_records, R);

b_total = b_total - b_total_bg;
b_total = b_total * 1e-4; % 从Gs转换为特斯拉

[Q, ~, ~] = qr(d_list'); % Note: must return 3 outputs even we only need Q
r = rank(d_list); % 构型判据
Q_bar = Q(:, r+1:end);
g = Q_bar' * ones(size(d_list, 2), 1);

for i = 1:size(b_total, 3)
    b_mag_local(:, i) = b_total(:, :, i) * Q_bar * g / (norm(g)^2);
    b_mag(:, i) = R_true(:, :, i) * b_mag_local(:, i);
end
%%
params = get_experiment_params();
params.sensor.d_list = d_list;
params.magnet.m_hat = m_hat;
params.magnet.m_pos = m_pos;
params.workspace.radius = radius;

%% 建立一个基于err的下降搜索来估磁矩大小、位置和方向
% 优化策略：仅优化 m_norm (大小) 和 m_pos (位置)，方向固定
% 变量映射：
% x(1:2)   : m_norm (2)
% x(3:8)   : m_pos (6)
% 总变量数: 8

% 1. 初始化参数
% m_norm: 初始猜测 60 A·m²
params.magnet.m_norm = [60, 60]; 
m_norm_init = params.magnet.m_norm;

% m_pos
m_pos_init = params.magnet.m_pos;

% 构造优化向量 x [8维]
x0 = [m_norm_init(:); ...
      m_pos_init(:)];

% 2. 设置约束
% m_norm: [40, 80]
% m_pos:  初始值 +/- 5cm
lb = [40; 40; ...
      m_pos_init(:) - 0.08];

ub = [400; 400; ...
      m_pos_init(:) + 0.08];

% 目标函数
objective = @(x) calculate_total_err(x, p_true, b_mag, params);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 3e4, 'StepTolerance', 1e-8);
[x_opt, min_err] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);

% 3. 解析结果
params.magnet.m_norm = x_opt(1:2)';
params.magnet.m_pos  = reshape(x_opt(3:8), 3, 2);

% m_hat 基本保持不变

fprintf('Optimized m_norm: %.2f, %.2f\n', params.magnet.m_norm(1), params.magnet.m_norm(2));
fprintf('Optimized m_pos 1: [%.4f, %.4f, %.4f]\n', params.magnet.m_pos(:,1));
fprintf('Optimized m_pos 2: [%.4f, %.4f, %.4f]\n', params.magnet.m_pos(:,2));
% fprintf('Optimized m_hat 1 (deg): az=%.2f, el=%.2f\n', rad2deg(az1_opt), rad2deg(el1_opt));
% fprintf('Optimized m_hat 2 (deg): az=%.2f, el=%.2f\n', rad2deg(az2_opt), rad2deg(el2_opt));

[rmse, b_p_sim] = calculate_total_err(x_opt, p_true, b_mag, params);
fprintf('Final RMSE: %.4f T\n', rmse);
%%
S = euclidean_similarity(b_p_sim);
figure;
imagesc(S);
axis equal tight;
colorbar;
colormap(jet);
title('Similarity Matrix');
xlabel('Column Index');
ylabel('Column Index');

%%
% figure
% for i = 1:4:size(b_mag, 2)
%     scatter3(p_true(1, i), p_true(2, i), p_true(3, i), 'k');
%     hold on
%     quiver3(p_true(1, i), p_true(2, i), p_true(3, i), b_mag(1, i), b_mag(2, i), b_mag(3, i), 1e1, 'r');
%     quiver3(p_true(1, i), p_true(2, i), p_true(3, i), b_p_sim(1, i), b_p_sim(2, i), b_p_sim(3, i), 1e1, 'b');
% end

% save('optimized_params.mat', 'params', 'test_points', 'b_total', 'b_total_bg');

for idx = 1:size(b_total, 3)

    params_current = params;
    params_current.ground_truth.p_true = test_points(idx).p_true;
    params_current.ground_truth.theta_true = test_points(idx).theta_true;
    [b_total_sim, ~, ~, ~, ~] = generate_magnetic_data(params_current);

    err(idx) = sqrt(mean(vecnorm(b_total(:,:,idx) - b_total_sim).^2));
end

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
function [b_total, p_true, R_true, m_hat, m_pos, test_points, radius] = process_matrix(data, R)

    % Rsa1 = eul2rotm([0.029868, 0.0392182, -2.4342], 'XYZ');
    % Tsa1 = RpToTrans(Rsa1,[1.42402, 0.400196, 0.0528185]');

    % Rsa2 = eul2rotm([-0.0121089, 0.00823125, 2.6721], 'XYZ');
    % Tsa2 = RpToTrans(Rsa2,[1.33131, -0.650143, 0.0313949]');  

    N = numel(data);
    p_true = nan(3, N);
    b_total = nan(3, 12, N);
    R_true = nan(3, 3, N);
    theta_true = nan(3, N);
    mag1_pos_mat = nan(3, N);
    mag2_pos_mat = nan(3, N);
    mag1_orn_mat = nan(4, N);
    mag2_orn_mat = nan(4, N);
    test_points = struct('p_true', {}, 'theta_true', {});

    for i = 1:N
        p_true(:, i) = data(i).poses.diana7.position;
        R_true(:, :, i) = quat2rotm([data(i).poses.diana7.orientation(4), data(i).poses.diana7.orientation(1), data(i).poses.diana7.orientation(2), data(i).poses.diana7.orientation(3)]);
        theta_true(:, i) = so3ToVec(MatrixLog3(R_true(:, :, i)));
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
    workspace_center = mean(p_true, 2);
    p_true = p_true - workspace_center;
    m_pos = [mean(mag1_pos_mat, 2), mean(mag2_pos_mat, 2)];
    m_pos = m_pos - workspace_center;
    
    radius = max(vecnorm(p_true)); % extra 0.02m margin
    mag1_orn = mean(mag1_orn_mat, 2);
    mag2_orn = mean(mag2_orn_mat, 2);

    Rmag1 = quat2rotm([mag1_orn(4), mag1_orn(1), mag1_orn(2), mag1_orn(3)]);
    Rmag2 = quat2rotm([mag2_orn(4), mag2_orn(1), mag2_orn(2), mag2_orn(3)]);

    % 磁矩由S极指向N极。测得输出侧为S极，则N极在反向。
    m_hat = [-Rmag1*[0;0;1], -Rmag2*[0;0;1]];
    for i = 1:N
        test_points(i).p_true = p_true(:, i);
        test_points(i).theta_true = theta_true(:, i);
    end
end

function [rmse, b_p_sim] = calculate_total_err(x, p_true, b_mag, params)

    % 解析优化变量
    m_norm = x(1:2)';
    m_pos  = reshape(x(3:8), 3, 2);
    
    % m_hat 固定 (方向已锁死)
    m_hat = params.magnet.m_hat;

    N = size(p_true, 2);
    errs = zeros(1, N);
    b_p_sim = zeros(size(b_mag));

    for i = 1:N
        [b_p, ~] = calcFieldAndGradient( ...
            p_true(:, i), ...
            m_pos, ...
            m_hat, ...
            m_norm);
            
        diff = b_p - b_mag(:, i);
        errs(i) = norm(diff)^2;
        b_p_sim(:, i) = b_p;
    end

    rmse = sqrt(mean(errs));
end


function draw_axes(p, R)
    quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), 1e-3, 'r');
    quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), 1e-3, 'g');
    quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), 1e-3, 'b');
end

function S = euclidean_similarity(X)
% EUCLIDEAN_SIMILARITY  Compute similarity based on Euclidean distance
% Input:
%   X : 3×N matrix
% Output:
%   S : N×N similarity matrix in (0,1]

    N = size(X, 2);
    S = zeros(N, N);

    for i = 1:N
        for j = 1:N
            d = norm(X(:,i) - X(:,j), 2);
            S(i,j) = 1 / (1 + d);
        end
    end
end
