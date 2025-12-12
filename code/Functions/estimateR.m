function [R_opt1, R_opt2] = estimateR(b_bar, B_bar, A_p, X, D_delta, B_delta)
    % 估计 R 使 R^T A R = X ，仅返回主解 R_opt1
    
    % --- 强制对称（数值上更稳） ---
    A_p = 0.5 * (A_p + A_p');
    X   = 0.5 * (X   + X');
    
    % --- 特征分解 ---
    [PA, LA] = eig(A_p);
    [PX, LX] = eig(X);
    
    % --- 排序特征值并重排特征向量 ---
    [~, IndA] = sort(diag(LA), 'descend');
    PA = PA(:, IndA);
    [~, IndX] = sort(diag(LX), 'descend');
    PX = PX(:, IndX);
    
    % --- 不对 PA / PX 强行调符号（避免时间跳变），符号只在最终 R 上调 ---
    % --- 枚举 6 个排列，选最优 ---
    perms3 = perms(1:3);
    bestCost = inf;
    bestP = eye(3);
    
    for k = 1:size(perms3,1)
        p = perms3(k,:);
        P = eye(3); 
        P = P(:, p);
    
        % 代价直接比较矩阵，不只比较特征值，数值更稳
        Y = P;
        R_cand = PA * Y * PX';
        cost = norm(R_cand' * A_p * R_cand - X, 'fro')^2;
        % cost = norm(R_cand' * A_p * R_cand * D_delta - B_delta, 'fro')^2;
    
        if cost < bestCost
            bestCost = cost;
            bestP = P;
        end
    end
    
    % --- 得到初步 R ---
    R_opt1 = PA * bestP * PX';

    R_opt1 = PA * diag(sign(diag(PX'*PA)))*PX';
    
    % --- 统一调整到 SO(3)，避免独立翻 PA/PX 的符号 ---
    if det(R_opt1) < 0
        % 翻转 R 的第一列，同时不破坏 R^T A R = X
        R_opt1(:,3) = -R_opt1(:,3);
    end
    
    % estimate R using SVD
    M = B_bar * b_bar';
    [U, ~, V] = svd(M);
    R_opt2 = V*diag([1,1,det(V*U')])*U';
    
end