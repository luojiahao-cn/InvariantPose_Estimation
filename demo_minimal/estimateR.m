function [R_opt1, R_opt2] = estimateR(b_bar, B_bar, A_p, X, D_delta, B_delta)
% ESTIMATER Estimate rotation matrix via eigen decomposition
% Solves R^T A R = X for R using eigenvalue matching

    A_p = 0.5 * (A_p + A_p');
    X   = 0.5 * (X   + X');

    [PA, LA] = eig(A_p);
    [PX, LX] = eig(X);

    [~, IndA] = sort(diag(LA), 'descend');
    PA = PA(:, IndA);
    [~, IndX] = sort(diag(LX), 'descend');
    PX = PX(:, IndX);

    perms3 = perms(1:3);
    bestCost = inf;
    bestP = eye(3);

    for k = 1:size(perms3,1)
        p = perms3(k,:);
        P = eye(3);
        P = P(:, p);

        Y = P;
        R_cand = PA * Y * PX';
        cost = norm(R_cand' * A_p * R_cand - X, 'fro')^2;

        if cost < bestCost
            bestCost = cost;
            bestP = P;
        end
    end

    R_opt1 = PA * bestP * PX';
    R_opt1 = PA * diag(sign(diag(PX'*PA)))*PX';

    if det(R_opt1) < 0
        R_opt1(:,3) = -R_opt1(:,3);
    end

    % SVD-based estimate
    M = B_bar * b_bar';
    [U, ~, V] = svd(M);
    R_opt2 = V*diag([1,1,det(V*U')])*U';
end
