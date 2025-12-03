function [R_opt1, R_opt2] = estimateR(b_bar, B_bar, A_p, X)
% estimate R using R^T A R = X
[PA, LA] = eig(A_p);
[PX, LX] = eig(X);
[lam, IndA] = sort(diag(LA), 'ascend');
PA = PA(:,IndA);
[sig, IndX] = sort(diag(LX), 'ascend');
PX = PX(:,IndX);

% ---- 规范为 SO(3) 的特征向量基 ----
if det(PA) < 0, PA(:,1) = -PA(:,1); end
if det(PX) < 0, PX(:,1) = -PX(:,1); end

% ---- 枚举 6 个置换，做最佳匹配 ----
perms3 = perms(1:3);
bestCost = inf; bestP = eye(3);
for k = 1:size(perms3,1)
    p = perms3(k,:);
    P = eye(3); P = P(:, p);                 % 置换矩阵
    cost = sum((lam(p) - sig).^2);           % 等价于最大化 sum(sig .* lam(p))
    if cost < bestCost
        bestCost = cost; bestP = P;
    end

end

% ---- 保证 Y ∈ SO(3) ----
Y = bestP;
if det(Y) < 0
    Y = Y * diag([1 1 -1]);                  % 用符号阵把行列式调成 +1
end

% ---- 回代得到最优解 R ----
R_opt1 = PA * Y * PX';
% R_opt1 = PA * diag(sign(diag(PX'*PA)))*PX';

% estimate R using SVD
M = B_bar * b_bar';
[U, ~, V] = svd(M);
R_opt2 = V*diag([1,1,det(V*U')])*U';
end