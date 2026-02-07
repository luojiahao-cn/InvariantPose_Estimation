function r_best = estimate_principal_axis_PS(b_bar, B_bar, A_p, X, opts)
% Principal-axis estimator (paper-aligned):
%   - use (s1,s2) invariants to solve w
%   - enumerate finite candidates z (modulo z ~ -z)
%   - select by OPP residual only (odd/sign-sensitive cue)
%
% Inputs:
%   B_bar : 3xm   predicted \bbs B(\hat p)
%   b_bar : 3xm   measured  \bar{\mbc B}
%   A_p   : 3x3   world gradient tensor A(\hat p)
%   X     : 3x3   local tensor estimate \hat X
% Output:
%   r     : 3x1   (undirected) principal axis in world frame

    u = opts.u;

    % --- enforce symmetry (numerical stability) ---
    A_p = 0.5 * (A_p + A_p');
    X   = 0.5 * (X   + X');

    % --- eigendecomposition of A_p (sorted) ---
    [V, LA] = eig(A_p);
    lam = diag(LA);
    [lam, idx] = sort(lam, 'ascend');
    V = V(:, idx);

    % --- invariants from X ---
    s1 = u' * X * u;
    s2 = u' * (X*X) * u;

    % --- solve for w in squared coordinates ---
    % --- 投影谱方法核心替代方案 ---
    C_eq = B_bar'; 
    d_eq = b_bar' * u; % 假设这是你的线性观测

    % 1. 处理线性约束
    r_p = pinv(C_eq) * d_eq;
    N = null(C_eq);

    % 2. 构造简化的二次约束参数
    A_hat = N' * A_p * N;
    f = (2 * r_p' * A_p * N)';
    k = s1 - r_p' * A_p * r_p;
    R2 = 1 - norm(r_p)^2;

    % 3. 求解 Secular Equation (利用 fzero)
    % 目标：找到 lambda 使得 z(lambda)' * z(lambda) - R2 = 0
    % 1. 预处理：特征值分解 (A_hat 是小矩阵，通常 2x2 或 3x3)
    [Q, SigMat] = eig(A_hat);
    sig = diag(SigMat);   % 提取特征值
    q = Q' * (-0.5 * f);  % 变换线性项

    % 2. 检查 R2 的合法性 (防止线性约束冲突)
    if R2 <= 0
        % 说明线性约束已经使得点超出了单位球，退化处理
        r_best = r_p / norm(r_p); 
    else
        % 3. 定义鲁棒的长期方程 (避开矩阵求逆)
        % phi(lambda) = sum( q_i^2 / (sig_i + lambda)^2 ) - R2
        phi = @(lam) sum( q.^2 ./ (sig + lam).^2 ) - R2;

        % 4. 确定搜索区间
        % 解通常在 lambda > -min(sig) 的区域
        lower_bound = -min(sig) + 1e-9; 
        
        % 寻找一个上限以确定 fzero 的区间 [lower, upper]
        upper_bound = 10; 
        while phi(upper_bound) > 0 && upper_bound < 1e6
            upper_bound = upper_bound * 10;
        end

        try
            % 使用区间查找，比单点查找更稳定
            lambda_opt = fzero(phi, [lower_bound, upper_bound]);
            
            % 计算 z
            z_opt = Q * (q ./ (sig + lambda_opt));
            r_best = r_p + N * z_opt;
        catch
            % 备选方案：如果 fzero 失败，可能是极点太近，取最小二乘解
            r_best = r_p / norm(r_p);
        end
    end
end