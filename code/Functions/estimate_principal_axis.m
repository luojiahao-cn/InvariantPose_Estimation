function r = estimate_principal_axis(b_bar, B_bar, A_p, X)
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

    u = [1;0;0];

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
    C = [ 1,        1,        1;
          lam(1),   lam(2),   lam(3);
          lam(1)^2, lam(2)^2, lam(3)^2 ];
    b = [1; s1; s2];

    % Unconstrained solve; if ill-conditioned, fall back to pinv
    if rcond(C) < 1e-12
        w = pinv(C) * b;
    else
        w = C \ b;
    end

    % Project onto simplex: w>=0, sum(w)=1
    w = project_to_simplex_3d(w);

    % --- enumerate candidates z modulo z ~ -z (<=4 candidates) ---
    sqrtw = sqrt(max(w, 0));

    S = dec2bin(0:7) - '0';    % 8x3
    S(S==0) = -1;              % -> ±1

    bestJ = inf;
    r_best = zeros(3,1);

    for ii = 1:size(S,1)
        sigma = S(ii,:).';
        z = sigma .* sqrtw;

        % eliminate z and -z: keep only one representative
        % rule: first nonzero component must be positive
        j = find(abs(z) > 1e-12, 1, 'first');
        if ~isempty(j) && z(j) < 0
            continue;
        end

        r_cand = V * z;
        nr = norm(r_cand);
        if nr < 1e-12
            continue;
        end
        r_cand = r_cand / nr;

        % --- OPP residual ONLY (sign-sensitive cue) ---
        J = norm(B_bar.' * r_cand - b_bar.' * u, 2)^2;

        if J < bestJ
            bestJ = J;
            r_best = r_cand;
        end
    end

    % output
    if norm(r_best) < 1e-12
        r = [1;0;0];  % fallback, should rarely happen
    else
        r = r_best / norm(r_best);
    end
end


% ====== helper: projection onto 3D simplex {w>=0, sum(w)=1} ======
function wproj = project_to_simplex_3d(w)
    w = w(:);

    v = sort(w, 'descend');
    cssv = cumsum(v) - 1;
    rho = find(v - cssv ./ (1:3).' > 0, 1, 'last');

    if isempty(rho)
        wproj = ones(3,1)/3;
        return;
    end

    theta = cssv(rho) / rho;
    wproj = max(w - theta, 0);

    s = sum(wproj);
    if s <= 0
        wproj = ones(3,1)/3;
    else
        wproj = wproj / s;
    end
end
