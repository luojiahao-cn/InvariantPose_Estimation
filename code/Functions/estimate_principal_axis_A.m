function [r_hat, info] = estimate_principal_axis_A(b_bar, B_bar, A_p, X, opts)
% estimate_principal_axis_A
% A-scheme: joint estimation on S^2 via MM (majorization-minimization).
%
% This implements the "MM joint descent" route:
%   minimize over ||r||=1:
%       sum_k alpha_k ( r' A^k r - s_k )^2  +  beta || B_bar' r - b_bar' u ||^2
% where s_k = u' X^k u.
%
% Inputs
%   B_bar : 3xm matrix, predicted projected field matrix \bar{B}(\hat p)
%   b_bar : 3xm matrix, measured projected/processed \bar{B}
%   A_p   : 3x3 symmetric world gradient tensor A(\hat p)
%   X     : 3x3 symmetric local estimated tensor \hat X
%   opts  : struct with optional fields:
%       .u            : 3x1 local principal axis (default [1;0;0])
%       .K            : vector of positive integers, powers to use (default [1 2 3])
%       .alpha        : weights for each k in K (default ones)
%       .beta         : weight for field-consistency term (default 1)
%       .maxIter      : max iterations (default 50)
%       .tol          : stopping tol on ||r_{t+1}-r_t|| (default 1e-8)
%       .lambdaTol    : bisection tol for lambda (default 1e-12)
%       .lambdaMax    : max lambda for bisection (default 1e12)
%       .proxEta      : proximal gain (default 0)  adds (eta/2)||r-r_t||^2
%       .checkSign    : if true, choose between r and -r by full cost (default true)
%       .r0           : initial r in world frame (3x1). If empty, uses A_p eig init.
%       .verbose      : 0/1 (default 0)
%
% Outputs
%   r_hat : 3x1 unit vector, estimated principal axis in world frame
%   info  : struct with fields:
%       .converged, .iters, .cost_hist, .step_hist, .lambda_hist
%
% Notes
%   - A_p and X are symmetrized for numerical stability.
%   - If the linear term is near zero, the subproblem reduces to an eigenvector pick.
%
% (c) MM refinement for principal-axis estimation.

% ---------------- options ----------------
if nargin < 5, opts = struct(); end
if ~isfield(opts,'u') || isempty(opts.u), opts.u = [1;0;0]; end
u = opts.u(:); u = u / max(norm(u), eps);

if ~isfield(opts,'K') || isempty(opts.K), opts.K = [1, 2]; end
K = opts.K(:)';

if ~isfield(opts,'alpha') || isempty(opts.alpha), opts.alpha = ones(size(K)); end
alpha = opts.alpha(:)'; 
if numel(alpha) ~= numel(K)
    error('opts.alpha must have same length as opts.K.');
end

if ~isfield(opts,'beta') || isempty(opts.beta), opts.beta = 1; end
beta = opts.beta;

if ~isfield(opts,'maxIter') || isempty(opts.maxIter), opts.maxIter = 50; end
maxIter = opts.maxIter;

if ~isfield(opts,'tol') || isempty(opts.tol), opts.tol = 1e-8; end
tol = opts.tol;

if ~isfield(opts,'lambdaTol') || isempty(opts.lambdaTol), opts.lambdaTol = 1e-12; end
lambdaTol = opts.lambdaTol;

if ~isfield(opts,'lambdaMax') || isempty(opts.lambdaMax), opts.lambdaMax = 1e12; end
lambdaMax = opts.lambdaMax;

if ~isfield(opts,'proxEta') || isempty(opts.proxEta), opts.proxEta = 0; end
proxEta = opts.proxEta;

if ~isfield(opts,'checkSign') || isempty(opts.checkSign), opts.checkSign = true; end
checkSign = opts.checkSign;

if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = 0; end
verbose = opts.verbose;

% ---------------- sanitize tensors ----------------
A_p = 0.5*(A_p + A_p');
X   = 0.5*(X   + X');

% ---------------- build field term matrices ----------------
% Field-consistency residual: || B_bar' r - b_bar' u ||^2
% => quadratic: r' G r - 2 g' r + const
G = B_bar * B_bar.';                 % 3x3
b_u = b_bar.' * u;                   % mx1
g = B_bar * b_u;                     % 3x1

% ---------------- precompute M_k = A^k and s_k = u' X^k u ----------------
M = cell(1, numel(K));
s = zeros(1, numel(K));

% Precompute A^k efficiently
A_pow = eye(3);
X_pow = eye(3);
kmax = max(K);
A_pows = cell(1, kmax);
X_pows = cell(1, kmax);

for kk = 1:kmax
    A_pow = A_pow * A_p;
    X_pow = X_pow * X;
    A_pows{kk} = 0.5*(A_pow + A_pow');  % keep symmetric
    X_pows{kk} = 0.5*(X_pow + X_pow');
end

for i = 1:numel(K)
    ki = K(i);
    M{i} = A_pows{ki};
    s(i) = u.' * X_pows{ki} * u;
end

% ---------------- initialization r0 ----------------
if isfield(opts,'r0') && ~isempty(opts.r0)
    r = opts.r0(:);
    r = r / max(norm(r), eps);
else
    % A simple, stable init: use eigenvector of A_p with largest |eigenvalue|
    [V, D] = eig(A_p);
    [~, idx] = max(abs(diag(D)));
    r = V(:, idx);
    r = r / max(norm(r), eps);
end

% ---------------- MM iterations ----------------
cost_hist   = zeros(maxIter+1,1);
step_hist   = zeros(maxIter,1);
lambda_hist = nan(maxIter,1);

cost_hist(1) = full_cost(r, M, K, alpha, s, beta, B_bar, b_u);

converged = false;

for t = 1:maxIter
    % build H_t = beta*G + 2 sum alpha_k (q_k(r_t)-s_k) M_k
    H = beta * G;
    for j = 1:numel(K)
        qjt = r.' * M{j} * r;
        H = H + 2 * alpha(j) * (qjt - s(j)) * M{j};
    end

    f = beta * g;

    % optional proximal term: (eta/2)||r - r_t||^2
    if proxEta > 0
        H = H + (proxEta/2) * eye(3);
        f = f + (proxEta/2) * r;
    end

    H = 0.5*(H + H');  % enforce symmetry

    % Solve subproblem: min_{||r||=1} r' H r - 2 f' r
    fnorm = norm(f);

    if fnorm < 1e-14
        % No linear term: pick eigenvector of smallest eigenvalue
        [U, D] = eig(H);
        [~, idxMin] = min(diag(D));
        r_new = U(:, idxMin);
        lambda = 0;
    else
        % Use eigen-basis + bisection on lambda:
        % (H + lambda I) r = f, ||r||=1
        [U, D] = eig(H);
        d = diag(D);
        c = U.' * f;

        % lower bound must keep denominators nonzero: lambda > -min(d)
        lam_lo = max(0, -min(d) + 1e-14);

        % function phi(lam) = ||r(lam)||^2 - 1
        phi = @(lam) sum((c.^2) ./ (d + lam).^2) - 1;

        % Find an upper bound where phi <= 0
        lam_hi = lam_lo;
        val_lo = phi(lam_lo);

        if val_lo <= 0
            % already satisfies ||r||<=1; use lam_lo
            lambda = lam_lo;
        else
            lam_hi = lam_lo + 1;
            while lam_hi < lambdaMax && phi(lam_hi) > 0
                lam_hi = 2 * lam_hi;
                if lam_hi == 0, lam_hi = 1; end
            end
            if lam_hi >= lambdaMax && phi(lam_hi) > 0
                % fallback: large lambda
                lambda = lambdaMax;
            else
                % bisection
                for it = 1:200
                    lambda = 0.5*(lam_lo + lam_hi);
                    if phi(lambda) > 0
                        lam_lo = lambda;
                    else
                        lam_hi = lambda;
                    end
                    if (lam_hi - lam_lo) < lambdaTol * max(1, lam_hi)
                        break;
                    end
                end
            end
        end

        % reconstruct r
        r_new = U * (c ./ (d + lambda));
        r_new = r_new / max(norm(r_new), eps);
    end

    % Optional sign check: choose between r and -r
    if checkSign
        if full_cost(-r_new, M, K, alpha, s, beta, B_bar, b_u) < full_cost(r_new, M, K, alpha, s, beta, B_bar, b_u)
            r_new = -r_new;
        end
    end

    step = norm(r_new - r);
    step_hist(t) = step;
    lambda_hist(t) = lambda;

    r = r_new;
    cost_hist(t+1) = full_cost(r, M, K, alpha, s, beta, B_bar, b_u);

    if verbose
        fprintf('[MM] iter %d: cost=%.6e, step=%.3e, lambda=%.3e\n', ...
            t, cost_hist(t+1), step, lambda);
    end

    if step < tol
        converged = true;
        break;
    end
end

r_hat = r / max(norm(r), eps);

% pack info
info = struct();
info.converged  = converged;
info.iters      = find(cost_hist~=0,1,'last') - 1;
if isempty(info.iters), info.iters = 0; end
info.cost_hist  = cost_hist(1:info.iters+1);
info.step_hist  = step_hist(1:max(0,info.iters));
info.lambda_hist= lambda_hist(1:max(0,info.iters));
info.K          = K;
info.alpha      = alpha;
info.beta       = beta;
info.u          = u;
info.proxEta    = proxEta;
end

% ---------------- helpers: full cost ----------------
function J = full_cost(rr, M, K, alpha, s, beta, B_bar, b_u)
    rr = rr(:); rr = rr / max(norm(rr), eps);
    J1 = 0;
    for j = 1:numel(K)
        q = rr.' * M{j} * rr;
        J1 = J1 + alpha(j) * (q - s(j))^2;
    end
    e = B_bar.' * rr - b_u;
    J2 = beta * (e.'*e);
    J = J1 + J2;
end