function [r_hat, info] = estimate_principal_axis_SDP_reduced(prob, opts)
% estimate_principal_axis_SDP
% Principal-axis estimation via SDP relaxation + recovery (paper-matching form)
%
% Paper's reduced SDP (variables: R, {tau_k}, t):
%   minimize   sum_k alpha_k * tau_k + tr(M*R) - 2*beta*t
%   s.t.       tr(R) = 1,  R >= 0
%              [1, tr(A^k R); tr(A^k R), tau_k] >= 0,  forall k
%              [g' R g, t; t, 1] >= 0,  t >= 0
%
% where  M = beta*G - 2 * sum_k alpha_k * s_k * A^k.
%
% Recovery:
%   if rank(R*)==1: r from factorization (up to sign)
%   else if g'R*g > 0: r_ellip = R*g / sqrt(g'R*g), then normalize to sphere
%   else: r = TopEigenVector(R)
%   finally: choose sign that minimizes original quartic objective
%
% Inputs:
%   prob.Ak    : cell array {A^k} (each 3x3 symmetric)
%   prob.alpha : vector alpha_k
%   prob.s     : vector s_k
%   prob.beta  : scalar beta
%   prob.G     : 3x3 symmetric
%   prob.g     : 3x1 vector
%
% opts fields (optional):
%   opts.quiet     (default true)
%   opts.rank_tol  (default 1e-6)  % for numerical rank / eigen gap
%   opts.eps_gRg   (default 1e-12) % degeneracy threshold for g'R*g
%
% Output:
%   r_hat : 3x1 unit vector
%   info  : struct with solver status, R*, tau*, t*, costs, etc.

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'quiet'),    opts.quiet = true; end
if ~isfield(opts, 'rank_tol'), opts.rank_tol = 1e-6; end
if ~isfield(opts, 'eps_gRg'),  opts.eps_gRg = 1e-12; end

% ---------------- unpack problem data ----------------
Ak    = prob.Ak;
alpha = prob.alpha(:);
s     = prob.s(:);
beta  = prob.beta;
G     = prob.G;
g     = prob.g(:);

K = numel(Ak);
assert(K == numel(alpha) && K == numel(s), 'Size mismatch: Ak/alpha/s');

% ---------------- build M (paper definition) ----------------
M = beta * G;
for k = 1:K
    M = M - 2 * alpha(k) * s(k) * Ak{k};
end

% ---------------- CVX: solve reduced SDP ----------------
if ~opts.quiet
    fprintf('[SDP] Starting CVX solver (reduced form)...\n');
end

cvx_begin sdp
    if opts.quiet
        cvx_quiet(true);
    end

    variable R(3,3) symmetric
    variable tau(K,1)
    variable t(1,1)

    minimize( alpha' * tau + trace(M * R) - 2 * beta * t )

    subject to
        trace(R) == 1;
        R >= 0;
        t >= 0;

        % tau_k >= (tr(A^k R))^2  via Schur complement
        for k = 1:K
            yk = trace(Ak{k} * R);
            [1,  yk;
             yk, tau(k)] >= 0;
        end

        % t^2 <= g' R g  via Schur complement
        q = quad_form(g, R); % equals g' * R * g
        [q, t;
         t, 1] >= 0;
cvx_end

% ---------------- recovery ----------------
% Eigen-decomposition for rank / fallback
[V, D] = eig(full(R));
eigvals = real(diag(D));
[eigvals_sorted, idx_sort] = sort(eigvals, 'descend');
v1 = V(:, idx_sort(1));
lambda1 = eigvals_sorted(1);
lambda2 = eigvals_sorted(min(2, numel(eigvals_sorted)));

% Numerical rank-1 check: eigen gap ratio
% (You can tune this; for 3x3 it is usually stable.)
is_rank1 = (lambda2 <= opts.rank_tol * max(lambda1, 1));

gRg = full(g' * R * g);

if is_rank1
    % If rank-1, r is principal eigenvector (sign undecided)
    r_tilde = v1 / max(norm(v1), eps);
else
    if gRg > opts.eps_gRg
        % Lemma-driven ellipsoidal maximizer, then project to sphere
        r_ellip = (R * g) / sqrt(gRg);
        r_tilde = r_ellip / max(norm(r_ellip), eps);
    else
        % Degenerate case: no reliable linear bias -> use TopEigenVector(R)
        r_tilde = v1 / max(norm(v1), eps);
    end
end

% ---------------- sign selection by evaluating original quartic objective ----------------
J_pos = quartic_cost(r_tilde, Ak, alpha, M, beta, g);
J_neg = quartic_cost(-r_tilde, Ak, alpha, M, beta, g);

if J_pos <= J_neg
    r_hat = r_tilde;
else
    r_hat = -r_tilde;
end

% ---------------- pack info ----------------
info = struct();
info.solver_status = cvx_status;
info.optval_sdp    = cvx_optval;
info.R             = full(R);
info.tau           = full(tau);
info.t             = full(t);
info.M             = M;

info.eigvals_R      = eigvals_sorted;
info.is_rank1        = is_rank1;
info.gRg             = gRg;

info.r_tilde         = r_tilde;
info.cost_pos        = J_pos;
info.cost_neg        = J_neg;
info.cost_est        = min(J_pos, J_neg);

end

% ---------------- helper: original quartic objective ----------------
function J = quartic_cost(r, Ak, alpha, M, beta, g)
r = r(:);
r = r / max(norm(r), eps);

% quartic sum: sum alpha_k * (r' A^k r)^2
Jq = 0;
for k = 1:numel(Ak)
    qk = r' * Ak{k} * r;
    Jq = Jq + alpha(k) * (qk^2);
end

% quadratic term: r' M r
J2 = r' * M * r;

% linear term: -2 beta g' r
J1 = -2 * beta * (g' * r);

J = Jq + J2 + J1;
end
