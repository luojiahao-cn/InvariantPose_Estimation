function [r_hat, info] = estimate_principal_axis_MM(prob, opts)
% estimate_principal_axis_MM
% principal-axis estimation on S^2 via structured MM + exact TRS.
%
% Solve (||r||=1):
%   min  sum_k alpha_k ( r' A^k r - s_k )^2  +  beta || B_bar' r - b_bar' u ||^2
%
% Inputs:
%   prob : Problem data from compute_principal_axis_prob.m
%   opts : Struct for solver parameters

if nargin < 2, opts = struct(); end

% ---------------- solver options ----------------
if ~isfield(opts,'maxIter') || isempty(opts.maxIter), opts.maxIter = 100; end
maxIter = opts.maxIter;

if ~isfield(opts,'tol') || isempty(opts.tol), opts.tol = 1e-8; end
tol = opts.tol;

if ~isfield(opts,'lambdaTol') || isempty(opts.lambdaTol), opts.lambdaTol = 1e-12; end
lambdaTol = opts.lambdaTol;

if ~isfield(opts,'lambdaMax') || isempty(opts.lambdaMax), opts.lambdaMax = 1e12; end
lambdaMax = opts.lambdaMax;

if ~isfield(opts,'proxEta') || isempty(opts.proxEta), opts.proxEta = 0; end
proxEta = opts.proxEta;

if ~isfield(opts,'checkSign') || isempty(opts.checkSign), opts.checkSign = false; end
checkSign = opts.checkSign;

if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = 0; end
verbose = opts.verbose;

if ~isfield(opts,'trsHardTol') || isempty(opts.trsHardTol), opts.trsHardTol = 1e-10; end
trsHardTol = opts.trsHardTol;

% ---------------- unpack prob ----------------
Ak    = prob.Ak;
s     = prob.s;
alpha = prob.alpha;
beta  = prob.beta;
G     = prob.G;
g     = prob.g;
r     = prob.r0;
K     = prob.K;
u     = prob.u;

% ---------------- MM iterations ----------------

cost_hist   = zeros(maxIter+1,1);
step_hist   = zeros(maxIter,1);
lambda_hist = nan(maxIter,1);

if ~isempty(prob.r_true)
    cost_true = full_cost(prob.r_true, Ak, alpha, s, beta, G, g);
else
    cost_true = NaN;
end

cost_hist(1) = full_cost(r, Ak, alpha, s, beta, G, g);
converged = false;
t_end = 0;

for t = 1:maxIter
    r_t = r;

    % ---- build structured MM surrogate: H^(t), f^(t) ----
    H = beta * G;
    f = beta * g;

    for j = 1:numel(K)
        a = Ak{j} * r_t;              % a_k^(t) = A^k r^(t)
        c = r_t.' * a;               % c_k^(t) = r^(t)' A^k r^(t)

        H = H + 4 * alpha(j) * (a * a.');              % rank-1 PSD
        f = f + alpha(j) * (2*c + s(j)) * a;           % linear term
    end


    % Proximal: (eta/2)||r - r_t||^2
    % Under ||r||=1, variable part is -eta r_t' r => f += (eta/2) r_t.
    if proxEta > 0
        f = f + (proxEta/2) * r_t;
    end

    H = 0.5*(H + H.');  % enforce symmetry

    % ---- Solve TRS subproblem exactly ----
    [r_new, lambda] = solve_trs_sphere(H, f, lambdaTol, lambdaMax, trsHardTol);

    % Optional sign check w.r.t. ORIGINAL cost
    if checkSign
        if full_cost(-r_new, Ak, alpha, s, beta, G, g) < ...
           full_cost( r_new, Ak, alpha, s, beta, G, g)
            r_new = -r_new;
        end
    end

    step = norm(r_new - r_t);
    step_hist(t) = step;
    lambda_hist(t) = lambda;

    r = r_new;
    cost_hist(t+1) = full_cost(r, Ak, alpha, s, beta, G, g);


    if verbose
        fprintf('[MM+TRS] iter %d: cost=%.6e, step=%.3e, lambda=%.6e\n', ...
            t, cost_hist(t+1), step, lambda);
    end

    t_end = t;

    if step < tol
        converged = true;
        break;
    end
end

r_hat = r / max(norm(r), eps);

% pack info
info = struct();
info.converged   = converged;
info.iters       = t_end;
info.cost_hist   = cost_hist(1:t_end+1);
info.step_hist   = step_hist(1:t_end);
info.lambda_hist = lambda_hist(1:t_end);
info.cost_final  = info.cost_hist(end);
info.K           = K;
info.alpha       = alpha;
info.beta        = beta;
info.u           = u;
info.proxEta     = proxEta;
info.cost_true   = cost_true;

end

% ---------------- helper: TRS on sphere (exact) ----------------
function [r, lambda] = solve_trs_sphere(H, f, lambdaTol, lambdaMax, hardTol)
% Solve: min_{||r||=1} r'Hr - 2 f'r
% KKT: (H + lambda I) r = f, ||r||=1, lambda > -lambda_min(H)
%
% Robust handling of the "hard case".

f = f(:);
fn = norm(f);

% If f is tiny, pick smallest eigenvector
if fn < 1e-14
    [U,D] = eig(H);
    [~,iMin] = min(diag(D));
    r = U(:,iMin);
    r = r / max(norm(r), eps);
    lambda = 0;
    return;
end

[U,D] = eig(H);
d = diag(D);
c = U.' * f;

% Identify minimal eigenspace
[dmin, idxMin] = min(d);
idxNull = find(abs(d - dmin) < 1e-12);

% Lower bound (strict) to keep invertible
lam_lo = -dmin + 1e-12;

phi = @(lam) sum((c.^2) ./ (d + lam).^2) - 1;

% Check hard-case condition: f orthogonal to min-eigenspace
if ~isempty(idxNull) && norm(c(idxNull)) < hardTol
    % Evaluate at lam_lo ~ -dmin: finite norm because null components are zero
    val = phi(lam_lo);
    if val <= 0
        % Hard case: lambda = -dmin, r = r0 + tau v_min
        lambda = -dmin;

        y = zeros(3,1);
        for i = 1:3
            if all(i ~= idxNull)
                y(i) = c(i) / (d(i) - dmin);
            end
        end

        r0 = U * y;
        nr0 = norm(r0);

        % choose tau so that ||r||=1
        tau = sqrt(max(0, 1 - nr0^2));
        v = U(:, idxNull(1)); % pick one basis vector in minimal eigenspace

        r_plus  = r0 + tau * v;
        r_minus = r0 - tau * v;

        % pick the better one for objective
        Jp = r_plus.'*H*r_plus - 2*f.'*r_plus;
        Jm = r_minus.'*H*r_minus - 2*f.'*r_minus;

        if Jm < Jp
            r = r_minus;
        else
            r = r_plus;
        end

        r = r / max(norm(r), eps);
        return;
    end
end

% Regular case: unique root in (lam_lo, +inf)
% Bracket an upper bound
lam_hi = max(lam_lo + 1, 1);
while lam_hi < lambdaMax && phi(lam_hi) > 0
    lam_hi = 2 * lam_hi;
end
if lam_hi >= lambdaMax && phi(lam_hi) > 0
    % Fallback (should be rare)
    lambda = lambdaMax;
else
    % Bisection
    for it = 1:200
        lambda = 0.5*(lam_lo + lam_hi);
        if phi(lambda) > 0
            lam_lo = lambda;
        else
            lam_hi = lambda;
        end
        if (lam_hi - lam_lo) < lambdaTol * max(1, abs(lam_hi))
            break;
        end
    end
end

% Reconstruct r
r = U * (c ./ (d + lambda));
% With correct lambda, ||r|| should be 1 (up to numerical error)
r = r / max(norm(r), eps);

end

% ---------------- helper: full cost ----------------
function J = full_cost(rr, Ak, alpha, s, beta, G, g)
rr = rr(:); rr = rr / max(norm(rr), eps);
J1 = 0;
for j = 1:numel(Ak)
    q = rr.' * Ak{j} * rr;
    J1 = J1 + alpha(j) * (q - s(j))^2;
end
J2 = beta * (rr.' * G * rr - 2 * g.' * rr);
J = J1 + J2;
end

