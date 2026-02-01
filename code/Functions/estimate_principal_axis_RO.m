function [r_hat, info] = estimate_principal_axis_RO(prob, opts)
% estimate_principal_axis_RO
% Riemannian Optimization (RO) on S^2 for principal-axis estimation.
%
% Solve (||r||=1):
%   min  sum_k alpha_k ( r' A^k r - s_k )^2  +  beta || B_bar' r - b_bar' u ||^2
%
% Inputs:
%   prob : Problem data from compute_principal_axis_prob.m
%   opts : Struct for solver parameters

if nargin < 2, opts = struct(); end

% ---------------- solver options ----------------
if ~isfield(opts,'maxIter') || isempty(opts,'maxIter'), opts.maxIter = 1000; end
maxIter = opts.maxIter;

if ~isfield(opts,'tol') || isempty(opts,'tol'), opts.tol = 1e-8; end
tol = opts.tol;

if ~isfield(opts,'checkSign') || isempty(opts,'checkSign'), opts.checkSign = false; end
checkSign = opts.checkSign;

if ~isfield(opts,'verbose') || isempty(opts,'verbose'), opts.verbose = 0; end
verbose = opts.verbose;

if ~isfield(opts,'step0') || isempty(opts,'step0'), opts.step0 = 1.0; end
step0 = opts.step0;

if ~isfield(opts,'ls_c1') || isempty(opts,'ls_c1'), opts.ls_c1 = 1e-4; end
ls_c1 = opts.ls_c1;

if ~isfield(opts,'ls_tau') || isempty(opts,'ls_tau'), opts.ls_tau = 0.5; end
ls_tau = opts.ls_tau;

if ~isfield(opts,'ls_maxSteps') || isempty(opts,'ls_maxSteps'), opts.ls_maxSteps = 25; end
ls_maxSteps = opts.ls_maxSteps;

if ~isfield(opts,'useBB') || isempty(opts,'useBB'), opts.useBB = true; end
useBB = opts.useBB;

if ~isfield(opts,'bbClamp') || isempty(opts,'bbClamp'), opts.bbClamp = [1e-6, 1e3]; end
bbClamp = opts.bbClamp(:)';

% ---------------- unpack prob ----------------
Ak    = prob.Ak;
s     = prob.s;
alpha = prob.alpha;
beta  = prob.beta;
G     = prob.G;
g     = prob.g;
r     = prob.r0;
u     = prob.u;
K     = prob.K;

% ---------------- RO iterations ----------------

cost_hist = zeros(maxIter+1, 1);
grad_hist = zeros(maxIter+1, 1);
step_hist = zeros(maxIter, 1);
ls_hist   = zeros(maxIter, 1);

cost_true = NaN;
if ~isempty(prob.r_true)
    cost_true = full_cost(prob.r_true, Ak, alpha, s, beta, G, g);
end

cost_hist(1) = full_cost(r, Ak, alpha, s, beta, G, g);
[gradR, ~]   = riemannian_grad(r, Ak, alpha, s, beta, G, g);
grad_hist(1) = norm(gradR);

converged = false;
t_end = 0;

% For BB step
r_prev = [];
g_prev = [];
step = step0;

for t = 1:maxIter
    r_t = r;

    % Riemannian gradient at r_t
    [gradR, gradE] = riemannian_grad(r_t, Ak, alpha, s, beta, G, g);
    ng = norm(gradR);

    if verbose
        fprintf('[RO] iter %d: cost=%.6e, ||grad_R||=%.3e, step=%.3e\n', ...
            t, cost_hist(t), ng, step);
    end

    if ng < tol
        converged = true;
        t_end = t-1;
        break;
    end

    dir = -gradR;

    % Armijo backtracking line-search on the manifold (via retraction)
    f0 = cost_hist(t);
    deriv0 = dot(gradR, dir); % = -||gradR||^2

    t_ls = step;
    accepted = false;
    nls = 0;

    for ls = 1:ls_maxSteps
        nls = ls;
        r_cand = retract_sphere(r_t + t_ls * dir);
        f_cand = full_cost(r_cand, Ak, alpha, s, beta, G, g);

        if f_cand <= f0 + ls_c1 * t_ls * deriv0
            accepted = true;
            break;
        end
        t_ls = ls_tau * t_ls;
    end


    if ~accepted
        % If line-search fails, fall back to a tiny step
        t_ls = max(bbClamp(1), ls_tau^ls_maxSteps * step);
        r_cand = retract_sphere(r_t + t_ls * dir);
    end

    r_new = r_cand;

    % Optional sign check w.r.t. original objective
    if checkSign
        if full_cost(-r_new, Ak, alpha, s, beta, G, g) < ...
           full_cost( r_new, Ak, alpha, s, beta, G, g)
            r_new = -r_new;
        end
    end

    step_hist(t) = norm(r_new - r_t);
    ls_hist(t) = nls;

    % Update history
    r = r_new;
    cost_hist(t+1) = full_cost(r, Ak, alpha, s, beta, G, g);
    [gradR_new, ~] = riemannian_grad(r, Ak, alpha, s, beta, G, g);

    grad_hist(t+1) = norm(gradR_new);

    % Barzilai-Borwein step-size (in ambient space, clamped)
    if useBB
        if ~isempty(r_prev)
            s_vec = r - r_prev;
            y_vec = gradE - g_prev;
            sty = dot(s_vec, y_vec);
            if abs(sty) > 1e-14
                stepBB = dot(s_vec, s_vec) / sty;
                % keep positive and clamped
                if isfinite(stepBB) && stepBB > 0
                    step = min(max(stepBB, bbClamp(1)), bbClamp(2));
                else
                    step = min(max(step, bbClamp(1)), bbClamp(2));
                end
            else
                step = min(max(step, bbClamp(1)), bbClamp(2));
            end
        end
    end

    % store previous for BB
    r_prev = r_t;
    g_prev = gradE;

    t_end = t;

    if step_hist(t) < tol
        converged = true;
        break;
    end
end

r_hat = r / max(norm(r), eps);

% pack info
info = struct();
info.converged  = converged;
info.iters      = t_end;
info.cost_hist  = cost_hist(1:t_end+1);
info.grad_hist  = grad_hist(1:t_end+1);
info.step_hist  = step_hist(1:t_end);
info.ls_hist    = ls_hist(1:t_end);
info.cost_final = info.cost_hist(end);
info.K          = K;
info.alpha      = alpha;
info.beta       = beta;
info.u          = u;
info.step0      = step0;
info.useBB      = useBB;
info.bbClamp    = bbClamp;
info.ls_c1      = ls_c1;
info.ls_tau     = ls_tau;
info.ls_maxSteps= ls_maxSteps;
info.cost_true  = cost_true;

end

% ---------------- helper: retraction on S^2 ----------------
function r = retract_sphere(x)
r = x(:);
r = r / max(norm(r), eps);
end

% ---------------- helper: Riemannian gradient on S^2 ----------------
function [gradR, gradE] = riemannian_grad(r, Ak, alpha, s, beta, G, g_vec)
% Euclidean gradient (ambient R^3)
rr = r(:);
rr = rr / max(norm(rr), eps);

g = zeros(3,1);

% sum_k alpha_k (r' Ak_k r - s_k)^2
for j = 1:numel(Ak)
    q = rr.' * Ak{j} * rr;      % scalar
    a = Ak{j} * rr;            % 3x1
    g = g + 4 * alpha(j) * (q - s(j)) * a;
end

% beta ||B_bar' r - b_u||^2 => Gradient is 2*beta*(G*r - g)
g = g + 2 * beta * (G * rr - g_vec);

gradE = g;

% Riemannian gradient: projection onto tangent space at r
% gradR = (I - r r') gradE
P = eye(3) - (rr * rr.');
gradR = P * gradE;
end

% ---------------- helper: full cost ----------------
function J = full_cost(rr, Ak, alpha, s, beta, G, g)
rr = rr(:); rr = rr / max(norm(rr), eps);
J1 = 0;
for j = 1:numel(Ak)
    q = rr.' * Ak{j} * rr;
    J1 = J1 + alpha(j) * (q - s(j))^2;
end
% Quadratic form for field consistency
J2 = beta * (rr.' * G * rr - 2 * g.' * rr);
J = J1 + J2;
end

