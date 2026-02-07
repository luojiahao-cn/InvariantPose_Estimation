function [r_hat, info] = estimate_principal_axis_SSL(prob, opts)
% estimate_principal_axis_SSL
% Principal-axis estimation on S^2 via Successive Semidefinite Lifting (SSL).
% 
% Solve (||r||=1):
%   min f(r) = sum_k alpha_k ( r' A^k r - s_k )^2  +  beta || B_bar' r - b_bar' u ||^2
%
% ---------------- solver options ----------------
if nargin < 2, opts = struct(); end
if ~isfield(opts,'maxIter') || isempty(opts.maxIter), opts.maxIter = 500; end
if ~isfield(opts,'tol') || isempty(opts.tol), opts.tol = 1e-8; end
if ~isfield(opts,'proxBeta') || isempty(opts.proxBeta), opts.proxBeta = 1e2; end % beta_prox
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = 0; end

% ---------------- unpack prob ----------------
Ak    = prob.Ak;
s     = prob.s;
alpha = prob.alpha;
beta  = prob.beta;
G     = prob.G;
g     = prob.g;
r     = prob.r0; % Initial estimate r_j
proxBeta = opts.proxBeta;

% ---------------- SSL iterations ----------------
maxIter = opts.maxIter;
cost_hist = zeros(maxIter+1,1);
step_hist = zeros(maxIter,1);

% Initial cost
cost_hist(1) = full_cost(r, Ak, alpha, s, beta, G, g);
converged = false;
t_end = 0;

for j = 1:maxIter
    r_j = r;

    % 1. 计算欧几里得梯度 grad_f(r_j)
    % grad_f = sum_k 4*alpha_k * (r' A^k r - s_k) * A^k * r + 2*beta * (G*r - g)
    grad_f = zeros(3,1);
    for k = 1:numel(Ak)
        A_k = Ak{k};
        scalar_term = r_j' * A_k * r_j - s(k);
        grad_f = grad_f + 4 * alpha(k) * scalar_term * (A_k * r_j);
    end
    grad_f = grad_f + 2 * beta * (G * r_j - g);

    % 2. 构造线性增长方向 m_j (Proximal Linearization)
    % m_j = -grad_f(r_j) + beta_prox * r_j
    m_j = -grad_f + proxBeta * r_j;

    % 3. 球面投影 (Normalization)
    % r_{j+1} = m_j / ||m_j||
    norm_m = norm(m_j);
    if norm_m < eps
        r_new = r_j; % Avoid division by zero
    else
        r_new = m_j / norm_m;
    end

    % Update and check convergence
    step = norm(r_new - r_j);
    step_hist(j) = step;
    r = r_new;
    cost_hist(j+1) = full_cost(r, Ak, alpha, s, beta, G, g);

    if opts.verbose
        fprintf('[SSL-S2] iter %d: cost=%.6e, step=%.3e\n', j, cost_hist(j+1), step);
    end

    t_end = j;
    if step < opts.tol
        converged = true;
        break;
    end
end

r_hat = r;

% pack info
info = struct();
info.converged = converged;
info.iters = t_end;
info.cost_hist = cost_hist(1:t_end+1);
info.step_hist = step_hist(1:t_end);
info.cost_final = info.cost_hist(end);
end

% ---------------- helper: full cost ----------------
function J = full_cost(rr, Ak, alpha, s, beta, G, g)
    J1 = 0;
    for j = 1:numel(Ak)
        q = rr' * Ak{j} * rr;
        J1 = J1 + alpha(j) * (q - s(j))^2;
    end
    % ||B'r - b'u||^2 simplifies to r'Gr - 2g'r + const
    J2 = beta * (rr' * G * rr - 2 * g' * rr);
    J = J1 + J2;
end