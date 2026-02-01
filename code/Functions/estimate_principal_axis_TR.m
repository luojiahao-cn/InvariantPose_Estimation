function [r_hat, info] = estimate_principal_axis_TR(prob, opts)
% estimate_principal_axis_TR
% Trust-region method on S^2 for principal-axis estimation.
%
% Same objective as estimate_principal_axis_RO:
%   min_{||r||=1}  sum_j alpha_j (r' A_j r - s_j)^2  +  beta*(r'G r - 2 g'r)

if nargin < 2, opts = struct(); end

% ---------- options ----------
if ~isfield(opts,'maxIter') || isempty(opts,'maxIter'), opts.maxIter = 200; end
if ~isfield(opts,'tol')     || isempty(opts,'tol'),     opts.tol = 1e-8; end
if ~isfield(opts,'verbose') || isempty(opts,'verbose'), opts.verbose = 0; end

% Trust-region parameters
if ~isfield(opts,'Delta0')  || isempty(opts,'Delta0'),  opts.Delta0 = 0.5; end
if ~isfield(opts,'DeltaMax')|| isempty(opts,'DeltaMax'),opts.DeltaMax = 2.0; end
if ~isfield(opts,'eta_accept')|| isempty(opts,'eta_accept'), opts.eta_accept = 0.1; end
if ~isfield(opts,'eta_good')|| isempty(opts,'eta_good'), opts.eta_good = 0.75; end
if ~isfield(opts,'gamma_dec')|| isempty(opts,'gamma_dec'), opts.gamma_dec = 0.25; end
if ~isfield(opts,'gamma_inc')|| isempty(opts,'gamma_inc'), opts.gamma_inc = 2.0; end

% TCG (Steihaug) parameters
if ~isfield(opts,'tcg_maxIter') || isempty(opts,'tcg_maxIter'), opts.tcg_maxIter = 50; end
if ~isfield(opts,'tcg_tol')     || isempty(opts,'tcg_tol'),     opts.tcg_tol = 1e-10; end

maxIter = opts.maxIter;
tol     = opts.tol;
verbose = opts.verbose;

Delta   = opts.Delta0;
DeltaMax= opts.DeltaMax;

% ---------- unpack prob ----------
Ak    = prob.Ak;
s     = prob.s;
alpha = prob.alpha;
beta  = prob.beta;
G     = prob.G;
g     = prob.g;
r     = prob.r0;

% ---------- history ----------
cost_hist = zeros(maxIter+1,1);
grad_hist = zeros(maxIter+1,1);
Delta_hist= zeros(maxIter+1,1);
rho_hist  = zeros(maxIter,1);
tcg_hist  = zeros(maxIter,1);

cost_hist(1) = full_cost(r, Ak, alpha, s, beta, G, g);
[gradR, gradE] = riemannian_grad(r, Ak, alpha, s, beta, G, g);
grad_hist(1) = norm(gradR);
Delta_hist(1)= Delta;

converged = false;
t_end = 0;

for t = 1:maxIter
    rr = r(:); rr = rr / max(norm(rr), eps);

    [gradR, gradE] = riemannian_grad(rr, Ak, alpha, s, beta, G, g);
    ng = norm(gradR);

    if verbose
        fprintf('[TR] iter %d: cost=%.6e, ||grad||=%.3e, Delta=%.3e\n', ...
            t, cost_hist(t), ng, Delta);
    end

    if ng < tol
        converged = true;
        t_end = t-1;
        break;
    end

    % ---- Solve TR subproblem via Steihaug-TCG in tangent space ----
    Hv = @(v) hessR_times(rr, v, Ak, alpha, s, beta, G, gradE); % Riemannian Hv
    [eta, n_cg] = steihaug_tcg(gradR, Hv, rr, Delta, opts.tcg_maxIter, opts.tcg_tol);
    tcg_hist(t) = n_cg;

    % candidate
    r_cand = retract_sphere(rr + eta);
    f0 = cost_hist(t);
    f_cand = full_cost(r_cand, Ak, alpha, s, beta, G, g);

    % predicted decrease: m(0)-m(eta) = - <g,eta> - 1/2 <eta, Heta>
    Heta = Hv(eta);
    pred = - dot(gradR, eta) - 0.5 * dot(eta, Heta);
    if pred <= 0
        % model not predictive, force rejection and shrink
        rho = -inf;
    else
        rho = (f0 - f_cand) / pred;
    end
    rho_hist(t) = rho;

    % ---- TR radius update ----
    if rho < opts.eta_accept
        % reject
        Delta = max(opts.gamma_dec * Delta, 1e-12);
        r_new = rr;
        f_new = f0;
    else
        % accept
        r_new = r_cand;
        f_new = f_cand;

        if rho > opts.eta_good && abs(norm(eta) - Delta) < 1e-8
            Delta = min(opts.gamma_inc * Delta, DeltaMax);
        end
    end

    r = r_new;
    cost_hist(t+1) = f_new;
    [gradR_new, ~] = riemannian_grad(r, Ak, alpha, s, beta, G, g);
    grad_hist(t+1) = norm(gradR_new);
    Delta_hist(t+1)= Delta;

    t_end = t;

    % extra stop: tiny update
    if norm(r_new - rr) < tol
        converged = true;
        break;
    end
end

r_hat = r(:) / max(norm(r), eps);

info = struct();
info.converged = converged;
info.iters = t_end;
info.cost_hist = cost_hist(1:t_end+1);
info.grad_hist = grad_hist(1:t_end+1);
info.Delta_hist= Delta_hist(1:t_end+1);
info.rho_hist  = rho_hist(1:t_end);
info.tcg_hist  = tcg_hist(1:t_end);
info.cost_final= info.cost_hist(end);
end

% ---------- helper: retraction ----------
function r = retract_sphere(x)
r = x(:);
r = r / max(norm(r), eps);
end

% ---------- helper: Riemannian gradient (reuse your version) ----------
function [gradR, gradE] = riemannian_grad(r, Ak, alpha, s, beta, G, g_vec)
rr = r(:);
rr = rr / max(norm(rr), eps);

g = zeros(3,1);
for j = 1:numel(Ak)
    q = rr.' * Ak{j} * rr;
    a = Ak{j} * rr;
    g = g + 4 * alpha(j) * (q - s(j)) * a;
end
g = g + 2 * beta * (G * rr - g_vec);

gradE = g;
P = eye(3) - (rr * rr.');
gradR = P * gradE;
end

% ---------- Euclidean Hessian-vector product ----------
function HvE = hessE_times(r, v, Ak, alpha, s, beta, G)
rr = r(:); rr = rr / max(norm(rr), eps);
vv = v(:);

HvE = zeros(3,1);

for j = 1:numel(Ak)
    Aj = Ak{j};
    q  = rr.' * Aj * rr;      % scalar
    a  = Aj * rr;             % Aj*r
    av = a.' * vv;            % (Aj*r)^T v
    HvE = HvE + 4*alpha(j) * ( 2*av*a + (q - s(j))*(Aj*vv) );
end

HvE = HvE + 2*beta*(G*vv);
end

% ---------- Riemannian Hessian-vector product on S^2 ----------
function HvR = hessR_times(r, v, Ak, alpha, s, beta, G, gradE)
rr = r(:); rr = rr / max(norm(rr), eps);
vv = v(:);

P = eye(3) - rr*rr.';
HvE = hessE_times(rr, vv, Ak, alpha, s, beta, G);

% Hess_S2[v] = P(HvE) - <r, gradE> v
HvR = P*HvE - (rr.'*gradE)*vv;

% ensure tangent (numerical)
HvR = P*HvR;
end

% ---------- Steihaug truncated CG for TR subproblem ----------
function [eta, iters] = steihaug_tcg(gradR, Hv, r, Delta, maxIter, tol)
% Solve: min <g,eta> + 1/2 <eta, H eta>  s.t. eta ⟂ r, ||eta|| <= Delta
% g = gradR, eta in tangent already if initialized as 0.

eta = zeros(3,1);
g = gradR(:);

% Residual for quadratic model: r0 = g + H*eta = g
res = g;
p = -res;

% Ensure p is tangent
P = eye(3) - r(:)*r(:).';
p = P*p;

rsold = dot(res,res);
iters = 0;

if sqrt(rsold) < tol
    return;
end

for k = 1:maxIter
    iters = k;
    Hp = Hv(p);
    pHp = dot(p, Hp);

    if pHp <= 0
        % negative curvature: go to boundary along p
        tau = tau_to_boundary(eta, p, Delta);
        eta = eta + tau*p;
        return;
    end

    alpha = rsold / pHp;
    eta_next = eta + alpha*p;

    if norm(eta_next) >= Delta
        % hit boundary
        tau = tau_to_boundary(eta, p, Delta);
        eta = eta + tau*p;
        return;
    end

    eta = eta_next;

    res = res + alpha*Hp;
    rsnew = dot(res,res);

    if sqrt(rsnew) < tol
        return;
    end

    beta_cg = rsnew / rsold;
    p = -res + beta_cg*p;
    p = P*p; % keep tangent

    rsold = rsnew;
end
end

function tau = tau_to_boundary(eta, p, Delta)
% Solve ||eta + tau p|| = Delta for tau>0
a = dot(p,p);
b = 2*dot(eta,p);
c = dot(eta,eta) - Delta^2;

% tau = (-b + sqrt(b^2 - 4ac)) / (2a)
disc = max(b*b - 4*a*c, 0);
tau = (-b + sqrt(disc)) / (2*a);
end

% ---------- full cost (reuse your version) ----------
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
