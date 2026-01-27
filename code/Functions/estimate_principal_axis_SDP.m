function [r_hat, info] = estimate_principal_axis_SDP(prob, opts)
% estimate_principal_axis_SDP
% principal-axis estimation via Semidefinite Relaxation (SDR)
% 
% Solve (||r||=1):
%   min  sum_k alpha_k ( r' A^k r - s_k )^2  +  beta || B_bar' r - b_bar' u ||^2
%
% Relaxation:
%   X = r*r', rank(X)=1 => X >= r*r' (SDR)
%
% Inputs:
%   prob : Problem data from compute_principal_axis_prob.m
%   opts : Struct for solver parameters

if nargin < 2, opts = struct(); end

% ---------------- unpack prob ----------------
M     = prob.M;
s     = prob.s;
alpha = prob.alpha;
beta  = prob.beta;
G     = prob.G;
g     = prob.g;
K     = prob.K;

% ---------------- CVX solver ----------------
fprintf('[SDP] Starting CVX solver...\n');

cvx_begin sdp quiet
    variable r(3)
    variable R(3,3) symmetric
    variable t(numel(K))

    % Objective: sum alpha_k * t_k + beta * (tr(G*R) - 2*g'*r)
    minimize( alpha(:)' * t + beta * (trace(G * R) - 2 * g' * r) )

    subject to
        trace(R) == 1;
        [1, r'; r, R] >= 0;
        
        for i = 1:numel(K)
            dev = trace(M{i} * R) - s(i);
            [1, dev; dev, t(i)] >= 0;
        end
cvx_end


% ---------------- 3. 结果提取与认证 ----------------
% 从 R 中提取 r。
[V, D] = eig(R);
[max_eig, idx] = max(diag(D));
r_hat = V(:, idx);

% 封装信息
info = struct();
info.solver_status = cvx_status;
info.optval = cvx_optval;
info.rank_X = rank(R, 1e-4); 
info.tightness = max_eig / trace(D); 
info.R = R;

% 计算目标函数值对比 (Cost Comparison)
info.cost_est = full_cost(r_hat, M, alpha, s, beta, G, g);
if ~isempty(prob.r_true)
    info.cost_true = full_cost(prob.r_true, M, alpha, s, beta, G, g);
else
    info.cost_true = NaN;
end

end

% ---------------- helper: full cost ----------------
function J = full_cost(rr, M, alpha, s, beta, G, g)
rr = rr(:); rr = rr / max(norm(rr), eps);
J1 = 0;
for i = 1:numel(M)
    q = rr' * M{i} * rr;
    J1 = J1 + alpha(i) * (q - s(i))^2;
end
J2 = beta * (rr' * G * rr - 2 * g' * rr);
J = J1 + J2;
end
