function [prob, opts] = compute_principal_axis_prob(b_bar, B_bar, A_p, X, opts)
% compute_principal_axis_prob
% Precompute objective function coefficients and standardize problem data
% for MM, RO, and SDP solvers.
%
% Objective:
%   min_{||r||=1} sum_k alpha_k ( r' A^k r - s_k )^2 + beta || B_bar' r - b_bar' u ||^2
%
% Input:
%   b_bar, B_bar, A_p, X : Measurements and Predicted matrices
%   opts : Struct with fields:
%       .u, .K, .alpha, .beta, .r0, .r_true, .maxIter, etc.
%
% Output:
%   prob : Standardized data structure
%       .Ak    : Cell array of symmetrized A^k
%       .s     : Vector of target invariants u' X^k u
%       .alpha : Weights for each power k in K
%       .G     : B_bar * B_bar'
%       .g     : B_bar * (b_bar' * u)
%       .beta  : Field consistency weight
%       .r0    : Normalized initial guess
%       .r_true: Normalized ground truth (if provided)
%   opts : Standardized options struct

if nargin < 5, opts = struct(); end

% 1. Standardize core parameters
if ~isfield(opts,'u') || isempty(opts.u), opts.u = [1;0;0]; end
u = opts.u(:); u = u / max(norm(u), eps);

if ~isfield(opts,'K') || isempty(opts.K), opts.K = [1, 2]; end
K = opts.K(:)';

if ~isfield(opts,'alpha') || isempty(opts.alpha), opts.alpha = [1, 1]; end
alpha = opts.alpha(:)';

if ~isfield(opts,'beta') || isempty(opts.beta), opts.beta = 1e3; end
beta = opts.beta;

% 2. Symmetrize base tensors
A_p = 0.5*(A_p + A_p');
X   = 0.5*(X   + X');

% 3. Precompute powers of A and benchmark invariants s_k
kmax = max(K);
A_pow = eye(3);
X_pow = eye(3);
A_pows = cell(1, kmax);
X_pows = cell(1, kmax);

for kk = 1:kmax
    A_pow = A_pow * A_p;
    X_pow = X_pow * X;
    A_pows{kk} = 0.5*(A_pow + A_pow');
    X_pows{kk} = 0.5*(X_pow + X_pow');
end

Ak = cell(1, numel(K));
s = zeros(1, numel(K));
for i = 1:numel(K)
    ki = K(i);
    Ak{i} = A_pows{ki};
    s(i) = u' * X_pows{ki} * u;
end

% 4. Precompute field terms G and g
G = B_bar * B_bar';
b_u = b_bar' * u;
g = B_bar * b_u;

% 5. Standardize initialization and ground truth
if isfield(opts,'r0') && ~isempty(opts.r0)
    r0 = opts.r0(:);
    r0 = r0 / max(norm(r0), eps);
else
    % Stable default: eigenvector of A_p with largest |eigenvalue|
    [V, D] = eig(A_p);
    [~, idx] = max(abs(diag(D)));
    r0 = V(:, idx);
    r0 = r0 / max(norm(r0), eps);
end

if isfield(opts,'r_true') && ~isempty(opts.r_true)
    r_true = opts.r_true(:);
    r_true = r_true / max(norm(r_true), eps);
else
    r_true = [];
end

% 6. Pack into prob structure
prob.Ak     = Ak;
prob.s      = s;
prob.K      = K;
prob.alpha  = alpha;
prob.G      = G;
prob.g      = g;
prob.beta   = beta;
prob.r0     = r0;
prob.r_true = r_true;
prob.u      = u; % Preserve for reference

end
