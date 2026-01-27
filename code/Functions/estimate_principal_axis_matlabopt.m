function [r_hat, info] = estimate_principal_axis_matlabopt(b_bar, B_bar, A_p, X, opts)
% estimate_principal_axis_matlabopt
% Build objective for directed principal axis r on S^2 and solve using MATLAB optimizers.
%
% IMPORTANT (your convention):
%   b_bar : 3xm measured projected/processed matrix  \bar{B}      (MEAS)
%   B_bar : 3xm predicted matrix from p-hat          \bar{B}(\hat p) (PRED)
%
% We parameterize r on S^2 by two angles:
%   r(theta,phi) = [cos(theta)*cos(phi); cos(theta)*sin(phi); sin(theta)]
% with theta in [-pi/2, pi/2], phi in (-pi, pi].
%
% Residuals (least squares):
%   Moments:  sqrt(alpha_k) * ( r' A^k r - s_k ),   s_k = u' X^k u
%   Field:    sqrt(beta) * ( B_bar' r - b_bar' u )   (PRED' * r - MEAS' * u)
%
% Inputs:
%   b_bar : 3xm measured matrix \bar{B} (MEAS)
%   B_bar : 3xm predicted matrix \bar{B}(\hat p) (PRED)
%   A_p   : 3x3 symmetric world gradient tensor A(\hat p)
%   X     : 3x3 symmetric local tensor estimate \hat X
%   opts  : struct options:
%       .u          : 3x1 local principal axis (default [1;0;0])
%       .K          : vector of powers (default [1 2])
%       .alpha      : weights for each k (default ones)
%       .beta       : field term weight (default 1)
%       .method     : 'lsqnonlin' (default) or 'fminunc'
%       .theta0phi0 : 2x1 initial angles [theta; phi]; if empty auto init
%       .r0         : 3x1 initial r; used if theta0phi0 empty (optional)
%       .maxIter    : default 200
%       .tolFun     : default 1e-12
%       .tolX       : default 1e-10
%       .verbose    : 0/1 (default 0)
%
% Outputs:
%   r_hat : 3x1 directed unit vector in world frame
%   info  : solver outputs + final cost + residuals

if nargin < 5, opts = struct(); end
if ~isfield(opts,'u') || isempty(opts.u), opts.u = [1;0;0]; end
u = opts.u(:); u = u / max(norm(u), eps);

if ~isfield(opts,'K') || isempty(opts.K), opts.K = [1 2]; end
K = opts.K(:)';

if ~isfield(opts,'alpha') || isempty(opts.alpha), opts.alpha = ones(size(K)); end
alpha = opts.alpha(:)';
if numel(alpha) ~= numel(K), error('opts.alpha must match opts.K length.'); end

if ~isfield(opts,'beta') || isempty(opts.beta), opts.beta = 1; end
beta = opts.beta;

if ~isfield(opts,'method') || isempty(opts.method), opts.method = 'lsqnonlin'; end
method = lower(opts.method);

if ~isfield(opts,'maxIter') || isempty(opts.maxIter), opts.maxIter = 200; end
if ~isfield(opts,'tolFun') || isempty(opts.tolFun), opts.tolFun = 1e-12; end
if ~isfield(opts,'tolX') || isempty(opts.tolX), opts.tolX = 1e-10; end
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = 0; end
verbose = opts.verbose;

% Symmetrize for stability
A_p = 0.5*(A_p + A_p');
X   = 0.5*(X   + X');

% Precompute A^k and s_k = u' X^k u
kmax = max(K);
A_pow = eye(3);
X_pow = eye(3);
A_pows = cell(1,kmax);
X_pows = cell(1,kmax);
for kk = 1:kmax
    A_pow = A_pow * A_p; A_pows{kk} = 0.5*(A_pow + A_pow');
    X_pow = X_pow * X;   X_pows{kk} = 0.5*(X_pow + X_pow');
end

M = cell(1,numel(K));
s = zeros(numel(K),1);
for i = 1:numel(K)
    ki = K(i);
    M{i} = A_pows{ki};
    s(i) = u.' * X_pows{ki} * u;
end

% Field residual pieces (IMPORTANT: use MEAS for u term, PRED for r term)
b_u = b_bar.' * u;   % mx1  (MEAS' * u)

% ---- parameterization helpers ----
    function r = angles_to_r(z)
        theta = z(1); phi = z(2);
        r = [cos(theta)*cos(phi);
             cos(theta)*sin(phi);
             sin(theta)];
        r = r / max(norm(r), eps);
    end

    function z = r_to_angles(r)
        r = r(:); r = r / max(norm(r), eps);
        theta = asin(max(-1,min(1,r(3))));
        phi = atan2(r(2), r(1));
        z = [theta; phi];
    end

% ---- residual function for lsqnonlin ----
    function res = residuals(z)
        r = angles_to_r(z);

        % Moment residuals
        res_m = zeros(numel(K),1);
        for j = 1:numel(K)
            res_m(j) = sqrt(alpha(j)) * (r.' * M{j} * r - s(j));
        end

        % Field residuals (vector): PRED' * r - MEAS' * u
        res_f = sqrt(beta) * (B_bar.' * r - b_u);

        res = [res_m; res_f];
    end

% ---- scalar objective for fminunc ----
    function J = objective(z)
        rr = residuals(z);
        J = rr.'*rr;
    end

% ---- initialization ----
if isfield(opts,'theta0phi0') && ~isempty(opts.theta0phi0)
    z0 = opts.theta0phi0(:);
else
    if isfield(opts,'r0') && ~isempty(opts.r0)
        r0 = opts.r0(:);
        r0 = r0 / max(norm(r0), eps);
    else
        % Default init: dominant eigenvector of A_p
        [V,D] = eig(A_p);
        [~,idx] = max(abs(diag(D)));
        r0 = V(:,idx);
        r0 = r0 / max(norm(r0), eps);
    end
    z0 = r_to_angles(r0);
end

% ---- solve ----
switch method
    case 'lsqnonlin'
        if ~license('test','optimization_toolbox')
            error('Optimization Toolbox required for lsqnonlin.');
        end
        options = optimoptions('lsqnonlin', ...
            'Display', ternary(verbose,'iter','off'), ...
            'MaxIterations', opts.maxIter, ...
            'FunctionTolerance', opts.tolFun, ...
            'StepTolerance', opts.tolX);
        
        

        [z_star, resnorm, resvec, exitflag, output] = lsqnonlin(@residuals, z0, [], [], options);
        r_hat = angles_to_r(z_star);

        info = struct();
        info.method   = 'lsqnonlin';
        info.exitflag = exitflag;
        info.output   = output;
        info.z0       = z0;
        info.z_star   = z_star;
        info.resnorm  = resnorm;
        info.resvec   = resvec;
        info.cost     = resnorm;

    case 'fminunc'
        options = optimoptions('fminunc', ...
            'Display', ternary(verbose,'iter','off'), ...
            'MaxIterations', opts.maxIter, ...
            'OptimalityTolerance', opts.tolFun, ...
            'StepTolerance', opts.tolX);

        [z_star, fval, exitflag, output] = fminunc(@objective, z0, options);
        r_hat = angles_to_r(z_star);

        info = struct();
        info.method   = 'fminunc';
        info.exitflag = exitflag;
        info.output   = output;
        info.z0       = z0;
        info.z_star   = z_star;
        info.cost     = fval;
        info.resvec   = residuals(z_star);

    otherwise
        error('Unknown opts.method. Use ''lsqnonlin'' or ''fminunc''.');
end

info.K     = K;
info.alpha = alpha;
info.beta  = beta;
info.u     = u;

end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
