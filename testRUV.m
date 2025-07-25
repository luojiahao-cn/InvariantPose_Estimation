clc,clear,close all
A = [
  -14.0833
    7.2774
   -8.5771
    7.2774
    3.9416
    9.2442
   -8.5771
    9.2442
   10.1416
];
A = reshape(A, 3, 3);
R_true = MatrixExp3(VecToso3(rand(3,1)));
X = R_true'*A*R_true;
C = rand(3,2);
B = R_true'*C;

% 添加扰动
A_ = A + 2 * rand * gallery('orthog',3,1);
X_ = X + 2 * rand * gallery('orthog',3,1);
B_ = B + 0.1 * rand(size(B));

R0 = estimateRTC(B_, C);
R1 = estmateRTAR(A_, X_);

R_opt = estmateRTAR(A_, X_);
R_true
for i = 1:10
    R_corr = estimateRTC(B_, R_opt'*C) % 在这边迭代的话就变成R0
    R_opt = R_opt * R_corr;
    % R_corr = estmateRTAR(R_opt'*A_*R_opt, X_)
    % R_opt = R_opt * R_corr;
end

R_opt

% [U, S, V] = svd(M);
% R_opt = U*V'

function R_opt = estimateRTC(B, C)
    M = B*C';
    [U, S, V] = svd(M);
    R_opt = V*diag([1,1,det(V*U')])*U';
end

function R_opt = estmateRTAR(A, X)
    [U, LA] = eig(A);
    [V, LX] = eigs(X);
    [~, Ind] = sort(diag(LA), 'ascend');
    U = U(:,Ind);
    [~, IndX] = sort(diag(LX), 'ascend');
    V = V(:,IndX);
    R_opt = U*diag(sign(diag(V'*U)))*V';
end