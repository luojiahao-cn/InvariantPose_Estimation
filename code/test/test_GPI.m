clc,clear,close all

B = rand(3,2)*rand(3,2)';

A = 10*rand(3,3);
A = (A + A')/2;
A = A - trace(A)/3*eye(3);

X = rand(3,3);
X = (X + X')/2;
X = X - trace(X)/3*eye(3);

R = MatrixExp3(VecToso3(rand(3,1)));

mu = 1;
f = @(R) trace(B'*R + mu * R'*A*R*X);
G = @(R) B+2*mu*A*R*X;

lA = -min(eig(A));
At = A + lA*eye(3);
lX = -min(eig(X));
Xt = X + lX*eye(3);

H = @(R) B+2*mu*At*R*Xt;
% trace(R'*At*R*Xt) - lA*trace(Xt) - lX*trace(At) + 3*lA*lX
R = MatrixExp3(VecToso3(pi*rand(3,1)));
[U, ~, V] = svd(G(R));
Rg = U*diag([1,1,det(U*V')])*V';

[U, ~, V] = svd(H(R));
Rh = U*diag([1,1,det(U*V')])*V';

f(Rg) - f(R)
f(Rh) - f(R)