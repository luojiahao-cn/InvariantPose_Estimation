syms s1 s2 s3 real
syms q11 q12 q13 q21 q22 q23 q31 q32 q33 real

% Sigma = diag([s1, s2, s3]);
Sigma = diag([3,2,1]);
% Q = [q11, q12, q13; q21, q22, q23; q31, q32, q33];
Q = MatrixExp3(VecToso3(2*pi*rand(3,1)));
Lambda = diag([1,1,-1]);
p1 = trace(Lambda*Sigma - Lambda*Q*Sigma*Q)
% trace(diag([0,0,1])*Q*Sigma*Q)
% trace(diag([0,0,1])*Sigma*Q*Q)
% det(Q*Lambda*Q)