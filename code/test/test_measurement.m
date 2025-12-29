clc,clear,close all

m_pos = [0; 0; 0];
m_hat = [0; 0; 1];
m_norm = 300e-4;

pos = [0.01; 0; 0];
theta = [0; 0; 0];

R = MatrixExp3(VecToso3(theta));

[B_t, gradB_t] = calcFieldAndGradient(pos, m_pos, m_hat, m_norm);

% measurement
b = R * B_t
gradb = R' * gradB_t * R