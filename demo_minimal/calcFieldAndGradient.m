%% Calculate magnetic field and gradient tensor
function [b_p, A_p] = calcFieldAndGradient(p, m_pos, m_hat, m_norm)
    b_p = zeros(3, 1);
    A_p = zeros(3, 3);
    for i = 1:size(m_pos,2)
        r = p - m_pos(:, i);
        [B, gradB] = dipole_b_and_gradb(r, m_hat(:, i), m_norm(i));
        b_p = b_p + B;
        A_p = A_p + gradB;
    end
end
