function [b, nabla_b, J] = dipole_magnetic_field(nu, px, k)
    mu0 = 4*pi*1e-7;
    hatma = nu(1:3);
    pa = nu(4:6);
    p = px - pa;
    hatp = p / norm(p);
    b = mu0 * k / (4*pi * norm(p)^3) * (3*(hatp*hatp') - eye(3)) * hatma;
    Z = eye(3) - 5 * (hatp * hatp');
    J = [mu0 * k / (4*pi * norm(p)^3) * (3*(hatp*hatp') - eye(3)),...
            3 * mu0 * k / (4*pi*norm(p)^4) * (hatp'*hatma*eye(3) + hatp*hatma' + Z * hatma*hatp')];
    nabla_b = J(:, 4:6);
end