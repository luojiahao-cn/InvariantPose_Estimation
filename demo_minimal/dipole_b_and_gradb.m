function [B, gradB] = dipole_b_and_gradb(position_vec, moment_unit, moment_mag, mu_0)
% DIPOLE_B_AND_GRADB Compute magnetic dipole field and gradient
%
% Input:
%   position_vec  - Observation point position [x, y, z] (m)
%   moment_unit   - Dipole moment unit direction [mx, my, mz] (normalized)
%   moment_mag    - Dipole moment magnitude (A m^2)
%   mu_0          - (optional) permeability of free space (T m/A), default 4pi x 10^-7
%
% Output:
%   B      - Magnetic field [Bx, By, Bz] (T)
%   gradB  - Magnetic field gradient (3x3) (T/m)
%
% Formula:
%   B = (mu_0/4pi) * [3(m.r)r/r^5 - m/r^3]
%   gradB = (3mu_0/4pi r^5)[(m.r)I + r x m + m x r - (5(m.r)/r^2)(r x r)]

if nargin < 4
    mu_0 = 4*pi*1e-7;
end

moment_vec = moment_unit * moment_mag;
r = position_vec;
r_norm = norm(r);

if r_norm < 1e-12
    error('Observation point too close to dipole');
end

m_dot_r = dot(moment_vec, r);
common_factor = mu_0/(4*pi);

inv_r3 = 1/(r_norm^3);
inv_r5 = 1/(r_norm^5);

B = common_factor * ( (3 * m_dot_r) * r * inv_r5 - moment_vec * inv_r3 );

if nargout > 1
    factor = (3 * common_factor) / (r_norm^5);
    term_5 = 5 * m_dot_r / r_norm^2;

    gradB = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            gradB(i, j) = factor * ( ...
                m_dot_r * (i == j) ...
                + r(i)*moment_vec(j) ...
                + moment_vec(i)*r(j) ...
                - term_5 * r(i)*r(j) ...
            );
        end
    end
end
end
