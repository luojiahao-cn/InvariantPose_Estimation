function so3mat = VecToso3(omega)
%VECTSO3 Convert a rotation vector to a 3x3 skew-symmetric matrix
%   so3mat = VECTSO3(omega) converts a 3x1 rotation vector omega to the
%   corresponding 3x3 skew-symmetric matrix so3mat in the Lie algebra so(3).
%
%   Input:
%       omega - 3x1 rotation vector (axis-angle representation)
%
%   Output:
%       so3mat - 3x3 skew-symmetric matrix
%
%   The mapping is:
%   omega = [omega_x; omega_y; omega_z] ->
%   so3mat = [  0  -omega_z  omega_y;
%             omega_z    0   -omega_x;
%            -omega_y  omega_x     0  ]

so3mat = [  0        -omega(3)   omega(2);
         omega(3)      0         -omega(1);
        -omega(2)   omega(1)       0      ];
end
