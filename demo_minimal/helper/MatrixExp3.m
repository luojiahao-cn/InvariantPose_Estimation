function R = MatrixExp3(so3mat)
%MATRIXEXP3 Compute the matrix exponential of a 3x3 skew-symmetric matrix
%   R = MATRIXEXP3(so3mat) computes the exponential of the skew-symmetric
%   matrix so3mat, mapping from the Lie algebra so(3) to the Lie group SO(3).
%
%   Input:
%       so3mat - 3x3 skew-symmetric matrix (so(3) element)
%
%   Output:
%       R - 3x3 rotation matrix (SO(3) element)
%
%   The exponential map for SO(3) is computed using Rodrigues' formula:
%   R = I + (sin(theta)/theta)*so3mat + ((1-cos(theta))/theta^2)*so3mat^2
%   where theta = norm(omega) is the rotation angle and omega is the
%   corresponding rotation vector (axis-angle representation).

% Extract the rotation vector from the skew-symmetric matrix
omega = VecToso3(so3mat);

% Compute the angle of rotation
theta = norm(omega);

% Handle the near-zero rotation case
if theta < 1e-10
    R = eye(3);
    return;
end

% Normalize to get the rotation axis
omega_hat = omega / theta;

% Compute sin(theta) and (1-cos(theta))
sin_theta = sin(theta);
cos_theta = cos(theta);

% Rodrigues' formula
R = eye(3) + sin_theta * so3mat / theta + (1 - cos_theta) * (so3mat / theta)^2;
end
