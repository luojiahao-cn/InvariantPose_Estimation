function test_points = generate_test_points(params, num_points)
% GENERATE_TEST_POINTS Generate multiple test points
% Input:
%   params     - Experiment parameter structure
%   num_points - Number of test points
% Output:
%   test_points - Struct array with:
%                 - p_true: 3x1 true position
%                 - theta_true: 3x1 true rotation vector

workspace_center = params.workspace.center;
workspace_radius = params.workspace.radius - 0.01;

test_points = struct('p_true', {}, 'theta_true', {});

i = 1;
while i <= num_points
    r = workspace_radius * (rand()^(1/3));
    theta = 2*pi*rand();
    phi = acos(rand());
    p_true_candidate = workspace_center + r * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)];
    if p_true_candidate(3) >= workspace_center(3)
        p_true = p_true_candidate;
        theta_true = pi * (2*rand(3,1) - 1);
        test_points(i).p_true = p_true;
        test_points(i).theta_true = theta_true;
        i = i + 1;
    end
end
end
