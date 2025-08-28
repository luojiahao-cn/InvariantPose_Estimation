function diagnose_magnetic_field()
% Diagnose magnetic field values to understand why contour lines are not visible

close all;
clear;
clc;

% Same parameters as main function
magnet_positions = [
    -2.0,  0.0,  0.0;   % Left magnet
     2.0,  0.0,  0.0;   % Right magnet
     0.0,  0.0,  2.5    % Top magnet
];

magnet_orientations = [
    cosd(170),  sind(170),  0;    % Left magnet: rotated 170 degrees
     1,  0,  0;                   % Right magnet: pointing right
     0,  1,  0                    % Top magnet: pointing up
];

magnetic_moment = 15.0;
magnet_radius = 0.3;

fprintf('=== Magnetic Field Diagnosis ===\n');

% Test calculation on a small grid first
[X, Y] = meshgrid(linspace(-4, 4, 20), linspace(-3, 3, 20));
Z = zeros(size(X));

B_magnitude = zeros(size(X));
valid_points = 0;
total_points = numel(X);

fprintf('Calculating field on %d test points...\n', total_points);

for i = 1:numel(X)
    point = [X(i), Y(i), Z(i)];
    
    % Check if inside any magnet
    inside_magnet = false;
    for j = 1:3
        if norm(point - magnet_positions(j,:)) < magnet_radius * 1.1
            inside_magnet = true;
            break;
        end
    end
    
    if ~inside_magnet
        B_field = calculate_magnetic_field_debug(point, magnet_positions, magnet_orientations, magnetic_moment);
        B_magnitude(i) = norm(B_field);
        valid_points = valid_points + 1;
    else
        B_magnitude(i) = NaN;
    end
end

% Remove NaN values for statistics
valid_B = B_magnitude(~isnan(B_magnitude));

fprintf('Valid calculation points: %d/%d\n', valid_points, total_points);
fprintf('Magnetic field statistics:\n');
fprintf('  Min: %.2e T\n', min(valid_B));
fprintf('  Max: %.2e T\n', max(valid_B));
fprintf('  Mean: %.2e T\n', mean(valid_B));
fprintf('  Median: %.2e T\n', median(valid_B));

% Suggest contour levels
suggested_levels = logspace(log10(min(valid_B)), log10(max(valid_B)), 8);
fprintf('Suggested contour levels:\n');
for i = 1:length(suggested_levels)
    fprintf('  %.2e\n', suggested_levels(i));
end

% Create a simple 2D plot to visualize the field
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1);
imagesc(linspace(-4, 4, 20), linspace(-3, 3, 20), log10(B_magnitude'));
colorbar;
title('Log10 of Magnetic Field Magnitude');
xlabel('X Position (m)');
ylabel('Y Position (m)');

subplot(2,1,2);
contour(X, Y, B_magnitude, suggested_levels, 'ShowText', 'on', 'LineWidth', 2);
title('Contour Lines with Suggested Levels');
xlabel('X Position (m)');
ylabel('Y Position (m)');
grid on;

% Test with the actual contour levels used in main function
test_levels = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2];
fprintf('\nTesting with original contour levels:\n');
for i = 1:length(test_levels)
    count = sum(valid_B >= test_levels(i));
    fprintf('  Level %.1e: %d points above this level\n', test_levels(i), count);
end

end

function B_field = calculate_magnetic_field_debug(point, magnet_positions, magnet_orientations, magnetic_moment)
% Same calculation as main function but with debug info

B_field = [0, 0, 0];
mu_0_over_4pi = 1e-7;

for i = 1:size(magnet_positions, 1)
    magnet_pos = magnet_positions(i, :);
    magnet_orient = magnet_orientations(i, :);
    
    r_vec = point - magnet_pos;
    r_magnitude = norm(r_vec);
    
    if r_magnitude < 0.1
        continue;
    end
    
    r_unit = r_vec / r_magnitude;
    m_dot_r = dot(magnet_orient, r_unit);
    
    % Enhanced field strength for horizontal configuration
    B_dipole = (mu_0_over_4pi * magnetic_moment * 800) / r_magnitude^3 * ...
               (3 * m_dot_r * r_unit - magnet_orient);
    
    B_field = B_field + B_dipole;
end

end
