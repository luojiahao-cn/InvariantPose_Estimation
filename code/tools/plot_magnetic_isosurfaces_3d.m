function plot_magnetic_isosurfaces_3d()
% PLOT_MAGNETIC_ISOSURFACES_3D - Visualize three permanent magnets and their magnetic field contour lines in 3D
% 
% This function creates a 3D visualization of three permanent magnets (left, right, and top)
% with their magnetic field contour lines (constant magnetic field strength lines).
% All magnets are configured with horizontal plane orientations for better visibility.
%
% Author: Generated for research publication
% Date: 2024

close all;
clear;
clc;

% Define magnet parameters
magnet_radius = 0.3;    % Radius of cylindrical magnets
magnet_height = 0.6;    % Height of cylindrical magnets
magnetic_moment = 15.0;  % Enhanced magnetic dipole moment strength

% Define magnet positions (left, right, top)
magnet_positions = [
    -2.0,  0.0,  0.0;   % Left magnet
     2.0,  0.0,  0.0;   % Right magnet
     0.0,  0.0,  2.5    % Top magnet
];

% Define magnet orientations (ALL in horizontal plane for maximum visibility)
magnet_orientations = [
    cosd(170),  sind(170),  0;    % Left magnet: rotated 170 degrees in XY plane
     1,  0,  0;                   % Right magnet: north pole pointing right (+X)
     0,  1,  0                    % Top magnet: north pole pointing up (+Y)
];

% Display magnet configuration
fprintf('=== Magnetic Field Contour Visualization ===\n');
fprintf('All magnets configured in horizontal plane for optimal visibility\n');
fprintf('Left magnet: Position (%.1f, %.1f, %.1f), Orientation (%.3f, %.3f, %.3f), Rotation: 170°\n', ...
        magnet_positions(1,:), magnet_orientations(1,:));
fprintf('Right magnet: Position (%.1f, %.1f, %.1f), Orientation (%.3f, %.3f, %.3f)\n', ...
        magnet_positions(2,:), magnet_orientations(2,:));
fprintf('Top magnet: Position (%.1f, %.1f, %.1f), Orientation (%.3f, %.3f, %.3f)\n', ...
        magnet_positions(3,:), magnet_orientations(3,:));

% Create figure
figure('Position', [100, 100, 1200, 900]);
hold on;
axis equal;
grid on;
view(45, 30);

% Set up 3D space
xlim([-4, 4]);
ylim([-3, 3]);
zlim([-1, 4]);

% Draw the three magnets as cylinders
fprintf('Drawing magnets...\n');
for i = 1:3
    pos = magnet_positions(i, :);
    orientation = magnet_orientations(i, :);
    
    % Create cylinder coordinates
    [theta, z] = meshgrid(linspace(0, 2*pi, 20), linspace(-magnet_height/2, magnet_height/2, 10));
    x_cyl = magnet_radius * cos(theta);
    y_cyl = magnet_radius * sin(theta);
    z_cyl = z;
    
    % Rotate cylinder for XY plane orientations
    if abs(orientation(3)) < 0.5  % XY plane orientation
        angle = atan2(orientation(2), orientation(1));
        cos_a = cos(angle);
        sin_a = sin(angle);
        
        % Rotate the cylinder
        x_rot = x_cyl * cos_a - y_cyl * sin_a;
        y_rot = x_cyl * sin_a + y_cyl * cos_a;
        z_rot = z_cyl;
    else
        x_rot = x_cyl;
        y_rot = y_cyl;
        z_rot = z_cyl;
    end
    
    % Translate to magnet position
    x_final = x_rot + pos(1);
    y_final = y_rot + pos(2);
    z_final = z_rot + pos(3);
    
    % Create color map for poles
    color_map = zeros(size(z_cyl, 1), size(z_cyl, 2), 3);
    
    % Determine pole orientation for horizontal plane magnets
    north_mask = z > 0;  % Simple pole assignment
    
    % Assign colors: red for north, blue for south
    color_map(:, :, 1) = north_mask;        % Red channel
    color_map(:, :, 3) = ~north_mask;       % Blue channel
    
    % Plot the cylinder
    surf(x_final, y_final, z_final, color_map, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
end

% Add magnet labels
text(magnet_positions(1,1)-0.7, magnet_positions(1,2)-0.5, magnet_positions(1,3)+0.8, ...
     'Left (170°)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'red');
text(magnet_positions(2,1)+0.7, magnet_positions(2,2)-0.5, magnet_positions(2,3)+0.8, ...
     'Right', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'blue');
text(magnet_positions(3,1), magnet_positions(3,2)+0.8, magnet_positions(3,3)+0.5, ...
     'Top', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'green');

% Generate magnetic field contour lines
fprintf('Calculating magnetic field contours...\n');

% Multiple horizontal planes for better visualization
planes = [
    struct('z', 0, 'color', [0.0, 0.0, 1.0], 'name', 'Z=0'),
    struct('z', 0.5, 'color', [0.0, 0.6, 1.0], 'name', 'Z=0.5'),
    struct('z', -0.5, 'color', [0.6, 0.0, 1.0], 'name', 'Z=-0.5'),
    struct('z', 1.0, 'color', [0.0, 1.0, 0.6], 'name', 'Z=1.0')
];

% Optimized contour levels for horizontal configuration
contour_levels = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2];

for p = 1:length(planes)
    plane = planes(p);
    
    % Create 2D grid for this plane
    [X, Y] = meshgrid(linspace(-4, 4, 100), linspace(-3, 3, 100));
    Z = ones(size(X)) * plane.z;
    
    % Calculate magnetic field magnitude
    B_magnitude = zeros(size(X));
    
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
            B_field = calculate_magnetic_field(point, magnet_positions, magnet_orientations, magnetic_moment);
            B_magnitude(i) = norm(B_field);
        else
            B_magnitude(i) = NaN;
        end
    end
    
    % Draw enhanced contour lines
    [C, h] = contour3(X, Y, Z, B_magnitude, contour_levels, ...
                      'LineWidth', 3.0, 'LineColor', plane.color);
    
    % Add contour labels for main plane
    if p == 1
        clabel(C, h, 'FontSize', 9, 'Color', plane.color, 'FontWeight', 'bold');
    end
    
    fprintf('Completed plane %s\n', plane.name);
end

% Enhance plot appearance
xlabel('X Position (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y Position (m)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Z Position (m)', 'FontSize', 14, 'FontWeight', 'bold');
title({'3D Magnetic Field Visualization: Horizontal Plane Configuration', ...
       'Left Magnet Rotated 170° - Enhanced Contour Lines'}, ...
      'FontSize', 15, 'FontWeight', 'bold');

% Create legend
h1 = patch([0 1 1 0], [0 0 1 1], [0 0 0 0], 'red', 'EdgeColor', 'none', 'Visible', 'off');
h2 = patch([0 1 1 0], [0 0 1 1], [0 0 0 0], 'blue', 'EdgeColor', 'none', 'Visible', 'off');
h3 = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 3.0);
h4 = plot3(NaN, NaN, NaN, 'Color', [0.0, 0.6, 1.0], 'LineWidth', 3.0);

legend([h1, h2, h3, h4], {'North Pole (Red)', 'South Pole (Blue)', ...
       'Main Plane Contours (Z=0)', 'Additional Plane Contours'}, ...
       'Location', 'northeast', 'FontSize', 11);

% Set lighting
lighting gouraud;
light('Position', [5, 5, 5]);
light('Position', [-5, -5, 5]);

% Set appearance
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 12);
grid on;
set(gca, 'GridAlpha', 0.3);

fprintf('=== Visualization Complete ===\n');
fprintf('Contour levels: ');
fprintf('%.1e ', contour_levels);
fprintf('\n');

end

function B_field = calculate_magnetic_field(point, magnet_positions, magnet_orientations, magnetic_moment)
% Enhanced magnetic field calculation for horizontal plane configuration

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