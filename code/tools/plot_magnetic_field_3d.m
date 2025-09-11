function plot_magnetic_field_3d()
% PLOT_MAGNETIC_FIELD_3D - Visualize three permanent magnets and their magnetic field lines in 3D
% 
% This function creates a 3D visualization of three permanent magnets (left, right, and top)
% with their magnetic field lines. Each magnet is represented as a cylinder with 
% red and blue colors indicating north and south poles respectively.
%
% Author: Generated for research publication
% Date: 2024

close all;
clear;
clc;

% Define magnet parameters
magnet_radius = 0.3;    % Radius of cylindrical magnets
magnet_height = 0.6;    % Height of cylindrical magnets
magnetic_moment = 1.0;  % Magnetic dipole moment strength

% Nature-inspired color scheme options
north_color = [0.902, 0.294, 0.208]; % NPG red (#E64B35)
south_color = [0.235, 0.329, 0.533]; % NPG blue (#3C5488)
field_line_color = [0.439, 0.502, 0.565]; % slate gray (#708090)
coupling_line_color = field_line_color;  


% Define magnet positions (left, right, top)
magnet_positions = [
    -2.0,  0.0,  0.0;   % Left magnet
     2.0,  0.0,  0.0;   % Right magnet
     0.0,  0.0,  2.5    % Top magnet
];

% Define magnet orientations (magnetic dipole directions)
magnet_orientations = [
     -1,  0,  0;    % Left magnet: north pole pointing right (+X)
    -1,  0,  0;    % Right magnet: north pole pointing left (-X)
     0,  0, -1     % Top magnet: north pole pointing down (-Z)
];

% Create figure
figure('Position', [100, 100, 1200, 900]);
hold on;
axis equal;
grid on;
view(45, 30);

% Set up 3D space
xlim([-4, 4]);
ylim([-3, 3]);
zlim([-2, 4]);

% Draw the three magnets as cylinders
for i = 1:3
    pos = magnet_positions(i, :);
    orientation = magnet_orientations(i, :);
    
    % Create cylinder coordinates
    [theta, z] = meshgrid(linspace(0, 2*pi, 20), linspace(-magnet_height/2, magnet_height/2, 10));
    x_cyl = magnet_radius * cos(theta);
    y_cyl = magnet_radius * sin(theta);
    z_cyl = z;
    
    % Rotate cylinder based on orientation
    if abs(orientation(1)) > 0.5  % X-oriented magnet
        % Rotate around Y axis
        x_rot = z_cyl;
        y_rot = y_cyl;
        z_rot = x_cyl;
    elseif abs(orientation(3)) > 0.5  % Z-oriented magnet
        % Keep as is
        x_rot = x_cyl;
        y_rot = y_cyl;
        z_rot = z_cyl;
    else
        % Default orientation
        x_rot = x_cyl;
        y_rot = y_cyl;
        z_rot = z_cyl;
    end
    
    % Translate to magnet position
    x_final = x_rot + pos(1);
    y_final = y_rot + pos(2);
    z_final = z_rot + pos(3);
    
    % Create color map for north (red) and south (blue) poles
    color_map = zeros(size(z_cyl, 1), size(z_cyl, 2), 3);
    
    % Determine which end is north based on orientation
    if abs(orientation(1)) > 0.5  % X-oriented magnet
        if orientation(1) > 0
            % North pole at positive X end (z > 0 in local coordinates)
            north_mask = z > 0;
        else
            % North pole at negative X end (z < 0 in local coordinates)
            north_mask = z < 0;
        end
    else  % Z-oriented magnet
        if orientation(3) < 0
            % North pole at negative Z end (z < 0 in local coordinates)
            north_mask = z < 0;
        else
            % North pole at positive Z end (z > 0 in local coordinates)
            north_mask = z > 0;
        end
    end
    
    % Assign colors using Nature color scheme
    color_map(:, :, 1) = north_mask * north_color(1) + (~north_mask) * south_color(1);  % Red channel
    color_map(:, :, 2) = north_mask * north_color(2) + (~north_mask) * south_color(2);  % Green channel
    color_map(:, :, 3) = north_mask * north_color(3) + (~north_mask) * south_color(3);  % Blue channel
    
    % Plot the cylinder as solid
    surf(x_final, y_final, z_final, color_map, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    
    % Add top and bottom caps to make it solid
    if abs(orientation(1)) > 0.5  % X-oriented magnet
        % Create circular caps
        theta_cap = linspace(0, 2*pi, 20);
        r_cap = linspace(0, magnet_radius, 5);
        [THETA_cap, R_cap] = meshgrid(theta_cap, r_cap);
        
        Y_cap = R_cap .* cos(THETA_cap);
        Z_cap = R_cap .* sin(THETA_cap);
        X_cap_top = ones(size(Y_cap)) * (magnet_height/2);
        X_cap_bottom = ones(size(Y_cap)) * (-magnet_height/2);
        
        X_cap_final_top = X_cap_top + pos(1);
        Y_cap_final_top = Y_cap + pos(2);
        Z_cap_final_top = Z_cap + pos(3);
        
        X_cap_final_bottom = X_cap_bottom + pos(1);
        Y_cap_final_bottom = Y_cap + pos(2);
        Z_cap_final_bottom = Z_cap + pos(3);
        
        % Color the caps based on pole orientation using Nature colors
        if orientation(1) > 0  % North pole pointing in +X direction
            % Top cap is at +X (north pole) - Nature north color
            surf(X_cap_final_top, Y_cap_final_top, Z_cap_final_top, ...
                 repmat(reshape(north_color, 1, 1, 3), size(Y_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
            % Bottom cap is at -X (south pole) - Nature south color
            surf(X_cap_final_bottom, Y_cap_final_bottom, Z_cap_final_bottom, ...
                 repmat(reshape(south_color, 1, 1, 3), size(Y_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        else  % North pole pointing in -X direction
            % Top cap is at +X (south pole) - Nature south color
            surf(X_cap_final_top, Y_cap_final_top, Z_cap_final_top, ...
                 repmat(reshape(south_color, 1, 1, 3), size(Y_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
            % Bottom cap is at -X (north pole) - Nature north color
            surf(X_cap_final_bottom, Y_cap_final_bottom, Z_cap_final_bottom, ...
                 repmat(reshape(north_color, 1, 1, 3), size(Y_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        end
    else  % Z-oriented magnet
        % Create circular caps
        theta_cap = linspace(0, 2*pi, 20);
        r_cap = linspace(0, magnet_radius, 5);
        [THETA_cap, R_cap] = meshgrid(theta_cap, r_cap);
        
        X_cap = R_cap .* cos(THETA_cap);
        Y_cap = R_cap .* sin(THETA_cap);
        Z_cap_top = ones(size(X_cap)) * (magnet_height/2);
        Z_cap_bottom = ones(size(X_cap)) * (-magnet_height/2);
        
        X_cap_final = X_cap + pos(1);
        Y_cap_final = Y_cap + pos(2);
        Z_cap_final_top = Z_cap_top + pos(3);
        Z_cap_final_bottom = Z_cap_bottom + pos(3);
        
        % Color the caps based on pole orientation using Nature colors
        if orientation(3) < 0  % North pole pointing in -Z direction (down)
            % Top cap is at +Z (south pole) - Nature south color
            surf(X_cap_final, Y_cap_final, Z_cap_final_top, ...
                 repmat(reshape(south_color, 1, 1, 3), size(X_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
            % Bottom cap is at -Z (north pole) - Nature north color
            surf(X_cap_final, Y_cap_final, Z_cap_final_bottom, ...
                 repmat(reshape(north_color, 1, 1, 3), size(X_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        else  % North pole pointing in +Z direction (up)
            % Top cap is at +Z (north pole) - Nature north color
            surf(X_cap_final, Y_cap_final, Z_cap_final_top, ...
                 repmat(reshape(north_color, 1, 1, 3), size(X_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
            % Bottom cap is at -Z (south pole) - Nature south color
            surf(X_cap_final, Y_cap_final, Z_cap_final_bottom, ...
                 repmat(reshape(south_color, 1, 1, 3), size(X_cap)), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        end
    end
end

% Generate magnetic field lines
field_line_length = 300;  % Increased steps
step_size = 0.01;  % Smaller step size for better accuracy

% Define starting points for field lines (around magnets)
starting_points = [];

% Field lines starting from magnet poles (north and south ends)
for i = 1:3
    pos = magnet_positions(i, :);
    orientation = magnet_orientations(i, :);
    
    % Generate points mainly from the pole faces
    if abs(orientation(1)) > 0.5  % X-oriented magnet
        % Determine north and south pole positions correctly
        if orientation(1) > 0  % North pole pointing right (+X)
            north_pos = pos + [magnet_height/2, 0, 0];   % Right end is north
            south_pos = pos + [-magnet_height/2, 0, 0];  % Left end is south
        else  % North pole pointing left (-X)
            north_pos = pos + [-magnet_height/2, 0, 0];  % Left end is north
            south_pos = pos + [magnet_height/2, 0, 0];   % Right end is south
        end
        
        % Generate points from north pole
        for j = 1:8
            angle = (j-1) * 2 * pi / 8;
            r = magnet_radius * 0.8 * sqrt(rand);
            y_offset = r * cos(angle);
            z_offset = r * sin(angle);
            starting_points = [starting_points; north_pos + [0, y_offset, z_offset]];
        end
        
        % Generate points from south pole
        for j = 1:4
            angle = (j-1) * 2 * pi / 4;
            r = magnet_radius * 0.6 * sqrt(rand);
            y_offset = r * cos(angle);
            z_offset = r * sin(angle);
            starting_points = [starting_points; south_pos + [0, y_offset, z_offset]];
        end
        
    else  % Z-oriented magnet
        % Determine north and south pole positions correctly
        if orientation(3) < 0  % North pole pointing down (-Z)
            north_pos = pos + [0, 0, -magnet_height/2];  % Bottom end is north
            south_pos = pos + [0, 0, magnet_height/2];   % Top end is south
        else  % North pole pointing up (+Z)
            north_pos = pos + [0, 0, magnet_height/2];   % Top end is north
            south_pos = pos + [0, 0, -magnet_height/2];  % Bottom end is south
        end
        
        % Generate points from north pole
        for j = 1:8
            angle = (j-1) * 2 * pi / 8;
            r = magnet_radius * 0.8 * sqrt(rand);
            x_offset = r * cos(angle);
            y_offset = r * sin(angle);
            starting_points = [starting_points; north_pos + [x_offset, y_offset, 0]];
        end
        
        % Generate points from south pole
        for j = 1:4
            angle = (j-1) * 2 * pi / 4;
            r = magnet_radius * 0.6 * sqrt(rand);
            x_offset = r * cos(angle);
            y_offset = r * sin(angle);
            starting_points = [starting_points; south_pos + [x_offset, y_offset, 0]];
        end
    end
end

% Add inter-magnet coupling field lines
% Left magnet (North pointing right) to Right magnet (North pointing left)
left_north = magnet_positions(1, :) + [magnet_height/2, 0, 0];   % Left magnet's north pole (right end)
right_south = magnet_positions(2, :) + [magnet_height/2, 0, 0];  % Right magnet's south pole (right end)
for i = 1:8
    % Points between left north and right south
    t = (i-1) / 7;  % Parameter from 0 to 1
    interp_point = (1-t) * left_north + t * right_south;
    % Add some perpendicular offset
    offset_y = 0.2 * sin(i * pi / 4);
    offset_z = 0.2 * cos(i * pi / 4);
    starting_points = [starting_points; interp_point + [0, offset_y, offset_z]];
end

% Right magnet (North pointing left) to Left magnet (South pointing right)
right_north = magnet_positions(2, :) + [-magnet_height/2, 0, 0];  % Right magnet's north pole (left end)
left_south = magnet_positions(1, :) + [-magnet_height/2, 0, 0];   % Left magnet's south pole (left end)
for i = 1:6
    t = (i-1) / 5;
    interp_point = (1-t) * right_north + t * left_south;
    offset_y = 0.15 * sin(i * pi / 3);
    offset_z = 0.15 * cos(i * pi / 3);
    starting_points = [starting_points; interp_point + [0, offset_y, offset_z]];
end

% Top magnet (North pointing down) to Left and Right magnets
top_north = magnet_positions(3, :) + [0, 0, -magnet_height/2];   % Top magnet's north pole (bottom end)

% Top to Left magnet's south pole
left_south = magnet_positions(1, :) + [-magnet_height/2, 0, 0];   % Left magnet's south pole (left end)
for i = 1:6
    t = (i-1) / 5;
    interp_point = (1-t) * top_north + t * left_south;
    offset_x = 0.1 * sin(i * pi / 3);
    offset_y = 0.1 * cos(i * pi / 3);
    starting_points = [starting_points; interp_point + [offset_x, offset_y, 0]];
end

% Top to Right magnet's north pole  
right_north = magnet_positions(2, :) + [-magnet_height/2, 0, 0];  % Right magnet's north pole (left end)
for i = 1:6
    t = (i-1) / 5;
    interp_point = (1-t) * top_north + t * right_north;
    offset_x = 0.1 * sin(i * pi / 3 + pi/2);
    offset_y = 0.1 * cos(i * pi / 3 + pi/2);
    starting_points = [starting_points; interp_point + [offset_x, offset_y, 0]];
end

% Add some field lines starting from space around magnets for long-range coupling
for i = 1:15
    angle1 = rand * 2 * pi;
    angle2 = rand * pi - pi/2;
    radius = 2.0 + rand * 1.5;
    x_start = radius * cos(angle1) * cos(angle2);
    y_start = radius * sin(angle1) * cos(angle2);
    z_start = radius * sin(angle2);
    starting_points = [starting_points; x_start, y_start, z_start];
end

% Calculate and plot field lines with different colors for different types
for i = 1:size(starting_points, 1)
    start_point = starting_points(i, :);
    
    % Calculate field line in both directions for better continuity
    field_line_forward = calculate_field_line(start_point, magnet_positions, magnet_orientations, ...
                                            magnetic_moment, field_line_length, step_size);
    field_line_backward = calculate_field_line(start_point, magnet_positions, magnet_orientations, ...
                                             magnetic_moment, field_line_length, -step_size);
    
    % Combine forward and backward field lines
    if size(field_line_backward, 1) > 1
        field_line_backward = flipud(field_line_backward(2:end, :));  % Remove duplicate start point
        field_line = [field_line_backward; field_line_forward];
    else
        field_line = field_line_forward;
    end
    
    if size(field_line, 1) > 10  % Require longer field lines for plotting
        % Check if this is a coupling line (starts between magnets)
        is_coupling_line = false;
        
        % Check if starting point is in the coupling region between magnets
        for j = 1:3
            for k = j+1:3
                magnet1_pos = magnet_positions(j, :);
                magnet2_pos = magnet_positions(k, :);
                midpoint = (magnet1_pos + magnet2_pos) / 2;
                dist_to_midpoint = norm(start_point - midpoint);
                dist_between_magnets = norm(magnet1_pos - magnet2_pos);
                
                if dist_to_midpoint < dist_between_magnets * 0.4
                    is_coupling_line = true;
                    break;
                end
            end
            if is_coupling_line
                break;
            end
        end
        
        if is_coupling_line
            % Use Nature coupling line color
            plot3(field_line(:, 1), field_line(:, 2), field_line(:, 3), ...
                  'Color', coupling_line_color, 'LineWidth', 1.2, 'LineStyle', '-');
        else
            % Use Nature field line color
            plot3(field_line(:, 1), field_line(:, 2), field_line(:, 3), ...
                  'Color', field_line_color, 'LineWidth', 1.2, 'LineStyle', '-');
        end
    end
end


% Set lighting and material properties
lighting gouraud;
material shiny;
light('Position', [5, 5, 5]);
light('Position', [-5, -5, 5]);

% Set background color
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');

% Set axis properties
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.5);

% Remove axis tick marks and numbers
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'ZTickLabel', []);

% Add subtle grid
grid on;
set(gca, 'GridAlpha', 0.3);
set(gca, 'GridLineStyle', '--');

end

function field_line = calculate_field_line(start_point, magnet_positions, magnet_orientations, ...
                                          magnetic_moment, max_steps, step_size)
% Calculate a magnetic field line starting from start_point
% 
% Inputs:
%   start_point: [x, y, z] starting position
%   magnet_positions: Nx3 array of magnet positions
%   magnet_orientations: Nx3 array of magnet orientation vectors
%   magnetic_moment: scalar magnetic moment strength
%   max_steps: maximum number of integration steps
%   step_size: integration step size
%
% Output:
%   field_line: Mx3 array of points along the field line

field_line = zeros(max_steps, 3);
current_point = start_point;
field_line(1, :) = current_point;

for step = 2:max_steps
    % Calculate magnetic field at current point
    B_field = calculate_magnetic_field(current_point, magnet_positions, ...
                                      magnet_orientations, magnetic_moment);
    
    % Check if field is too weak or point is too far
    B_magnitude = norm(B_field);
    if B_magnitude < 1e-8 || norm(current_point) > 8  % Relaxed conditions
        field_line = field_line(1:step-1, :);
        break;
    end
    
    % Check if too close to any magnet
    min_distance = inf;
    for i = 1:size(magnet_positions, 1)
        distance = norm(current_point - magnet_positions(i, :));
        min_distance = min(min_distance, distance);
    end
    
    if min_distance < 0.25  % Reduced threshold for closer approach
        field_line = field_line(1:step-1, :);
        break;
    end
    
    % Normalize field direction and take step
    field_direction = B_field / B_magnitude;
    current_point = current_point + step_size * field_direction;
    field_line(step, :) = current_point;
end

% Remove unused rows
field_line = field_line(1:step-1, :);

end

function B_field = calculate_magnetic_field(point, magnet_positions, magnet_orientations, magnetic_moment)
% Calculate the total magnetic field at a given point due to all magnets
%
% Inputs:
%   point: [x, y, z] position where to calculate field
%   magnet_positions: Nx3 array of magnet positions
%   magnet_orientations: Nx3 array of magnet orientation vectors
%   magnetic_moment: scalar magnetic moment strength
%
% Output:
%   B_field: [Bx, By, Bz] magnetic field vector

B_field = [0, 0, 0];
mu_0_over_4pi = 1e-7;  % mu_0/(4*pi) simplified

for i = 1:size(magnet_positions, 1)
    magnet_pos = magnet_positions(i, :);
    magnet_orient = magnet_orientations(i, :);
    
    % Vector from magnet to field point
    r_vec = point - magnet_pos;
    r_magnitude = norm(r_vec);
    
    if r_magnitude < 0.1  % Avoid singularity, increased threshold
        continue;
    end
    
    r_unit = r_vec / r_magnitude;
    
    % Magnetic dipole field formula with enhanced strength
    % B = (mu_0 / 4*pi) * (1/r^3) * [3*(m·r̂)r̂ - m]
    m_dot_r = dot(magnet_orient, r_unit);
    
    B_dipole = (mu_0_over_4pi * magnetic_moment * 100) / r_magnitude^3 * ...  % Enhanced by factor of 100
               (3 * m_dot_r * r_unit - magnet_orient);
    
    B_field = B_field + B_dipole;
end

end
