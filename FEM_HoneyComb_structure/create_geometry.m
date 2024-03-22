clear;
close all;
fclose all;
clc;

%%
% Define parameters
h = 0.5; % Height of the walls
% Parameters for the hexagonal ring
r_inner = 1; % Inner radius
t = 0.2; % Wall thickness
r_outer = r_inner + t; % Outer radius
N = 6; % Number of sides

theta = linspace(0, 2*pi, N+1); % Angles for vertices
theta(end) = []; % Remove the duplicate end point

% Inner edge vertices
inner_x = r_inner * cos(theta);
inner_y = r_inner * sin(theta);

% Outer edge vertices
outer_x = r_outer * cos(theta);
outer_y = r_outer * sin(theta);

z_bottom = zeros(1, 2*N);
z_top    = h * ones(1, 2*N);

% Plotting the bottom and top frames of the hexagon
figure;
% Order: bottom inner, bottom outer, top inner, top outer
vertices = [inner_x, outer_x, inner_x, outer_x; ...
            inner_y, outer_y, inner_y, outer_y; ...
            z_bottom, z_top];

% Scatter plot for bottom and top frames
scatter3(vertices(1,:), vertices(2,:), vertices(3,:), 'LineWidth', 2);
% Annotating each point with its index number
numVertices = size(vertices, 2); % Total number of vertices
for idx = 1:numVertices
    % Adjust the position of the text for better visibility
    textOffset = 0.1; % Offset for text annotation to avoid overlapping with the point
    text(vertices(1,idx) + textOffset, vertices(2,idx) + textOffset, vertices(3,idx), num2str(idx), ...
        'FontSize', 8, 'HorizontalAlignment', 'center');
end

% Formatting the plot
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Honeycomb Wall Structure with Vertex Indices');
axis equal;
grid on;
view(3); % 3D view
hold off;

%% 
% Combine vertices for 3D representation
% Side faces connect bottom and top vertices
% Each side face is defined by four vertices: [bottomInner, bottomOuter, topOuter, topInner]
% side faces
faces = nan(4, 24);
% The index here is the current point, and the value at that index is the next point
next_point_mapping = ...
    [2, 3, 4, 5, 6, 1, 8, 9, 10, 11, 12, 7 ...
    14, 15, 16, 17, 18, 13];
for i = 1:6 % inner and bottom
    next_point = next_point_mapping(i);
    faces(:, i) = [i next_point next_point+12 i+12];
    faces(:, i+12) = [i next_point next_point+6 i+6];
end
for i = 7:12 % outer
    next_point = next_point_mapping(i);
    faces(:, i) = [i next_point next_point+12 i+12];
end
for i = 13:18 % top
    next_point = next_point_mapping(i);
    faces(:, i+6) = [i next_point next_point+6 i+6];
end

% Plot
figure;
patch('Vertices', vertices.', 'Faces', faces.', 'FaceColor', 'b');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%%
% vertices should be a 3-by-N matrix (N is the number of vertices)
% faces should be an M-by-4 matrix for quadrilaterals (M is the number of faces),
% or M-by-3 for triangular faces, specifying vertex indices

% Create a PDE model
model = createpde();

% Use geometryFromMesh to create the geometry from vertices and faces
gm = geometryFromMesh(model, vertices, faces);

% Assign the geometry to the model
model.Geometry = gm;

% Visualizing the geometry
figure;
pdeplot3D(model);
title('3D Visualization of Hexagonal Prism');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;