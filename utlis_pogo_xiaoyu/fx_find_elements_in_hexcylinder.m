function indices = fx_find_elements_in_hexcylinder(cylinder_centers, depth, side_length, X, Y, Z)
% Calculate the vertices of the hexagon in the XY plane
% Assuming one of the vertices is aligned with the X-axis
theta = (0:5) * 2 * pi / 6; % Angles for the vertices
radius = side_length / sqrt(3); % Distance from center to vertex for a regular hexagon
% Initialize inHex as a logical array with false values
inHex = false(length(X), 1);

for i = 1:size(cylinder_centers, 1)
    vertices_x = cylinder_centers(i, 1) + int32(radius * cos(theta));
    vertices_y = cylinder_centers(i, 2) + int32(radius * sin(theta));
    % Extend the depth (z-axis) range of the hexagonal cylinder
    z_start = cylinder_centers(i, 3);
    z_end = cylinder_centers(i, 3) + depth;
    % Vectorized check for elements within the hexagonal cylinder
    inHex = inHex | inpolygon(X, Y, vertices_x, vertices_y) & Z >= z_start & Z <= z_end;
end
% Find the indices of the elements within the hexagonal cylinder
indices = int32(find(inHex));

end