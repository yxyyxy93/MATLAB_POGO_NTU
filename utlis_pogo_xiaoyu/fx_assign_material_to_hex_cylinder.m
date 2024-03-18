function model = fx_assign_material_to_hex_cylinder(model, center, cylinder_center, depth, side_length, material_name)
    % Calculate the vertices of the hexagon in the XY plane
    % Assuming one of the vertices is aligned with the X-axis
    theta = (0:5) * 2 * pi / 6; % Angles for the vertices
    radius = side_length / sqrt(3); % Distance from center to vertex for a regular hexagon
    vertices_x = cylinder_center(1) + radius * cos(theta);
    vertices_y = cylinder_center(2) + radius * sin(theta);
    % Extend the depth (z-axis) range of the hexagonal cylinder
    z_start = cylinder_center(3);
    z_end = cylinder_center(3) + depth;
    % Iterate through all points in the model to check if they are within the hexagonal cylinder
    for i = 1:size(center, 1)
        x = center(i, 1);
        y = center(i, 2);
        z = center(i, 3);
        % Check if the point is within the hexagon in the XY plane and within the depth range
        inHex = inpolygon(x, y, vertices_x, vertices_y) && z >= z_start && z <= z_end;
        % Assign material to elements within the hexagonal cylinder
        if inHex
            model.matTypeRefs(i, 1) = material_name; % Assuming material_name is correctly formatted for assignment
        end
    end
end