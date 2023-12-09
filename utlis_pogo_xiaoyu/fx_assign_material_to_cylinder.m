function model = fx_assign_material_to_cylinder(model, cylinder_start_z, cylinder_end_z, radius, material_name, center, cylinder_center)
    % Assign multiple cylinders
    for j = 1:length(cylinder_start_z)
        % Calculate squared distance from cylinder axis (only XY plane)
        squaredDistance = (center(:, 1) - cylinder_center(j, 1)).^2 + (center(:, 2) - cylinder_center(j, 2)).^2;

        % Check if points are within the cylinder
        inCylinder = squaredDistance <= radius(j)^2 & ...
                     center(:, 3) >= cylinder_start_z(j) & ...
                     center(:, 3) <= cylinder_end_z(j);

        % Assign material to elements in the cylinder
        if j < length(material_name)
            model.matTypeRefs(inCylinder, 1) = material_name(j); % solid media
        else
            model.matTypeRefs(inCylinder, 1) = material_name(end); % solid media
        end
    end
end


