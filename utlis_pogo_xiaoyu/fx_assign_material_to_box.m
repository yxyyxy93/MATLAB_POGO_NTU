function model = fx_assign_material_to_box(model, X_box_start, Y_box_start, Z_box_start, ...
    X_box_end, Y_box_end, Z_box_end, material_name, center)

% Define a new material
% assign multiple boxes
for j = 1:length(X_box_start)
    dims = size(center);

    % check if 2D
    if length(dims)==2
        inBox = center(:, 1) >= X_box_start(j) & center(:, 1) <= X_box_end(j) & ...
            center(:, 2) >= Y_box_start(j) & center(:, 2) <= Y_box_end(j);
        % check 3D or higher
    else
        % Find elements in the box
        inBox = center(:, 1) >= X_box_start(j) & center(:, 1) <= X_box_end(j) & ...
            center(:, 2) >= Y_box_start(j) & center(:, 2) <= Y_box_end(j) & ...
            center(:, 3) >= Z_box_start(j) & center(:, 3) <= Z_box_end(j);
    end

    % Assign material to elements in the box
    if j<length(material_name)
        model.matTypeRefs(inBox, 1) = material_name(j); % solid media
    else
        model.matTypeRefs(inBox, 1) = material_name(end); % solid media
    end

end
