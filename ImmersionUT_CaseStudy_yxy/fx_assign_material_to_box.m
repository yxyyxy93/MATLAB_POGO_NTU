function model = fx_assign_material_to_box(model, X_box_start, Y_box_start, Z_box_start, X_box_end, Y_box_end, Z_box_end, ...
    material_name)
% Define a new material
% X_box_start, Y_box_start, Z_box_start, X_box_end, Y_box_end, Z_box_end be vectors ..

% calculate centers of all elements at once
center_x = mean(model.nodePos(1,model.elNodes),1);
center_y = mean(model.nodePos(2,model.elNodes),1);
center_z = mean(model.nodePos(3,model.elNodes),1);

% convert centers to column vectors
center_x = center_x(:);
center_y = center_y(:);
center_z = center_z(:);

% create a matrix of all centers for comparison
centers = [center_x center_y center_z];

% loop through each box
for j = 1:length(X_box_start)
    % create a matrix for comparison
    box_start = [X_box_start(j), Y_box_start(j), Z_box_start(j)];
    box_end = [X_box_end(j), Y_box_end(j), Z_box_end(j)];

    % find the indices of elements within the box
    in_box = all(bsxfun(@ge, centers, box_start) & bsxfun(@le, centers, box_end), 2);

    % assign the material to those elements
    model.matTypeRefs(in_box) = material_name; % solid media
end

end
