function model = fx_assign_material_to_box(model, X_box_start, Y_box_start, Z_box_start, X_box_end, Y_box_end, Z_box_end, ...
    material_name, center)
% Define a new material
% X_box_start, Y_box_start, Z_box_start, X_box_end, Y_box_end, Z_box_end be vectors ..

% % Preallocate array for faster execution
% elementCount = length(model.elNodes(1,:));
% center = zeros(elementCount, 3);
% 
% % Calculate all centers at once
% for i = 1:4
%     center = center + model.nodePos(:, model.elNodes(i, :))';
% end
% center = center / 4;

% assign multiple boxes
for j = 1:length(X_box_start)
    % Find elements in the box
    inBox = center(:, 1) > X_box_start(j) & center(:, 1) < X_box_end(j) & ...
            center(:, 2) > Y_box_start(j) & center(:, 2) < Y_box_end(j) & ...
            center(:, 3) > Z_box_start(j) & center(:, 3) < Z_box_end(j);
    
    % Assign material to elements in the box
    model.matTypeRefs(inBox, 1) = material_name; % solid media
end

end
