function model = fx_assign_material_to_rectangle(model, X_rect_start, Y_rect_start, X_rect_end, Y_rect_end, ...
    material_name)
% Define a new material
% X_rect_start, Y_rect_start, X_rect_end, Y_rect_end be vectors ..
for i=1:length(model.elNodes(1,:))
    center_x = (model.nodePos(1,model.elNodes(1,i))+model.nodePos(1,model.elNodes(2,i))+...
        model.nodePos(1,model.elNodes(3,i))+model.nodePos(1,model.elNodes(4,i)))/4;

    center_y = (model.nodePos(2,model.elNodes(1,i))+model.nodePos(2,model.elNodes(2,i))+...
        model.nodePos(2,model.elNodes(3,i))+model.nodePos(2,model.elNodes(4,i)))/4;

    % assign multiply rectangulars
    for j = 1:length(X_rect_start)
        if center_x >= X_rect_start(j) && center_x <= X_rect_end(j) ...
                && center_y >= Y_rect_start(j) && center_y <= Y_rect_end(j)
            model.matTypeRefs(i,1) = material_name; % solid media
        end
    end
end

end