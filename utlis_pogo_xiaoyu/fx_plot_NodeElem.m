function fx_plot_NodeElem(new_elems, new_nodes, faceColor, edgeColor)
    % FX_PLOT_NODEELEM Plots a 3D finite element mesh using 'trisurf'.
    %
    % Syntax:
    % fx_plot_NodeElem(new_elems, new_nodes)
    % fx_plot_NodeElem(new_elems, new_nodes, faceColor, edgeColor)
    %
    % Inputs:
    % new_elems - Element connectivity matrix
    % new_nodes - Node coordinates matrix (Nx3)
    % faceColor - (Optional) Face color of the elements
    % edgeColor - (Optional) Edge color of the elements
    
    % Check if face and edge colors are provided, and set defaults if not
    if nargin < 3 || isempty(faceColor)
        faceColor = 'cyan'; % Default face color
    end
    if nargin < 4 || isempty(edgeColor)
        edgeColor = 'black'; % Default edge color
    end
    
    % Plotting the 3D finite element mesh
    figure;
    trisurf(new_elems, new_nodes(:,1), new_nodes(:,2), new_nodes(:,3), ...
        'FaceColor', faceColor, 'EdgeColor', edgeColor);
    axis equal; % Maintain aspect ratio
    title('3D Finite Element Mesh');
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    
    % Optionally, enhance viewing experience
    lighting phong; % Improving the lighting
    camlight headlight; % Adding a light source
    material dull; % Adjusting the material properties
end
