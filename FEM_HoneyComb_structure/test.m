clear;
close all;
fclose all;
clc;

% Define parameters for the hexagonal prism wall
zBottomOuter = 0;
zTopOuter = 4e-1;
cx = 0; % Center X coordinate
cy = 0; % Center Y coordinate
r = 4e-1; % Radius of the hexagon
wall_thickness = 2e-1;

outer_radius = r + (wall_thickness / 2);
inner_radius = r - (wall_thickness / 2);
ratio = inner_radius / outer_radius;
% Display the ratio as a fraction
[numerator, denominator] = rat(ratio);
disp(['Ratio of internal to outer hexagon: ', ...
    num2str(numerator), '/', num2str(denominator)]);

nPoints_inner = numerator; % Number of points for interpolation on the inner prism
nPoints_outer = denominator; % Number of points for interpolation on the outer prism

% ************* create centers array ******************
% Calculations for adjacent hexagons
dx = 1.5 * r; % Horizontal distance between centers
dy = sqrt(3) * r / 2; % Vertical distance for the shift
% Calculations for adjacent hexagons
% Initialize the centers array for a 4x4 structure
rows = 3;
cols = 3;
centers = zeros(rows*cols, 2); % Updated for a 4x4 grid
% Populate the centers array
index = 1;
for row = 0:rows - 1% rows
    for col = 0:cols - 1% columns
        cx_temp = cx + col * dx;
        cy_temp = cy - 2*row*dy - mod(col,2)*dy; % Adjust vertical position
        centers(index, :) = [cx_temp, cy_temp];
        index = index + 1;
    end
end
% Plot each hexagon
figure;
hold on;
axis equal;
for i = 1:size(centers, 1)
    fx_plotHexagon(centers(i, 1), centers(i, 2), r);
end
title('Honeycomb Structure');
xlabel('X coordinate (m)');
ylabel('Y coordinate (m)');
grid on;
% ********************************************************
gm1 = fx_createHexagonalPrismWall(zBottomOuter, zTopOuter, cx, cy, r, wall_thickness, nPoints_inner, nPoints_outer);
% Generate a mesh for the geometry
gm1 = generateMesh(gm1, "GeometricOrder", "linear");
% Plot the geometry with the face labels
figure;
pdegplot(gm1, 'FaceLabels', "on", 'FaceAlpha', 0.5);
figure;
pdemesh(gm1);

% [gm, elements, nodes] = fx_createHoneycombStructure(zBottomOuter, ...
%     zTopOuter, centers, r, wall_thickness, nPoints_inner, nPoints_outer);

gm_voided = fx_HoneyCombwithPlate(zTopOuter, centers, r, wall_thickness);

% % a 3D plot using 'trisurf':
% figure;
% trisurf(elements, nodes(:,1), nodes(:,2), nodes(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'black');
% axis equal;  % To maintain the aspect ratio
% title('3D Finite Element Mesh');
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');

% % Scatter plot for bottom and top frames
% figure,
% scatter3(gm.Vertices(:, 1), gm.Vertices(:, 2), gm.Vertices(:, 3), 'LineWidth', 2);
% % Annotating each point with its index number
% numVertices = size(gm.Vertices, 1); % Total number of vertices
% for idx = 1:numVertices
%     % Adjust the position of the text for better visibility
%     textOffset = 0.001; % Offset for text annotation to avoid overlapping with the point
%     hold on;
%     text(gm.Vertices(idx, 1) + textOffset, gm.Vertices(idx, 2)  + textOffset, ...
%         gm.Vertices(idx, 3) , num2str(idx), 'FontSize', 8, 'HorizontalAlignment', 'center');
% end
 
% %%
% % Open a new file
% fid = fopen('mesh.geo', 'w');
% % Write the points
% for i = 1:size(nodes, 1)
%     fprintf(fid, 'Point(%d) = {%f, %f, %f};\n', i, nodes(i, 1), nodes(i, 2), nodes(i, 3));
% end
% % Initialize a counter for unique line tags
% lineTag = 1;
% % Write the lines or triangles depending on the element connectivity
% for i = 1:size(elements, 1)
%     % Assume elements(i, :) are the vertex indices of a triangular element
%     fprintf(fid, 'Line(%d) = {%d, %d};\n', lineTag, elements(i, 1), elements(i, 2));
%     fprintf(fid, 'Line(%d) = {%d, %d};\n', lineTag+1, elements(i, 2), elements(i, 3));
%     fprintf(fid, 'Line(%d) = {%d, %d};\n', lineTag+2, elements(i, 3), elements(i, 1));
%     fprintf(fid, 'Line Loop(%d) = {%d, %d, %d};\n', i, lineTag, lineTag+1, lineTag+2);
%     fprintf(fid, 'Plane Surface(%d) = {%d};\n', i, i);
%     lineTag = lineTag + 3; % Increment lineTag by 3 for the next set of lines
% end
% % Close the file
% fclose(fid);
% 
% % Run Gmsh to mesh the geometry and output a .msh file
% system('gmsh -3 mesh.geo -o mesh.msh -format msh2');
% 
% %% 
% % Open the .msh file
% fid = fopen('mesh.msh', 'r');
% % Initialize variables to hold node and element data
% gmsh_nodes = [];
% gmsh_elements = [];
% % Read the file line by line
% while ~feof(fid)
%     line = strtrim(fgetl(fid));
%     % Check for the nodes section
%     if strcmp(line, '$Nodes')
%         numNodes = str2double(fgetl(fid)); % The next line contains the number of nodes
%         gmsh_nodes = zeros(numNodes, 4); % Preallocate the nodes array [nodeID, x, y, z]
%         for i = 1:numNodes
%             nodeData = str2num(fgetl(fid)); % Read node data
%             gmsh_nodes(i, :) = nodeData;
%         end
%     end
%     % Check for the elements section
%     if strcmp(line, '$Elements')
%         numElements = str2double(fgetl(fid)); % The next line contains the number of elements
%         gmsh_elements = zeros(numElements, 4); % Preallocate for simplicity, adjust according to your element type
%         idx = 1;
%         for i = 1:numElements
%             elemData = str2num(fgetl(fid)); % Read element data
%             if length(elemData)>=8 
%                 % Assuming all elements are of the same type, adjust indexing as necessary
%                 gmsh_elements(idx, :) = elemData(end-3:end);
%                 idx = idx + 1;
%             end
%         end
%     end
% end
% % Close the file
% fclose(fid);
% % a 3D plot using 'trisurf':
% figure;
% trisurf(gmsh_elements, gmsh_nodes(:,2), gmsh_nodes(:,3), gmsh_nodes(:,4), 'FaceColor', 'cyan', 'EdgeColor', 'black');
% axis equal;  % To maintain the aspect ratio
% title('3D Finite Element Mesh');
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');

%% add plates 
% % Assuming dimensions for plates
% plateThickness = 1e-1; % Example thickness of each plate
% % Calculate the bounds
% xBounds = [min(nodes(:,1)), max(nodes(:,1))];
% yBounds = [min(nodes(:,2)), max(nodes(:,2))];
% zBounds = [min(nodes(:,3)), max(nodes(:,3))];
% 
% numX = 10;
% numY = 10;
% % Mesh grid for plate surfaces
% [xp, yp] = meshgrid(linspace(xBounds(1), xBounds(2), numX), linspace(yBounds(1), yBounds(2), numY));
% xp = xp(:);
% yp = yp(:);
% 
% % Top plate bottom surface Z-coordinate (subtracting thickness from the top surface)
% zTopPlateBottom = zBounds(2);
% zTopPlateTop = zTopPlateBottom + plateThickness;
% % Bottom plate top surface Z-coordinate (adding thickness to the bottom surface)
% zBottomPlateTop = zBounds(1);
% zBottomPlateBottom = zBottomPlateTop - plateThickness;
% % Extend surface points to include thickness for the top plate
% zTopPlatePointsTop = repmat(zTopPlateTop, size(xp));
% zTopPlatePointsBottom = repmat(zTopPlateBottom, size(xp));
% % Extend surface points to include thickness for the bottom plate
% zBottomPlatePointsTop = repmat(zBottomPlateTop, size(xp));
% zBottomPlatePointsBottom = repmat(zBottomPlateBottom, size(xp));
% 
% % Combine into full arrays for the top plate
% topPlateVolumePoints = [xp, yp, zTopPlatePointsTop; xp, yp, zTopPlatePointsBottom];
% % Combine into full arrays for the bottom plate
% bottomPlateVolumePoints = [xp, yp, zBottomPlatePointsBottom; xp, yp, zBottomPlatePointsTop];
% 
% % Optional: Creating alphaShape for visualization or further operations
% alphaRadius = max([xBounds(2)-xBounds(1), yBounds(2)-yBounds(1), plateThickness]) * 0.1; % Adjust as needed
% topPlateShapeVolume = alphaShape(topPlateVolumePoints(:,1), topPlateVolumePoints(:,2), topPlateVolumePoints(:,3), alphaRadius);
% bottomPlateShapeVolume = alphaShape(bottomPlateVolumePoints(:,1), bottomPlateVolumePoints(:,2), bottomPlateVolumePoints(:,3), alphaRadius);
% 
% [elements_top, nodes_top] = boundaryFacets(topPlateShapeVolume);
% [elements_bottom, nodes_bottom] = boundaryFacets(bottomPlateShapeVolume);
% 
% new_nodes = [nodes; nodes_bottom];
% new_elems = [elements; elements_bottom];
% gm = fegeometry(new_nodes, new_elems);
% figure;
% pdegplot(gm, 'FaceLabels', "on", 'FaceAlpha', 0.5);
% 
% % Round the points to the nearest tolerance to group close points together
% tolerance = zTopOuter/1e2; % Adjust this based on your requirements
% roundedPoints = round(new_nodes / tolerance) * tolerance;
% [uniquePoints, ia] = unique(roundedPoints, 'rows', 'stable');
% % Use the index vector ia to select the unique points from the original allPrismPoints
% uniqueAllPrismPoints = allPrismPoints(ia, :);
% 
% fx_plot_NodeElem(new_elems, new_nodes);

% % Optionally, visualize to verify
% figure;
% subplot(1, 2, 1);
% plot(topPlateShapeVolume);
% title('Top Plate Volume');
% 
% subplot(1, 2, 2);
% plot(bottomPlateShapeVolume);
% title('Bottom Plate Volume');
 
%% Generate the mesh with adjusted element sizes
gm_voided = generateMesh(gm_voided, 'GeometricOrder', 'linear', Hmax=2e-1, Hmin=1e-1);

% Plot the geometry with the face labels
figure;
pdegplot(gm_voided, 'FaceLabels', "on", 'FaceAlpha', 0.5);
figure;
pdemesh(gm_voided);

% Determine the z-coordinate for cutting
z_cut = (max(gm_voided.Mesh.Nodes(3,:)) + min(gm_voided.Mesh.Nodes(3,:))) / 2;

% Identify nodes below the cutting plane
nodes_below = gm_voided.Mesh.Nodes(3,:) <= z_cut;

% Identify elements with all vertices below the cutting line
elements_to_keep = all(nodes_below(gm_voided.Mesh.Elements), 1);

% Create new mesh data
newElements = gm_voided.Mesh.Elements(:, elements_to_keep);
newNodes = gm_voided.Mesh.Nodes(:, unique(newElements(:)));

% Adjust element indices to match the new node indexing
[~, loc] = ismember(newElements, unique(newElements));
newElements = reshape(loc, size(newElements));

% Visualize the bottom half mesh
figure;
pdemesh(newNodes, newElements);

%% solve in MATLAB
% model1 = createpde('structural','static-solid');
% geometryFromMesh(model1, (gm.Mesh.Nodes), (gm.Mesh.Elements));
% 
% % Assign material properties (example: structural steel)
% structuralProperties(model1, 'YoungsModulus', 210e9, ...
%     'PoissonsRatio', 0.3, ...
%     'MassDensity', 7850);
% 
% % Apply boundary conditions
% bottomFaceID = 1; % Example ID, replace with the correct one
% % Apply fixed support at the bottom face
% structuralBC(model1,'Face',bottomFaceID,'Constraint','fixed');
% 
% topFaceID = 2;
% % Apply pressure load on the top face
% structuralBoundaryLoad(model1, 'Face', topFaceID, 'Pressure', 1e4);
% 
% % Solve the model
% result = solve(model1);
% 
% % Post-processing: Visualize the displacement magnitude
% figure;
% pdeplot3D(model1, 'ColorMapData', result.Displacement.Magnitude);
% title('Displacement Magnitude');

%%
nx = 0.5e2;
ny = 0.5e2;
nz = 0.5e2;
dx = 1e-5;
dy = 1e-5;
dz = 1e-5;
[ model ] = genTetGrid3D( nx, ny, nz, dx, dy, dz);

% assign the honeycomb nodes
model.nodePos = gm1.Mesh.Nodes;
model.elNodes = gm1.Mesh.Elements;

model.elTypes{1}.name       = 'C3D4';
model.elTypes{1}.paramsType = 0;
model.nDims                 = 3; % default 3D
model.nDofPerNode           = 3;
model.elTypeRefs            = ones(length(model.elNodes(1,:)),1);
model.matTypeRefs           = ones(length(model.elNodes(1,:)), 1); % water medium

%% Stimulation signal
frequency    = 0.5e6;% unit: Hz
cycles       = 3;
timedelay    = 0.5e-5;
timestep     = 1e-8;
endtime      = 2e-5;
phase        = -pi/2;
windowlength = cycles/(frequency);

clear tb_signal;
tb_signal(:,1) = (0:timestep:endtime)';
tb_signal(:,2) = 0;

% Generate toneburst signal
signal = ...
    0.5 * (1 - cos(2*pi*frequency/cycles*tb_signal(1:round(windowlength/timestep+1),1)))...
    .*sin((2*pi*frequency*tb_signal(1:round(windowlength/timestep+1),1))+phase);

% Adjust delay
tb_signal(round(timedelay/timestep+1):round((timedelay+windowlength)/timestep+1), 2) = signal;

% for dispalying the signal
plot(tb_signal(:,1), tb_signal(:,2));

%% Settings except for material
model.prec    = 8;      % Precision
model.runName = 'Job';
model.nt      = 1e5;
model.dt      = timestep;
% Material settings
% For isotropi49c material
%-----------------------------------------------------%
model.matTypes{1,1}.paramsType  = 0;
model.matTypes{1,1}.paramValues = [69e9, 0.33, 2700]; % Aluminium

%% Generator
model.shots{1, 1}.ntSig = length(tb_signal);
model.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);
Element_size = dx;

% point generation
node_index_generator = find( ...
    model.nodePos(3, :) == min(model.nodePos(3,:)) ...
    );

model.shots{1, 1}.sigs{1, 1}.sigType    = 0; % 0 - force, 1 - displacement
model.shots{1, 1}.sigs{1, 1}.isDofGroup = 0;
model.shots{1, 1}.sigs{1, 1}.dofSpec    = ones(length(node_index_generator),1) * 3;
model.shots{1, 1}.sigs{1, 1}.nodeSpec   = node_index_generator';
model.shots{1, 1}.sigs{1, 1}.sigAmps    = ones(length(model.shots{1}.sigs{1}.dofSpec),1)*1e-13;
model.shots{1, 1}.sigs{1, 1}.sig        = tb_signal(:,2)';

%% Receiver
node_index_receiver = node_index_generator;

model.measSets{1, 1}.name       = 'main';
model.measSets{1, 1}.isDofGroup = 0;
model.measSets{1, 1}.measDof    = 3*ones(length(node_index_receiver),1);
model.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver, 1, 1),length(node_index_receiver)*1,1);
model.measFreq  = 1;
model.measStart = 1;
model.fieldStoreIncs  = round((1:1:30) / 30 * model.nt)';
node_num              = length(model.nodePos);
model.fieldStoreNodes = round((1:1e3) / 1e3 * node_num)';

%%
% fx_display_model(model, 1);

%% Save pogo-inp file
PogoFilename = 'test_honeycomb';
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO

disp(".pogo-inp saved");
