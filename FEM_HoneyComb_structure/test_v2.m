clear;
close all;
fclose all;
clc;

% Define parameters for the hexagonal prism wall
zBottomOuter = 0;
zTopOuter = 1e-2;
cx = 0; % Center X coordinate
cy = 0; % Center Y coordinate
r = 1e-2; % Radius of the hexagon
wall_thickness = 1e-3;

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

gm_voided = fx_HoneyCombwithPlate(zTopOuter, rows, cols, r, wall_thickness);

%% Generate the mesh with adjusted element sizes
gm_voided = generateMesh(gm_voided, 'GeometricOrder', 'linear', Hmax=1e-4, Hmin=1e-5);
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
