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
gm_voided = generateMesh(gm_voided, 'GeometricOrder', 'linear', Hmax=1e-3, Hmin=1e-5);
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
model.nodePos = gm_voided.Mesh.Nodes;
model.elNodes = gm_voided.Mesh.Elements;

% assign the honeycomb nodes --- for visualize the inside !!
model.nodePos = newNodes;
model.elNodes = newElements;

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

%% add absorbing boundary
% Only isotropic materials supported for SRM (stiffness reduction method)
xyz_min = min(model.nodePos, [], 2);
xyz_max = max(model.nodePos, [], 2);
x_min = xyz_min(1);
y_min = xyz_min(2);
x_max = xyz_max(1);
y_max = xyz_max(2);

nAbsVals = 60;
abs_size = wall_thickness/2;
xLims    = [x_min x_min+abs_size x_max-abs_size x_max];
ylims    = [y_min y_min+abs_size y_max-abs_size y_max];
zLims    = [-102 -100 100 102];
model    = addAbsBound(model, xLims, ylims, zLims, nAbsVals, [], []);

%%
fx_display_model(model, 1);

%% Save pogo-inp file
PogoFilename = 'test_honeycomb';
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO

disp(".pogo-inp saved");


%% run pogo
% Command for running pogoBlockGreedy3d and pogoSolve3d
inputFile = sprintf('%s.pogo-inp', PogoFilename);
runBlockCommand = sprintf('pogoBlockGreedy3d %s', inputFile);
runSolveCommand = sprintf('pogoSolve3d %s.pogo-inp --setFieldSaveOff --setOutputFile "%s" <<< o', PogoFilename, PogoFilename);
% Execute the block command
[status_block, cmdout_block] = system(runBlockCommand);
% Check if the block command was successful before proceeding
if status_block == 0
    disp('pogoBlockGreedy3d executed successfully.');
else
    disp('Error in pogoBlockGreedy3d:');
    disp(cmdout_block);
end
% Execute the solve command if block command was successful
if status_block == 0
    [status_solve, cmdout_solve] = system(runSolveCommand);
    % Check if the solve command was successful
    if status_solve == 0
        disp('pogoSolve3d executed successfully.');
    else
        disp('Error in pogoSolve3d:');
        disp(cmdout_solve);
    end
end
% Command for removing the block file
removeBlockCommand = sprintf('rm %s.pogo-block', PogoFilename);
% Execute the remove block file command if solve command was successful
if status_solve == 0
    [status_rm, cmdout_rm] = system(removeBlockCommand);
    
    % Check if the remove block file command was successful
    if status_rm == 0
        disp('Block file removed successfully.');
    else
        disp('Error removing block file:');
        disp(cmdout_rm);
    end
end

%% read the file 
% Specify the file name
file_name = 'test_honeycomb.pogo-hist';

% Check if the file exists in the current folder
if ~isfile(file_name)
    disp('The specified file does not exist in the current directory.');
else
    % Display the selected file
    fprintf('Selected file: %s\n', file_name);
    % Load the .pogo-hist file
    h = loadPogoHist(file_name);
    % Time vector
    t = h.startMeas + (0:h.nt-1) * h.dt;
    % Take all the wave signals
    ascans = h.sets.main.histTraces;
    % Average of the A-scans
    ascan = mean(ascans, 2);
    % Envelope and phase of the A-scan
    inam  = abs(hilbert(ascan));
    inph  = angle(hilbert(ascan));
    % Plotting the averaged A-scan with its envelope
    figure,
    subplot(2, 1, 1);
    plot(t, ascan);
    hold on;
    plot(t, inam);
    title('Averaged A-scan and Envelope');
    xlabel('Time');
    ylabel('Amplitude');
    % Plotting the phase of the A-scan
    subplot(2, 1, 2);
    plot(t, inph, "Color", 'magenta');
    title('Phase of Averaged A-scan');
    xlabel('Time');
    ylabel('Phase');
    % Plotting the mean of all A-scans (without delay adjustment)
    figure,
    plot(t, mean(ascans, 2));
    title('Mean of All A-scans');
    xlabel('Time');
    ylabel('Amplitude');
end


