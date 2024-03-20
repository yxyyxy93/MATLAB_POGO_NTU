clear;
close all;
clc
fclose all;

PogoFilename = 'test_3D';

% Stimulation signal
frequency    = 5e6;% unit: Hz
cycles       = 3;
timedelay    = 2e-7;
timestep     = 1e-9;
endtime      = 1e-6;
phase        = -pi/2;
windowlength = cycles/(frequency);

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

%% Generating cubic and meshes
% Grid definition
nx = 1e2;
ny = 1e2;
nz = 0.2e2;
dx = 2e-5; 
dy = 2e-5;
dz = 2e-5;
model = genGrid3D(nx, ny, nz, dx, dy, dz);
model.elTypes{1}.name       = 'C3D8R';
model.elTypes{1}.paramsType = 0;
model.nDims                 = 3; % default 3D 
model.nDofPerNode           = 3;
model.elTypeRefs            = ones(length(model.elNodes(1,:)),1);
model.matTypeRefs = ones(length(model.elNodes(1,:)), 1); % water medium

%% ********** define a center out of function for acceleration
% Preallocate array for faster execution
tic;
% Get the node positions for all elements at once
nodes_pos_all      = model.nodePos(:, model.elNodes(:));
% Reshape the array to separate each element's nodes
nodes_pos_reshaped = reshape(nodes_pos_all, size(model.nodePos, 1), size(model.elNodes, 1), []);
% Compute the centroids by taking the mean across the second dimension
centroids          = squeeze(mean(nodes_pos_reshaped, 2))';
% centroids          = squeeze(nodes_pos_reshaped(:,1,:))';
X = int32(centroids(:, 1)*1e6);
Y = int32(centroids(:, 2)*1e6);
Z = int32(centroids(:, 3)*1e6);
clear nodes_pos_all;
clear nodes_pos_reshapes;
% X_unique = unique(X);
% X_len = max(X) - min(X);
% Y_len = max(Y) - min(Y);
toc;

%% Settings except for material
model.prec    = 8;      % Precision
model.runName = 'Job';
model.nt      = 4e3;
model.dt      = timestep;

%% Material settings
% For isotropi49c material
%-----------------------------------------------------%
cWater = 1500; rhoWater = 1000; visc = 0;

% 
model.matTypes{1,1}.paramsType  = 5; % water media
model.matTypes{1,1}.paramValues = [cWater, rhoWater, visc];
% ********** uncomment below to change it to solid elastic material type
% virtual water
model.matTypes{2,1}.paramsType  = 0;
E = cWater^2 * rhoWater;
model.matTypes{2,1}.paramValues = [E, 0, rhoWater];

%% Generator
model.shots{1, 1}.ntSig = length(tb_signal);
model.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);

Element_size = dx;

% Plane wave
diameter = 5e-4;
z_loc    = min(model.nodePos(3,:)); % minimum position along z axis

center = [    ...
    mean(model.nodePos(1, :)) ...
    mean(model.nodePos(2, :)) ...
    z_loc];

% Plane wave
dis_to_center = (model.nodePos(1, :) - center(1)).^2 + ...
    (model.nodePos(2, :) - center(2)).^2 + ...
    (model.nodePos(3, :) - center(3)).^2;

node_index_generator = find(   ...
    model.nodePos(3, :) == z_loc & ...
    dis_to_center <= (diameter/2).^2);

% % point generation
% node_index_generator = find( ...
%     model.nodePos(1, :) == min(abs(model.nodePos(1,:))) & ...
%     model.nodePos(2, :) == min(abs(model.nodePos(2,:))) & ...
%     model.nodePos(3, :) == min(model.nodePos(3,:)) ...
%     );

% % set the point generator to solid elastic material
% node_surface = find(model.nodePos(3, :) >= min(model.nodePos(3,:)) & ...
%     model.nodePos(3, :) < min(model.nodePos(3,:)) + dz * 50);
% model.matTypeRefs(node_surface) = 2;

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
model.measSets{1, 1}.measDof    = repmat((1:model.nDims)',length(node_index_receiver),1);
model.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model.nDims,1),length(node_index_receiver)*model.nDims,1);
model.measFreq  = 1;
model.measStart = 1; 
model.fieldStoreIncs  = round((1:1:30) / 30 * model.nt)';
node_num              = length(model.nodePos);
model.fieldStoreNodes = round((1:1e3) / 1e3 * node_num)';

%%
d_plate = 2 * dz;

centerX = 0;
centerY = 0;
centerZ = min(model.nodePos(3,:));  % Starting at the bottom
hexCenter = [centerX, centerY, centerZ];
hexDepth = nz * dz - d_plate;
hexSideLength = 5e-4; % Example side length, adjust based on your requirements
material_name = 2;

% assign materials
xyz_min = min(centroids, [], 1);
xyz_max = max(centroids, [], 1);
x_min = xyz_min(1);
y_min = xyz_min(2);
z_min = xyz_min(3);
x_max = xyz_max(1);
y_max = xyz_max(2);
z_max = xyz_max(3);

figure, plot(model.matTypeRefs);

model = fx_assign_material_to_box(model, x_min, y_min, z_max-d_plate, ...
    x_max, y_max, z_max, material_name, centroids);
figure, plot(model.matTypeRefs);

% Define your model's central point offset
center_offset = [0, -1e-3];
num_rows = 4; % The number of rows you want
hex_centers = fx_calculate_honeycomb_centers(num_rows, hexSideLength, center_offset, centerZ);

inHex_ele = fx_find_elements_in_hexcylinder(int32(hex_centers*1e6), ...
    int32(hexDepth*1e6), hexSideLength*1e6, X, Y, Z);

[ mDel ] = deleteEls(model, inHex_ele);


% for i = 1:size(hex_centers, 1)
%     hexCenter = hex_centers(i, :);
%     model = fx_assign_material_to_hex_cylinder(model, center, [hexCenter centerZ], hexDepth, hexSideLength, material_name);
% end

% test functions
clc;
close all;
flag_plotelements = 1;
fx_display_model(mDel, flag_plotelements);

% %% Save pogo-inp file
% savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO
% 
% disp(".pogo-inp saved");
% 
% close all;
