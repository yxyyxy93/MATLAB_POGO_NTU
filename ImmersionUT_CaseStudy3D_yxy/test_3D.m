clear;
close all;
clc
fclose all;

PogoFilename = 'test_3D'; 

% Stimulation signal
frequency = 5e6;% unit: Hz
cycles    = 3;
timedelay = 2e-7;
timestep  = 2e-9;
endtime   = 1e-6;
phase     = -pi/2;
filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

%% Generating cubic and meshes
% Grid definition
X_cubic = 2e-3;
Y_cubic = 2e-3;
Z_cubic = 8e-3;
X_mesh  = 2e-5;
Y_mesh  = 2e-5;
Z_mesh  = 2e-5;

% Node position vectors
X_nodePos = 0: X_mesh: X_cubic;
Y_nodePos = 0: Y_mesh: Y_cubic;
Z_nodePos = 0: Z_mesh: Z_cubic;

% Lengths
lx = length(X_nodePos);
ly = length(Y_nodePos);
lz = length(Z_nodePos);

tic;

% 3D grids
[X, Y, Z] = ndgrid(X_nodePos, Y_nodePos, Z_nodePos);

% Node positions
model.nodePos = [X(:)'; Y(:)'; Z(:)'];

% Element nodes
[X, Y, Z] = ndgrid(1:(lx-1), 1:(ly-1), 1:(lz-1));
idx = @(x, y, z) x + (y-1)*lx + (z-1)*lx*ly;

model.elNodes = [
    reshape(idx(X(:),   Y(:),   Z(:)),   1, []); 
    reshape(idx(X(:)+1, Y(:),   Z(:)),   1, []); 
    reshape(idx(X(:)+1, Y(:)+1, Z(:)),   1, []);
    reshape(idx(X(:),   Y(:)+1, Z(:)),   1, []);
    reshape(idx(X(:),   Y(:),   Z(:)+1), 1, []); 
    reshape(idx(X(:)+1, Y(:),   Z(:)+1), 1, []);
    reshape(idx(X(:)+1, Y(:)+1, Z(:)+1), 1, []);
    reshape(idx(X(:),   Y(:)+1, Z(:)+1), 1, [])
    ];

toc;

model.elTypes{1}.name       = 'C3D8R';
model.elTypes{1}.paramsType = 0;
model.nDims                 = 3; % default 2D 
model.nDofPerNode           = 3;
model.elTypeRefs            = ones(length(model.elNodes(1,:)),1);

model.matTypeRefs = ones(length(model.elNodes(1,:)), 1); % water medium

%% Define the CFRP region
X_box_start = 0;    % Starting X position of the box
Y_box_start = 0;    % Starting Y position of the box
Z_box_start = 4e-3; 

X_box_end   = X_cubic; % Ending X position of the box
Y_box_end   = Y_cubic; % Ending Y position of the box
Z_box_end   = 6e-3;

space = (Z_box_end - Z_box_start)/24;

material_name = 2;

%********** define a center out of function for acceleration
% Preallocate array for faster execution
elementCount = length(model.elNodes(1,:));
center       = zeros(elementCount, 3);
% Calculate all centers at once
for i = 1: size(model.elNodes, 1)
    center = center + model.nodePos(:, model.elNodes(i, :))';
end
center = center / size(model.elNodes, 1);
%******************

model         = fx_assign_material_to_box(model, ...
    X_box_start, Y_box_start, Z_box_start, ...
    X_box_end, Y_box_end, Z_box_end, material_name, center);


%% Settings except for material
model.prec    = 8;      % Precision
model.runName = 'Job';
model.nt      = 10e3;
model.dt      = timestep;

%% Material settings
% For isotropi49c material
%-----------------------------------------------------%
cWater = 1500; rhoWater = 1000; visc = 0;

model.matTypes{1,1}.paramsType  = 5; % water media
model.matTypes{1,1}.paramValues = [cWater, rhoWater, visc];

% For isotropi49c material
%-----------------------------------------------------%
model.matTypes{2,1}.paramsType  = 0; % CFRP
model.matTypes{2,1}.paramValues = [13.47e9, 0, 1588]; % virtual parameter!

% % virtual water
% model.matTypes{1,1}.paramsType  = 0; 
% E = cWater^2 * rhoWater;
% model.matTypes{1,1}.paramValues = [E, 0, rhoWater];

%% Generator
model.shots{1, 1}.ntSig = length(tb_signal);
model.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);

Element_size = X_mesh;

% Plane wave
fd       = 25.6e-3 / 4; % m
diameter = 6.35e-3 / 4;

center = [                 ...
    X_cubic/2 ...
    Y_cubic/2 ...
    0];

% Plane wave
dis_to_center = (model.nodePos(1, :) - center(1)).^2 + ...
    (model.nodePos(2, :) - center(2)).^2 + ...
    (model.nodePos(3, :) - center(3)).^2;

node_index_generator = find(   ...
    model.nodePos(3, :) == 0 & ...
    dis_to_center <= (diameter/2).^2);

[focused_waves, delays] = fx_focused_wave(fd, center, timestep, ...
    model.nodePos(1:3, node_index_generator), cWater, tb_signal(:, 2));

for i = 1: size(focused_waves, 1)
    model.shots{1, 1}.sigs{i, 1}.sigType    = 0; % 0 - force, 1 - displacement
    model.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
    model.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
    model.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
    model.shots{1, 1}.sigs{i, 1}.sigAmps    = ones(length(model.shots{1}.sigs{i}.dofSpec),1)*1e-13;
    model.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
    %     model.shots{1, 1}.sigs{i, 1}.sig        = tb_signal(:, 2)';
end

save("focusing_delays.txt", 'delays', '-ascii');

%% Boundary
X_lim_up  = 0;
X_lim_low = 0;
Y_lim_up  = Y_cubic;
Y_lim_low = 0;
Z_lim_up  = Z_cubic;
Z_lim_low = 0;
node_index_generator_yz0 = find(     ...
    (model.nodePos(1,:)>=X_lim_low)& ...
    (model.nodePos(1,:)<=X_lim_up)&  ...
    (model.nodePos(2,:)>=Y_lim_low)& ...
    (model.nodePos(2,:)<=Y_lim_up)&  ...
    (model.nodePos(3,:)>=Z_lim_low)& ...
    (model.nodePos(3,:)<=Z_lim_up));

X_lim_up  = X_cubic;
X_lim_low = X_cubic;
node_index_generator_yzx = find(     ...
    (model.nodePos(1,:)>=X_lim_low)& ...
    (model.nodePos(1,:)<=X_lim_up)&  ...
    (model.nodePos(2,:)>=Y_lim_low)& ...
    (model.nodePos(2,:)<=Y_lim_up)&  ...
    (model.nodePos(3,:)>=Z_lim_low)& ...
    (model.nodePos(3,:)<=Z_lim_up));

X_lim_up  = X_cubic;
X_lim_low = 0;
Y_lim_up  = 0;
Y_lim_low = 0;
Z_lim_up  = Z_cubic;
Z_lim_low = 0;
node_index_generator_xz0 = find(     ...
    (model.nodePos(1,:)>=X_lim_low)& ...
    (model.nodePos(1,:)<=X_lim_up)&  ...
    (model.nodePos(2,:)>=Y_lim_low)& ...
    (model.nodePos(2,:)<=Y_lim_up)&  ...
    (model.nodePos(3,:)>=Z_lim_low)& ...
    (model.nodePos(3,:)<=Z_lim_up));

Y_lim_up  = Y_cubic;
Y_lim_low = Y_cubic;
node_index_generator_xzy = find(     ...
    (model.nodePos(1,:)>=X_lim_low)& ...
    (model.nodePos(1,:)<=X_lim_up)&  ...
    (model.nodePos(2,:)>=Y_lim_low)& ...
    (model.nodePos(2,:)<=Y_lim_up)&  ...
    (model.nodePos(3,:)>=Z_lim_low)& ...
    (model.nodePos(3,:)<=Z_lim_up));

model.fixNodes = [            ...
    node_index_generator_yz0  ...
    node_index_generator_yzx  ...
    node_index_generator_xz0  ...
    node_index_generator_xzy];

hold on;
scatter3(...
    model.nodePos(1, model.fixNodes), ...
    model.nodePos(2, model.fixNodes), ...
    model.nodePos(3, model.fixNodes));

model.fixDof = [                                ...
    ones(length(node_index_generator_yz0),1)*1; ...
    ones(length(node_index_generator_yzx),1)*1; ...
    ones(length(node_index_generator_xz0),1)*2; ...
    ones(length(node_index_generator_xzy),1)*2];

%% absorbing boundary
% Absorbing Regions
% Only isotropic materials supported for SRM (stiffness reduction method)

xLims    = [];
yLims    = [];
zLims    = [Z_cubic Z_cubic-1e-3 100 102];
nAbsVals = 60;
c0       = 1500 ;
freq     = frequency;
model    = addAbsBound(model, xLims , yLims , zLims , [], [], c0, freq);

%% Receiver
node_index_receiver = node_index_generator;


model.measSets{1, 1}.name       = 'main';
model.measSets{1, 1}.isDofGroup = 0;
model.measSets{1, 1}.measDof    = repmat((1:model.nDims)',length(node_index_receiver),1);
model.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model.nDims,1),length(node_index_receiver)*model.nDims,1);
model.measFreq       = 1;
model.measStart      = 1; 
model.fieldStoreIncs  = round((1:1:60) / 60 * model.nt)';
node_num              = length(model.nodePos);
model.fieldStoreNodes = round((1:1e6)/1e6 * node_num)'; 

%% Save pogo-inp file
% model = rmfield(model,'grain_Orientations');
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO

disp(".pogo-inp saved");

% pic = reshape(model.matTypeRefs(:,1),X_cubic/X_mesh,Y_cubic/Y_mesh)-1;
% figure;
% imshow(pic');

close all;
