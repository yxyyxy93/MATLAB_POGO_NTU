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
nz = 1e2;

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

% set the point generator to solid elastic material
node_surface = find(model.nodePos(3, :) >= min(model.nodePos(3,:)) & ...
    model.nodePos(3, :) < min(model.nodePos(3,:)) + dz * 50);

model.matTypeRefs(node_surface) = 2;

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
model.fieldStoreNodes = round((1:1e6) / 1e6 * node_num)';


%% Save pogo-inp file
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO

disp(".pogo-inp saved");

close all;
