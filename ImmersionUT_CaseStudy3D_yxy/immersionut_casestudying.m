
clear;
close all;
clc
fclose all;

PogoFilename = 'immersionUT_CaseStudy_yxy'; 

% Stimulation signal
frequency = 5e6;% unit: Hz
cycles    = 3;
timedelay = 2e-7;
timestep  = 2e-10;
endtime   = 1e-6;
phase     = -pi/2;
filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

%% Generating cubic and meshes
% Grid definition
X_cubic = 1e-3;
Y_cubic = 1e-3;
Z_cubic = 1e-3;
X_mesh  = 1e-5;
Y_mesh  = 1e-5;
Z_mesh  = 1e-5;

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

model.elNodes = [reshape(idx(X(:),   Y(:),   Z(:)),   1, []); 
                 reshape(idx(X(:)+1, Y(:),   Z(:)),   1, []); 
                 reshape(idx(X(:)+1, Y(:)+1, Z(:)),   1, []);
                 reshape(idx(X(:),   Y(:)+1, Z(:)),   1, []);
                 reshape(idx(X(:),   Y(:),   Z(:)+1), 1, []); 
                 reshape(idx(X(:)+1, Y(:),   Z(:)+1), 1, []);
                 reshape(idx(X(:)+1, Y(:)+1, Z(:)+1), 1, []);
                 reshape(idx(X(:),   Y(:)+1, Z(:)+1), 1, [])];

toc;

model.elTypes{1}.name       = 'C3D8R';
model.elTypes{1}.paramsType = 0;
model.nDims                 = 3; % default 2D 
model.nDofPerNode           = 3;
model.elTypeRefs            = ones(length(model.elNodes(1,:)),1);

% model.grain_Orientations(:,1) = [0 0 0]'; %<100>
% model.grain_Orientations(:,2) = [0 45 0]'; %<110>
% model.grain_Orientations(:,1) = [0 0 0]'; %<100>

% model.grain_Orientations(:,2) = [45 0 0]'; %<110>
% lamellar_width1 = 5e-3;%<110>
% lamellar_width2 = 5e-3;%<100>

model.matTypeRefs = ones(length(model.elNodes(1,:)), 1); % water medium

%% Select elements that in the range of solid parts and change their matTypeRefs
% to 2, which refers to the elastic media.

% Define the CFRP region
X_box_start   = 0;       % Starting X position of the box
Y_box_start   = 0;       % Starting Y position of the box
Z_box_start   = 5e-4;
X_box_end     = X_cubic; % Ending X position of the box
Y_box_end     = Y_cubic; % Ending Y position of the box
Z_box_end     = 8e-4;
material_name = 2;
model         = fx_assign_material_to_box(model, ...
    X_box_start, Y_box_start, Z_box_start, X_box_end, Y_box_end, Z_box_end, material_name);

% Define interplies
material_name = 3;
Z_box_start2  = Z_box_start-10e-6: 230e-6: Z_box_end; 
X_box_start2  = 0 * ones(1, length(Z_box_start2));       % Starting X position of the box
Y_box_start2  = 0 * ones(1, length(Z_box_start2));       % Starting Y position of the box

Z_box_end2 = Z_box_start: 230e-6: Z_box_end + 10e-6; % Ending Y position of the rectangle
X_box_end2 = X_cubic * ones(1, length(Z_box_end2)); % Ending X position of the box
Y_box_end2 = Y_cubic * ones(1, length(Z_box_end2)); % Ending X position of the box

model      = fx_assign_material_to_box(model, ...
    X_box_start2, Y_box_start2, Z_box_start2, X_box_end2, Y_box_end2, Z_box_end2, material_name);

%% Settings except for material
model.prec    = 8;      % Precision
model.runName = 'Job';
model.nt      = 1e4;
model.dt      = 0.5e-9;

%% Material settings
% For isotropi49c material
%-----------------------------------------------------%
cWater = 1500; rhoWater = 1000; visc = 0;

model.matTypes{1,1}.paramsType  = 5; % water media
model.matTypes{1,1}.paramValues = [cWater, rhoWater, visc];

% % virtual water
% model.matTypes{1,1}.paramsType  = 0; 
% E = cWater^2 * rhoWater;
% model.matTypes{1,1}.paramValues = [E, 0, rhoWater];

model.matTypes{2,1}.paramsType  = 0; % CFRP
% model.matTypes{2,1}.paramValues = [200e9,0.3,7800];
model.matTypes{2,1}.paramValues = [13.47e9, 0, 1588]; % virtual parameter!

model.matTypes{3,1}.paramsType  = 0; % resin-rich interply
model.matTypes{3,1}.paramValues = [3.7e9, 0.4, 1270];

%-----------------------------------------------------%
% For anisotropic material
%-----------------------------------------------------%
% model.grain_Orientations=zeros(3,100);
% C11=243.1e9;C12=138.1e9;C44=121.9e9; %Fe
% % C11=234.6e9;C12=145.4e9;C44=126.2e9; %Inconel

% ES=[C11 C12 C11 C12 C12 C11 0 0 0 C44 0 0 0 0 C44 0 0 0 0 0 C44]; % 21 elastic constant for cubic FE

% SM_origin=[ES(1) ES(2) ES(4) ES(7) ES(11) ES(16);
%            ES(2) ES(3) ES(5) ES(8) ES(12) ES(17);
%            ES(4) ES(5) ES(6) ES(9) ES(13) ES(18);
%            ES(7) ES(8) ES(9) ES(10) ES(14) ES(19);
%            ES(11) ES(12) ES(13) ES(14) ES(15) ES(20);
%            ES(16) ES(17) ES(18) ES(19) ES(20) ES(21);];

% for i=1:length(model.grain_Orientations(1,:))
%     SM_rotated(:,:)=StiffnessMatrixRotate2(SM_origin,-model.grain_Orientations(1,i)*pi/180,-model.grain_Orientations(2,i)*pi/180,-model.grain_Orientations(3,i)*pi/180);
%     ES_rotated(i,:)=SM_rotated(:)';
% end
% ES_rotated(:,30)=[];ES_rotated(:,23:24)=[];ES_rotated(:,16:18)=[];ES_rotated(:,9:12)=[];ES_rotated(:,2:6)=[];

% matTypeRefs=model.matTypeRefs;
% model = rmfield(model,'matTypeRefs');
% model.matTypeRefs(1:length(model.elNodes),:)=matTypeRefs;

% Poissonratio=zeros(length(model.grain_Orientations),1);
% Poissonratio(:)=0.3;
% Density=zeros(length(model.grain_Orientations),1);
% Density(:)=7800;
% Density(:)=8260;
% for i=1:length(model.grain_Orientations(1,:))
% model.matTypes{i,1}.paramsType=2;%2 for anistropic
% model.matTypes{i,1}.paramValues=[ES_rotated(i,:),Density(i)];
% end

% Absorbing Regions
% Only isotropic materials supported for SRM (stiffness reduction method)
% xLims    = [];
% yLims    = [];
% zLims    = [Z_cubic Z_cubic-2e-3 100 102];
% nAbsVals = 60;
% c0       = 1500 ;
% freq     = frequency;
% model    = addAbsBound( model, xLims , yLims , zLims , [], [], c0, freq);

%% Generator
model.shots{1, 1}.ntSig = length(tb_signal);
model.shots{1, 1}.dtSig = tb_signal(2,1)-tb_signal(1,1);
node_index_generator    = [];

% % Single transducer
% for i        = 16:18
%     Element_size = 10e-6;
%     X_pos        = 0;
%     Y_pos        = 0.002075+(i-1)*0.31e-3;
%     % X_lim_up   = max(model.nodePos(1,:));
%     % X_lim_low  = min(model.nodePos(1,:));
%     X_lim_up     = X_pos+Element_size/2;
%     X_lim_low    = X_pos-Element_size/2;
%     Y_lim_up     = Y_pos+Element_size/2;
%     Y_lim_low    = Y_pos-Element_size/2;
%     node_index_generator2=find(...
%         (model.nodePos(1,:)>=X_lim_low)&...
%         (model.nodePos(1,:)<=X_lim_up)&...
%         (model.nodePos(2,:)>=Y_lim_low)&...
%         (model.nodePos(2,:)<=Y_lim_up));
%     node_index_generator = [node_index_generator node_index_generator2];
% end

Element_size = X_mesh;

% Plane wave
X_lim_up  = 8e-4+Element_size/2;
X_lim_low = 2e-4-Element_size/2;
Y_lim_up  = 8e-4+Element_size/2;
Y_lim_low = 2e-4-Element_size/2;
Z_lim_up  = 0+Element_size/2;
Z_lim_low = 0-Element_size/2;

node_index_generator = find( ...
    (model.nodePos(1,:) >= X_lim_low) & ...
    (model.nodePos(1,:) <= X_lim_up)  & ...
    (model.nodePos(2,:) >= Y_lim_low) & ...
    (model.nodePos(2,:) <= Y_lim_up)  & ...
    (model.nodePos(3,:) >= Z_lim_low) & ...
    (model.nodePos(3,:) <= Z_lim_up));

% Calculate local spacing
% [x_coor_sorted,index_sorted]=sort(model.nodePos(1,node_index_generator));
% LocalSpacing=(x_coor_sorted(3:end)-x_coor_sorted(1:end-2))/2;
% %             [y_coor_sorted,index_sorted]=sort(model.nodePos(2,node_index_generator));
% %             LocalSpacing=(y_coor_sorted(3:end)-y_coor_sorted(1:end-2))/2;
% LocalSpacing=5e-5;
% node_index_generator=node_index_generator(index_sorted(2:end-1));

hold on;
scatter3(model.nodePos(1,node_index_generator), ...
    model.nodePos(2,node_index_generator), ...
    model.nodePos(3,node_index_generator));

fd     = 25.4e-3; % m
% diameter = 6.35e-3; 
center = [(X_lim_up+X_lim_low)/2 (Y_lim_low+Y_lim_up)/2 (Z_lim_low+Z_lim_up)/2];
focused_waves = fx_focused_wave(fd, center, timestep, ...
    model.nodePos(1:3, node_index_generator), cWater, tb_signal(:,2));

figure;
plot(focused_waves(1, :)');
hold on;
plot(focused_waves(2001, :)');

for i = 1: size(focused_waves, 1)
    model.shots{1, 1}.sigs{i, 1}.sigType    = 0; % 0 - force, 1 - displacement
    model.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
    model.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
    model.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
    model.shots{1, 1}.sigs{i, 1}.sigAmps    = ones(length(model.shots{1}.sigs{i}.dofSpec), 1)*1e-13;
    model.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
end

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
scatter3(model.nodePos(1, model.fixNodes), ...
        model.nodePos(2, model.fixNodes), ...
        model.nodePos(3, model.fixNodes));

model.fixDof = ones(length(model.fixNodes),1)*1;

%% Receiver
node_index_receiver = node_index_generator;

% node_index_receiver=[];
% X_lim = linspace(16.4e-3,54.2e-3,64);
% for i = 1:length(X_lim)
% Element_size = 50e-6;
% X_pos = X_lim(i);
% Y_pos = 20e-3;
% X_lim_up = X_pos+Element_size/2;
% X_lim_low = X_pos-Element_size/2;
% Y_lim_up = Y_pos+Element_size/2;
% Y_lim_low = Y_pos-Element_size/2;
% node_index_receiver2=find(...
%                 (model.nodePos(1,:)>=X_lim_low)&...
%                 (model.nodePos(1,:)<=X_lim_up)&...
%                 (model.nodePos(2,:)>=Y_lim_low)&...
%                 (model.nodePos(2,:)<=Y_lim_up));
% node_index_receiver=[node_index_receiver node_index_receiver2];
% end
% X_lim_up = X_cubic/2+10e-6*3;
% X_lim_low = X_cubic/2-10e-6*3;
% Y_lim_up = 0;
% Y_lim_low = 0;
% Z_lim_up = Z_cubic/2 + 25e-6*3;
% Z_lim_low = Z_cubic/2 - 25e-6*3;
% node_index_receiver = find(...
%             (model.nodePos(1,:)>=X_lim_low)&...
%             (model.nodePos(1,:)<=X_lim_up)&...
%             (model.nodePos(2,:)>=Y_lim_low)&...
%             (model.nodePos(2,:)<=Y_lim_up)&...
%             (model.nodePos(3,:)>=Z_lim_low)&...
%             (model.nodePos(3,:)<=Z_lim_up)...
%             );
%
% X_lim_up = X_cubic/2+10e-6*3;
% X_lim_low = X_cubic/2-10e-6*3;
% Y_lim_up = Y_cubic;
% Y_lim_low = Y_cubic;
% Z_lim_up = Z_cubic/2+25e-6*3;
% Z_lim_low = Z_cubic/2-25e-6*3;
% node_index_receiver2 = find(...
%             (model.nodePos(1,:)>=X_lim_low)&...
%             (model.nodePos(1,:)<=X_lim_up)&...
%             (model.nodePos(2,:)>=Y_lim_low)&...
%             (model.nodePos(2,:)<=Y_lim_up)&...
%             (model.nodePos(3,:)>=Z_lim_low)&...
%             (model.nodePos(3,:)<=Z_lim_up)...
%             );
%
% node_index_receiver=[node_index_receiver node_index_receiver2];
hold on;
scatter(model.nodePos(1,node_index_receiver), model.nodePos(2,node_index_receiver));

model.measSets{1, 1}.name       = 'main';
model.measSets{1, 1}.isDofGroup = 0;
model.measSets{1, 1}.measDof    = repmat((1:model.nDims)',length(node_index_receiver),1);
model.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model.nDims,1),length(node_index_receiver)*model.nDims,1);
model.measFreq       = 1;
model.measStart      = 1;
model.fieldStoreIncs = round((1:2:60)/60*model.nt)';

%% Save pogo-inp file
% model = rmfield(model,'grain_Orientations');
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO

disp(".pogo-inp saved");

% pic = reshape(model.matTypeRefs(:,1),X_cubic/X_mesh,Y_cubic/Y_mesh)-1;
% figure;
% imshow(pic');

close all;
