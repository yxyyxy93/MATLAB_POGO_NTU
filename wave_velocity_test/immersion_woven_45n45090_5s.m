% clear;
% clc;
% close all;
% fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

%% Stimulation signal
frequency = 15e6; % unit: Hz
cycles    = 3;
timedelay = 1e-7;
timestep  = 0.5e-9;
endtime   = 0.5e-6;
phase     = -pi/2;
filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

%% Open the file of the texgen model
% Create the full file path
full_path = "/home/xiaoyu/TexGen_work/CrossPlyWoven_45n45090/#(45-45)_10um.inp";
tic;
data_45n45 = fx_read_inp_file(full_path);
toc;
full_path = "/home/xiaoyu/TexGen_work/CrossPlyWoven_45n45090/#(090)_10um.inp";
tic;
data_090 = fx_read_inp_file(full_path);
toc;

%%
% define the basci pogo model
mesh_size = 10e-6;

max_data_node_pos = data_45n45.nodes(:, 1);
min_data_node_pos = data_45n45.nodes(:, 2);
% calcualte the height for Matrix
diff_data_node_pos = max_data_node_pos - min_data_node_pos;

%% subdivide
fields_090   = fieldnames(data_090.elsets);
fields_45n45 = fieldnames(data_45n45.elsets);

%%
nx = round(diff_data_node_pos(1)/mesh_size   + 1); % adjust accoording to the data
ny = round(diff_data_node_pos(2)/mesh_size   + 1);
nz = round(diff_data_node_pos(3)/mesh_size*20*2);
dx = mesh_size;
dy = mesh_size;
dz = mesh_size;

center = (max_data_node_pos+min_data_node_pos)/2;
cx     = center(1);
cy     = center(2);
cz     = center(3);

model1 = genGrid3D(nx, ny, nz, dx, dy, dz, cx, cy, cz);

% avoid the random error from float-point number minus
model1.nodePos     = round(model1.nodePos, 5);
model1.elTypeRefs  = ones(length(model1.elNodes(1,:)), 1);
model1.matTypeRefs = ones(length(model1.elNodes(1,:)), 1);

% calcualte the height for Matrix
layer_height_dz = round((max_data_node_pos(3) - min_data_node_pos(3))/dz);

%% asign the matTypeRefs
% the element number of one z location layar
z_elnum = (nx-1) * (ny-1);
layerup_num = z_elnum * 200; %

% [#(-45/45)#(090)]_5
for j = 0:2:9
    model1.matTypeRefs(data_45n45.elsets.(fields_45n45{2})     + layerup_num + j*z_elnum*layer_height_dz) = 2;
    for i = 3:length(fields_45n45)/2+1
        model1.matTypeRefs(data_45n45.elsets.(fields_45n45{i}) + layerup_num + j*z_elnum*layer_height_dz) = 4;
    end
    for i = length(fields_45n45)/2+2:length(fields_45n45)
        model1.matTypeRefs(data_45n45.elsets.(fields_45n45{i}) + layerup_num + j*z_elnum*layer_height_dz) = 6;
    end
    j = j + 1;
    model1.matTypeRefs(data_090.elsets.(fields_090{2}) + layerup_num + j*z_elnum*layer_height_dz) = 2;
    for idx = 3:length(fields_090)/2+1
        model1.matTypeRefs(data_090.elsets.(fields_090{idx}) + layerup_num + j*z_elnum*layer_height_dz) = 3;
    end
    for idx = length(fields_090)/2+2:length(fields_090)
        model1.matTypeRefs(data_090.elsets.(fields_090{idx}) + layerup_num + j*z_elnum*layer_height_dz) = 5;
    end
end

% [#(090)#(-45/45)]_5
for j = 10:2:19
    model1.matTypeRefs(data_090.elsets.(fields_090{2}) + layerup_num + j*z_elnum*layer_height_dz) = 2;
    for idx = 3:length(fields_090)/2+1
        model1.matTypeRefs(data_090.elsets.(fields_090{idx}) + layerup_num + j*z_elnum*layer_height_dz) = 3;
    end
    for idx = length(fields_090)/2+2:length(fields_090)
        model1.matTypeRefs(data_090.elsets.(fields_090{idx}) + layerup_num + j*z_elnum*layer_height_dz) = 5;
    end
    j = j + 1;
    model1.matTypeRefs(data_45n45.elsets.(fields_45n45{2})     + layerup_num + j*z_elnum*layer_height_dz) = 2;
    for i = 3:length(fields_45n45)/2+1
        model1.matTypeRefs(data_45n45.elsets.(fields_45n45{i}) + layerup_num + j*z_elnum*layer_height_dz) = 4;
    end
    for i = length(fields_45n45)/2+2:length(fields_45n45)
        model1.matTypeRefs(data_45n45.elsets.(fields_45n45{i}) + layerup_num + j*z_elnum*layer_height_dz) = 6;
    end
end

% tof = (19*8 * dz / 2900) / timestep * 2;
% Settings except for material
model1.prec = 8;        % Precision    model1.runName = 'Job';
model1.nt   = 15e3;
model1.dt   = timestep;
% % Material settings
% % virtual water
% cWater = 1500; rhoWater = 1000;
% model1.matTypes{1,1}.paramsType  = 0;
% E = cWater^2 * rhoWater;
% model1.matTypes{1,1}.paramValues = [E, 0, rhoWater];
% 

% make water layer the same materials
model1.matTypes{1,1}.paramsType  = 0;                  % resin-rich interply
model1.matTypes{1,1}.paramValues = [3.7e9, 0.4, 1270];

v = sqrt(3.7e9*(1-0.4)/1270/(1+0.4)/(1-2*0.4))
% % For isotropi49c material
% %-----------------------------------------------------%
% model1.matTypes{1,1}.paramsType  = 5; % water media
% model1.matTypes{1,1}.paramValues = [cWater, rhoWater, visc];

model1.matTypes{2,1}.paramsType  = 0;                  % resin-rich interply
model1.matTypes{2,1}.paramValues = [3.7e9, 0.4, 1270];
%-----------------------------------------------------%
% For anisotropic material
%-----------------------------------------------------%
% load anisotropic material properties from the text file
mat_paras = readmatrix([currentFileFolder '/anisotroic_material_prop_0_45_90_135.txt']);

% x direction is the main direciton
model1.matTypes{3}.paramsType  = 2; % 2 for anistropic
model1.matTypes{3}.paramValues = mat_paras(1, :);

% pi/4 direction is the main direciton
model1.matTypes{4}.paramsType  = 2; % 2 for anistropic
model1.matTypes{4}.paramValues = mat_paras(2, :);

% y direction is the main direciton
model1.matTypes{5}.paramsType  = 2; % 2 for anistropic
model1.matTypes{5}.paramValues = mat_paras(3, :);

% 3*pi/4 direction is the main direciton
model1.matTypes{6}.paramsType  = 2; % 2 for anistropic
model1.matTypes{6}.paramValues = mat_paras(4, :);

% Generator
model1.shots{1, 1}.ntSig = length(tb_signal);
model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);
z_loc = min(model1.nodePos(3, :));
node_index_generator = find(   ...
    round(model1.nodePos(3, :), 8) == round(z_loc + 20*dz, 8));
% % plain wave
model1.shots{1, 1}.sigs{1, 1}.sigType    = 1; % 0 - force, 1 - displacement
model1.shots{1, 1}.sigs{1, 1}.isDofGroup = 0;
model1.shots{1, 1}.sigs{1, 1}.dofSpec    = ones(length(node_index_generator),1)*3;
model1.shots{1, 1}.sigs{1, 1}.nodeSpec   = node_index_generator';
model1.shots{1, 1}.sigs{1, 1}.sigAmps    = ones(length(model1.shots{1}.sigs{1}.dofSpec),1)*1e-13; % / sqrt(mesh_size);
model1.shots{1, 1}.sigs{1, 1}.sig        = tb_signal(:,2)';
% % Receiver
node_index_receiver1 = find(   ...
    round(model1.nodePos(3, :), 8) == round(z_loc + 200*dz, 8));
node_index_receiver2 = find(   ...
    round(model1.nodePos(3, :), 8) == round(z_loc + (200 +(j+1)*layer_height_dz)*dz, 8));
disp(length(node_index_generator));
disp(length(node_index_receiver1));
disp(length(node_index_receiver2));
dist = mean(model1.nodePos(3, node_index_receiver1))- ...
    mean(model1.nodePos(3, node_index_receiver2));
disp(['distance between 2 receivers: ' num2str(dist)]);
%
model1.measSets{1, 1}.name       = 'main';
model1.measSets{1, 1}.isDofGroup = 0;
model1.measSets{1, 1}.measDof    = 3 * ones(length([node_index_receiver1 node_index_receiver2]),1);
model1.measSets{1, 1}.measNodes  = [node_index_receiver1 node_index_receiver2];

model1.measFreq  = 1;
model1.measStart = 1;
%
node_num               = size(model1.nodePos, 2);
model1.fieldStoreIncs  = round((1:2:60)/60 * model1.nt)';
model1.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';

%% review
% plot the elsets with the model basic dataset
flag_plotelem = 1;
fx_display_model(model1, flag_plotelem);

%%
tic;
% Get the node positions for all elements at once
nodes_pos_all      = model1.nodePos(:, model1.elNodes(:));
% Reshape the array to separate each element's nodes
nodes_pos_reshaped = reshape(nodes_pos_all, size(model1.nodePos, 1), size(model1.elNodes, 1), []);
% Compute the centroids by taking the mean across the second dimension
centroids          = squeeze(mean(nodes_pos_reshaped, 2))';
% centroids          = squeeze(nodes_pos_reshaped(:,1,:))';

X = int32(centroids(:, 1)*1e6);
Y = int32(centroids(:, 2)*1e6);
Z = int32(centroids(:, 3)*1e6);

X_len = max(X) - min(X);
Y_len = max(Y) - min(Y);
toc;


%% switch to the the scanning point
clc;

chunk_size = 1650;

format long e;
disp(center);
% cut the model
elsDelete = find(X<=int32(center(1)*1e6)-chunk_size | X>=int32(center(1)*1e6)+chunk_size  ...
    | Y<=int32(center(2)*1e6)-chunk_size | Y>=int32(center(2)*1e6)+chunk_size );
tic;
[ mDel ] = deleteEls(model1, elsDelete);
disp(['model cut points:' num2str(length(elsDelete))]);
toc;

%% Boundary
xyz_max = max(mDel.nodePos, [], 2);
x_max   = xyz_max(1);
y_max   = xyz_max(2);
z_max   = xyz_max(3);
xyz_min = min(mDel.nodePos, [], 2);
x_min   = xyz_min(1);
y_min   = xyz_min(2);
z_min   = xyz_min(3);

X_lim_up  = x_min;
X_lim_low = x_min;
Y_lim_up  = y_max;
Y_lim_low = y_min;
Z_lim_up  = z_max;
Z_lim_low = z_min;
node_index_generator_yz0 = find(     ...
    (mDel.nodePos(1,:)>=X_lim_low)& ...
    (mDel.nodePos(1,:)<=X_lim_up)&  ...
    (mDel.nodePos(2,:)>=Y_lim_low)& ...
    (mDel.nodePos(2,:)<=Y_lim_up)&  ...
    (mDel.nodePos(3,:)>=Z_lim_low)& ...
    (mDel.nodePos(3,:)<=Z_lim_up));

X_lim_up  = x_max;
X_lim_low = x_max;
node_index_generator_yzx = find(     ...
    (mDel.nodePos(1,:)>=X_lim_low)& ...
    (mDel.nodePos(1,:)<=X_lim_up)&  ...
    (mDel.nodePos(2,:)>=Y_lim_low)& ...
    (mDel.nodePos(2,:)<=Y_lim_up)&  ...
    (mDel.nodePos(3,:)>=Z_lim_low)& ...
    (mDel.nodePos(3,:)<=Z_lim_up));

X_lim_up  = x_max;
X_lim_low = x_min;
Y_lim_up  = y_min;
Y_lim_low = y_min;
Z_lim_up  = z_max;
Z_lim_low = z_min;
node_index_generator_xz0 = find(     ...
    (mDel.nodePos(1,:)>=X_lim_low)& ...
    (mDel.nodePos(1,:)<=X_lim_up)&  ...
    (mDel.nodePos(2,:)>=Y_lim_low)& ...
    (mDel.nodePos(2,:)<=Y_lim_up)&  ...
    (mDel.nodePos(3,:)>=Z_lim_low)& ...
    (mDel.nodePos(3,:)<=Z_lim_up));

Y_lim_up  = y_max;
Y_lim_low = y_max;
node_index_generator_xzy = find(     ...
    (mDel.nodePos(1,:)>=X_lim_low)& ...
    (mDel.nodePos(1,:)<=X_lim_up)&  ...
    (mDel.nodePos(2,:)>=Y_lim_low)& ...
    (mDel.nodePos(2,:)<=Y_lim_up)&  ...
    (mDel.nodePos(3,:)>=Z_lim_low)& ...
    (mDel.nodePos(3,:)<=Z_lim_up));

mDel.fixNodes = [            ...
    node_index_generator_yz0  ...
    node_index_generator_yzx  ...
    node_index_generator_xz0  ...
    node_index_generator_xzy];

mDel.fixDof = [                               ...
    ones(length(node_index_generator_yz0),1)*1; ...
    ones(length(node_index_generator_yzx),1)*1; ...
    ones(length(node_index_generator_xz0),1)*2; ...
    ones(length(node_index_generator_xzy),1)*2];


%% review
% plot the elsets with the model basic dataset
flag_plotelem = 1;
fx_display_model(mDel, flag_plotelem);


%%
node_num             = size(mDel.nodePos, 2);
mDel.fieldStoreIncs  = round((1:2:60)/60 * mDel.nt)';
mDel.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';
% save % Save pogo-inp file
PogoFilename = [currentFileFolder '/' '_45n45090_'];
save([PogoFilename '_dist.mat'], "dist");
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), mDel, 1, 15);  % new version POGO
disp(".pogo-inp saved");


