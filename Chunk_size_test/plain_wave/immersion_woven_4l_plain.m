% clear;
% clc;
% close all;
% fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

%%
% Display a file selection dialog
% [file_name, file_path] = uigetfile( ...
%     '*.pogo-inp; *.inp', 'Select a  file');
% full_path = fullfile(file_path, file_name);

%% Open the file
% Create the full file path
full_path = "/home/xiaoyu/TexGen_work/orientation_test/090_nogap_noaddheight_20230917.inp";
tic;
data = fx_read_inp_file(full_path);
toc;

%%
% have to define with the same dimention. thus we can reset the coordinate
% by matrix +-.

nx = 601; % adjust accoording to the data
ny = 601;
nz = 440;

dx = 1e-5;
dy = 1e-5;
dz = 1e-5;

% put the xy center to the same center as the nodes form texgen
center = mean(data.nodes);
cx     = center(1);
cy     = center(2);

model1 = genGrid3D(nx, ny, nz, dx, dy, dz, cx, cy, 0);
% avoid the random error from float-point number minus
model1.nodePos = round(model1.nodePos, 5);

model1.elTypeRefs  = ones(length(model1.elNodes(1,:)), 1);
model1.matTypeRefs = ones(length(model1.elNodes(1,:)), 1);

%% calcualte the height for Matrix

yarn_max        = data.nodes(3, 1);
yarn_min        = data.nodes(3, 2);
layer_height_dz = (yarn_max- yarn_min)/dz;

%% asign the matTypeRefs

% the element number of one z location layar
el_num  = size(model1.elTypeRefs, 1);
z_elnum = (nx-1) * (ny-1);

layerup_num = z_elnum * 210; %

for i = 0:7
    model1.matTypeRefs(data.elsets.Matrix + layerup_num + i*z_elnum*layer_height_dz) = 2;
    model1.matTypeRefs(data.elsets.Yarn0  + layerup_num + i*z_elnum*layer_height_dz) = 3;
    model1.matTypeRefs(data.elsets.Yarn1  + layerup_num + i*z_elnum*layer_height_dz) = 3;
    model1.matTypeRefs(data.elsets.Yarn2  + layerup_num + i*z_elnum*layer_height_dz) = 3;
    model1.matTypeRefs(data.elsets.Yarn3  + layerup_num + i*z_elnum*layer_height_dz) = 3;
    model1.matTypeRefs(data.elsets.Yarn4  + layerup_num + i*z_elnum*layer_height_dz) = 4;
    model1.matTypeRefs(data.elsets.Yarn5  + layerup_num + i*z_elnum*layer_height_dz) = 4;
    model1.matTypeRefs(data.elsets.Yarn6  + layerup_num + i*z_elnum*layer_height_dz) = 4;
    model1.matTypeRefs(data.elsets.Yarn7  + layerup_num + i*z_elnum*layer_height_dz) = 4;
end

%% plot the elsets with the model basic dataset
flag_plotelem = 1;

fx_display_model(model1, flag_plotelem);

%%
% Stimulation signal
frequency = 10e6; % unit: Hz
cycles    = 3;
timedelay = 2e-7;
timestep  = 0.5e-9;
endtime   = 1e-6;
phase     = -pi/2;
filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

tof = (19*8 * dz / 2900) / timestep * 2;

%% Settings except for material
model1.prec    = 8;        % Precision
model1.runName = 'Job';
model1.nt      = 10e3;
model1.dt      = timestep;

%% Material settings
% For isotropi49c material
%-----------------------------------------------------%
cWater = 1500; rhoWater = 1000; visc = 0;

% virtual water
model1.matTypes{1,1}.paramsType  = 0;
E = cWater^2 * rhoWater;
model1.matTypes{1,1}.paramValues = [E, 0, rhoWater];

% % For isotropi49c material
% %-----------------------------------------------------%
% model1.matTypes{1,1}.paramsType  = 5; % water media
% model1.matTypes{1,1}.paramValues = [cWater, rhoWater, visc];

model1.matTypes{2,1}.paramsType  = 0;                  % resin-rich interply
model1.matTypes{2,1}.paramValues = [3.7e9, 0.4, 1270];

model1.matTypes{3,1}.paramsType  = 0;                  % CFRP
model1.matTypes{3,1}.paramValues = [13.47e9, 0, 1588]; % virtual parameter!
model1.matTypes{4,1}.paramsType  = 0;                  % CFRP
model1.matTypes{4,1}.paramValues = [13.47e9, 0, 1588]; % virtual parameter!

%-----------------------------------------------------%
% For anisotropic material
%-----------------------------------------------------%
% load anisotropic material properties from the text file
mat_paras = readmatrix([currentFileFolder '/anisotroic_material_prop_0_90.txt']);

% x direction is the main direciton
model1.matTypes{3}.paramsType  = 2; % 2 for anistropic
model1.matTypes{3}.paramValues = mat_paras(1, :);

% pi/4 direction is the main direciton
model1.matTypes{4}.paramsType  = 2; % 2 for anistropic
model1.matTypes{4}.paramValues = mat_paras(2, :);

% % y direction is the main direciton
% model1.matTypes{5}.paramsType  = 2; % 2 for anistropic
% model1.matTypes{5}.paramValues = mat_paras(3, :);

%% Generator
model1.shots{1, 1}.ntSig = length(tb_signal);
model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);
% Plane wave
fd       = 25.6e-3 / 10; % m
diameter = 12.7e-3 / 10;

% Plane wave
z_loc     = min(model1.nodePos(3, :));
center = [    ...
    mean(model1.nodePos(1, :)) ...
    mean(model1.nodePos(2, :)) ...
    z_loc];

dis_to_center = ((model1.nodePos(1, :) - center(1)).^2 + ...
    (model1.nodePos(2, :) - center(2)).^2 + ...
    (model1.nodePos(3, :) - center(3)).^2);

node_index_generator = find(   ...
    model1.nodePos(3, :) == z_loc & ...
    dis_to_center <= (diameter/2).^2);

model1.shots{1, 1}.sigs{1, 1}.sigType    = 1; % 0 - force, 1 - displacement
model1.shots{1, 1}.sigs{1, 1}.isDofGroup = 0;
model1.shots{1, 1}.sigs{1, 1}.dofSpec    = ones(length(node_index_generator),1)*3;
model1.shots{1, 1}.sigs{1, 1}.nodeSpec   = node_index_generator';
model1.shots{1, 1}.sigs{1, 1}.sigAmps    = ones(length(model1.shots{1}.sigs{1}.dofSpec),1)*1e-13;
model1.shots{1, 1}.sigs{1, 1}.sig        = tb_signal(:,2)';

%% Receiver
% node_index_receiver = node_index_generator;
center(3) = z_loc + dz;
% Plane wave
dis_to_center = (model1.nodePos(1, :) - center(1)).^2 + ...
    (model1.nodePos(2, :) - center(2)).^2 + ...
    (model1.nodePos(3, :) - center(3)).^2;

node_index_receiver = find(   ...
    model1.nodePos(3, :) == round(z_loc + dz, 5) & ...
    dis_to_center <= (diameter/2).^2);

model1.measSets{1, 1}.measDof   = repmat((1:model1.nDims)',length(node_index_receiver),1);
model1.measSets{1, 1}.measNodes = reshape(repmat(node_index_receiver,model1.nDims,1),length(node_index_receiver)*model1.nDims,1);

model1.measSets{1, 1}.name       = 'main';
model1.measSets{1, 1}.isDofGroup = 0;
model1.measSets{1, 1}.measDof    = repmat((1:model1.nDims)',length(node_index_receiver),1);
model1.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model1.nDims,1), ...
    length(node_index_receiver)*model1.nDims,1);
model1.measFreq  = 1;
model1.measStart = 1;

node_num               = size(model1.nodePos, 2);
model1.fieldStoreIncs  = round((1:2:60)/60 * model1.nt)';
model1.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';

%%
tic;
% Get the node positions for all elements at once
nodes_pos_all      = model1.nodePos(:, model1.elNodes(:));
% Reshape the array to separate each element's nodes
nodes_pos_reshaped = reshape(nodes_pos_all, size(model1.nodePos, 1), size(model1.elNodes, 1), []);
% Compute the centroids by taking the mean across the second dimension
centroids          = squeeze(mean(nodes_pos_reshaped, 2))';

X = centroids(:, 1);
Y = centroids(:, 2);
Z = centroids(:, 3);

X_len = max(X) - min(X);
Y_len = max(Y) - min(Y);
toc;

%% cut and save
% Save pogo-inp file

% chunk_sizes = [30e-4 26e-4 22e-4 20e-4 18e-4 16.4e-4];
chunk_sizes = [16.4e-4 12.4e-4 8.4e-4 6.4e-4]; % without abs boundary
% chunk_sizes = [6.4e-4 5.4e-4]; % without abs boundary

clc;

for i = 1:numel(chunk_sizes)
    chunk_size = chunk_sizes(i);
    % cut the model
    elsDelete = find(X<cx-chunk_size | X>cx+chunk_size ...
        | Y<cy-chunk_size | Y>cy+chunk_size);

    tic;
    [ mDel ] = deleteEls(model1, elsDelete);
    disp('model cut');
    toc;

    % update center
    center = [    ...
        mean(mDel.nodePos(1, :)) ...
        mean(mDel.nodePos(2, :)) ...
        z_loc];

    % Plane wave
    dis_to_center = (mDel.nodePos(1, :) - center(1)).^2 + ...
        (mDel.nodePos(2, :) - center(2)).^2 + ...
        (mDel.nodePos(3, :) - center(3)).^2;

    node_index_generator = find(   ...
        mDel.nodePos(3, :) == z_loc & ...
        dis_to_center <= (diameter/2).^2);

    mDel.shots{1, 1}.sigs{1, 1}.sigType    = 1; % 0 - force, 1 - displacement
    mDel.shots{1, 1}.sigs{1, 1}.isDofGroup = 0;
    mDel.shots{1, 1}.sigs{1, 1}.dofSpec    = ones(length(node_index_generator),1)*3;
    mDel.shots{1, 1}.sigs{1, 1}.nodeSpec   = node_index_generator';
    mDel.shots{1, 1}.sigs{1, 1}.sigAmps    = ones(length(mDel.shots{1}.sigs{1}.dofSpec),1)*1e-13;
    mDel.shots{1, 1}.sigs{1, 1}.sig        = tb_signal(:,2)';

    % Receiver
    % node_index_receiver = node_index_generator;
    center(3) = z_loc + dz;
    % Plane wave
    dis_to_center = (mDel.nodePos(1, :) - center(1)).^2 + ...
        (mDel.nodePos(2, :) - center(2)).^2 + ...
        (mDel.nodePos(3, :) - center(3)).^2;

    node_index_receiver = find(   ...
        mDel.nodePos(3, :) == round(z_loc + dz, 5) & ...
        dis_to_center <= (diameter/2).^2);

    mDel.measSets{1, 1}.name       = 'main';
    mDel.measSets{1, 1}.isDofGroup = 0;
    mDel.measFreq  = 1;
    mDel.measStart = 2e3;

    mDel.measSets{1, 1}.measDof   = 3 * ones(length(node_index_receiver),1);
    mDel.measSets{1, 1}.measNodes = node_index_receiver;
    node_num             = size(mDel.nodePos, 2);
    mDel.fieldStoreIncs  = round((1:2:60)/60 * mDel.nt)';
    mDel.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';
    % check point 2
    disp(length(node_index_generator));
    disp(length(node_index_receiver));

    %% absorbing boundary
    % Absorbing Regions
    % Only isotropic materials supported for SRM (stiffness reduction method)
    xyz_min = min(mDel.nodePos, [], 2);
    xyz_max = max(mDel.nodePos, [], 2);
    x_min = xyz_min(1);
    y_min = xyz_min(2);
    z_min = xyz_min(3);
    x_max = xyz_max(1);
    y_max = xyz_max(2);
    z_max = xyz_max(3);
    % remove anisotropic boundary
    % ********** define a center out of function for acceleration
    % % Preallocate array for faster execution
    % elementCount = length(mDel.elNodes(1,:));
    % center       = zeros(elementCount, 3); % 3D model
    % % Calculate all centers at once
    % for j = 1: size(mDel.elNodes, 1)
    %     center = center + mDel.nodePos(:, mDel.elNodes(j, :))';
    % end
    % center       = center / size(mDel.elNodes, 1);
    % %
    % cut_size     = 5e-4;
    % X_box_starts = [x_min          x_min          x_min          x_max-cut_size]; % south sides
    % X_box_ends   = [x_min+cut_size x_max          x_max          x_max     ]; % west sides
    % Y_box_starts = [y_min          y_min          y_max-cut_size y_min     ]; % east sides
    % Y_box_ends   = [y_max          y_min+cut_size y_max          y_max     ]; % north sides
    % mDel         = fx_assign_material_to_box(mDel, ...
    %     X_box_starts, Y_box_starts, ones(length(X_box_starts), 1)*z_min, ...
    %     X_box_ends, Y_box_ends, ones(length(X_box_starts), 1)*z_max, 1, center);
    % % review
    % % plot the elsets with the model basic dataset
    % flag_plotelem = 1;
    % fx_display_model(mDel, flag_plotelem);
    % % sides
    % nAbsVals = 50;
    % c0       = 1500;
    % freq     = frequency;
    % abs_size = cut_size;
    % xLims    = [x_min x_min+abs_size x_max-abs_size   x_max];
    % ylims    = [y_min y_min+abs_size y_max-abs_size   y_max];
    % zLims    = [-102  -100           z_max-abs_size   z_max];
    % mDel     = addAbsBound(mDel, xLims, ylims, zLims, nAbsVals, [], [], freq);

    %%
    PogoFilename = [currentFileFolder '/woven_test_chucksize_' num2str(chunk_size*1e6)];
    savePogoInp(sprintf([PogoFilename,'.pogo-inp']), mDel, 1, 15);  % new version POGO
    disp(".pogo-inp saved");

end


%% save the basic model
save([currentFileFolder '/base_model_4d_only090'], 'model1', 'timestep', 'tb_signal', 'dz', 'center', '-v7.3');


