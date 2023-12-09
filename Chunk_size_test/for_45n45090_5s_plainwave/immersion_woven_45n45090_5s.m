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
nz = round(diff_data_node_pos(3)/mesh_size*20 + 160);
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
layerup_num = z_elnum * 80; %

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
model1.nt   = 10e3;
model1.dt   = timestep;
% % Material settings
% virtual water
cWater = 1500; rhoWater = 1000;
model1.matTypes{1,1}.paramsType  = 0;
E = cWater^2 * rhoWater;
model1.matTypes{1,1}.paramValues = [E, 0, rhoWater];

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

% change to isotropic materials
model1.matTypes{3}.paramsType  = 0;
model1.matTypes{3}.paramValues = [7.74e9, 0.36, 1588];
model1.matTypes{4}.paramsType  = 0;
model1.matTypes{4}.paramValues = [7.74e9, 0.36, 1588];
model1.matTypes{5}.paramsType  = 0;
model1.matTypes{5}.paramValues = [7.74e9, 0.36, 1588];
model1.matTypes{6}.paramsType  = 0;
model1.matTypes{6}.paramValues = [7.74e9, 0.36, 1588];

sqrt(E_equal*(1-nu_equal)/1588/(1+nu_equal)/(1-2*nu_equal))

%
node_num               = size(model1.nodePos, 2);
model1.fieldStoreIncs  = round((1:2:60)/60 * model1.nt)';
model1.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';

%% review
% plot the elsets with the model basic dataset
% flag_plotelem = 1;
% fx_display_model(model1, flag_plotelem);

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

%% Plane wave
fd       = 25.6e-3; % m
diameter = 6.35e-3/15;

% update center
z_loc = min(model1.nodePos(3, :)) + 30*dz;

center = [    ...
    mean(model1.nodePos(1, :)) - 1e-3...
    mean(model1.nodePos(2, :)) - 1e-3 ...
    z_loc];

% Plane wave
dis_to_center = (model1.nodePos(1, :) - center(1)).^2 + ...
    (model1.nodePos(2, :) - center(2)).^2 + ...
    (model1.nodePos(3, :) - center(3)).^2;
node_index_generator = find(   ...
    int32(model1.nodePos(3, :)*1e6) == int32(z_loc*1e6) & ...
    dis_to_center <= (diameter/2).^2);
[focused_waves, delays] = fx_focused_wave(fd, center, timestep, ...
    model1.nodePos(1:3, node_index_generator), cWater, tb_signal(:,2));
for i = 1: size(focused_waves, 1)
    model1.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
    model1.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
    model1.shots{1, 1}.sigs{i, 1}.sigType    = 1; % 0 - force, 1 - displacement
    model1.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
    model1.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
    model1.shots{1, 1}.sigs{i, 1}.sigAmps    = 1e-13;
end

model1.shots{1, 1}.ntSig = length(tb_signal);
model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);

% Receiver
% node_index_receiver = node_index_generator;
center(3) = z_loc + dz;
% Plane wave
dis_to_center = (model1.nodePos(1, :) - center(1)).^2 + ...
    (model1.nodePos(2, :) - center(2)).^2 + ...
    (model1.nodePos(3, :) - center(3)).^2;

node_index_receiver = find(   ...
    model1.nodePos(3, :) == round(z_loc + dz, 5) & ...
    dis_to_center <= (diameter/2).^2);

model1.measSets{1, 1}.name       = 'main';
model1.measSets{1, 1}.isDofGroup = 0;
model1.measFreq  = 1;
model1.measStart = 0.5e3;

model1.measSets{1, 1}.measDof   = 3 * ones(length(node_index_receiver),1);
model1.measSets{1, 1}.measNodes = node_index_receiver;
node_num               = size(model1.nodePos, 2);
model1.fieldStoreIncs  = round((1:2:60)/60 * model1.nt)';
model1.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';

%% cut and save
% Save pogo-inp file
 
chunk_sizes = [15 10 8 7] * 1e2; % without abs boundary

clc;

for i = 1:numel(chunk_sizes)
    chunk_size = chunk_sizes(i);
    % cut the model
    elsDelete = find(X<int32(center(1)*1e6)-chunk_size | X>int32(center(1)*1e6)+chunk_size ...
        | Y<int32(center(2)*1e6)-chunk_size | Y>int32(center(2)*1e6)+chunk_size);

    tic;
    [ mDel ] = deleteEls(model1, elsDelete);
    disp('model cut');
    toc;

    % check point 2
    disp(length(node_index_generator));
    disp(length(node_index_receiver));

    % absorbing boundary
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

    % % remove anisotropic boundary
    % ********** define a center out of function for acceleration
    % Preallocate array for faster execution
    elementCount = length(mDel.elNodes(1,:));
    node_center  = zeros(elementCount, 3); % 3D model
    % Calculate all centers at once
    for j = 1: size(mDel.elNodes, 1)
        node_center = node_center + mDel.nodePos(:, mDel.elNodes(j, :))';
    end
    node_center = node_center / size(mDel.elNodes, 1);
    %
    cut_size     = 1e-4;
    X_box_starts = [x_min          x_min          x_min          x_max-cut_size]; % south sides
    X_box_ends   = [x_min+cut_size x_max          x_max          x_max         ]; % west sides
    Y_box_starts = [y_min          y_min          y_max-cut_size y_min         ]; % east sides
    Y_box_ends   = [y_max          y_min+cut_size y_max          y_max         ]; % north sides
    % mDel         = fx_assign_material_to_box(mDel, ...
    %     X_box_starts, Y_box_starts, ones(length(X_box_starts), 1)*z_min, ...
    %     X_box_ends, Y_box_ends, ones(length(X_box_starts), 1)*z_max, 1, node_center);

    % mDel         = fx_assign_material_to_box(mDel, ...
    %     X_box_starts, Y_box_starts, ones(length(X_box_starts), 1)*(z_min+80*dz), ...
    %     X_box_ends, Y_box_ends, ones(length(X_box_starts), 1)*(z_min+80*dz), 2, node_center);

    % review
    % plot the elsets with the model basic dataset
    flag_plotelem = 1;
    fx_display_model(mDel, flag_plotelem);

    % sides
    freq     = frequency;
    nAbsVals = 50;
    abs_size = cut_size;
    xLims    = [x_min x_min+abs_size x_max-abs_size x_max];
    ylims    = [y_min y_min+abs_size y_max-abs_size y_max];
    % abs_size = 1e-4;
    % xLims    = [ ];
    % ylims    = [ ];
    zLims    = [z_min z_min+2*abs_size z_max-2*abs_size z_max];
    mDel     = addAbsBound(mDel, xLims, ylims, zLims, nAbsVals, [], [], freq);

    % Boundary
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

    PogoFilename = [currentFileFolder '/woven_test_chucksize_' num2str(chunk_size*1e6)];
    savePogoInp(sprintf([PogoFilename,'.pogo-inp']), mDel, 1, 15);  % new version POGO
    disp(".pogo-inp saved");

end

%% save the basic model
save([currentFileFolder '/base_model_45n45090'], 'model1', 'timestep', 'tb_signal', 'dz', 'center', '-v7.3');


