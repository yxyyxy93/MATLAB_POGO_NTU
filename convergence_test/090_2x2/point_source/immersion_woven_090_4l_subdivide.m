clear;
clc;
close all;
fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

%% Open the file of the texgen model
% Create the full file path
full_path = "//home/xiaoyu/TexGen_work/convergence_test_model/#090_2x2";

files = dir(full_path);
pattern  = '(?<=_)(\d+)(?=um)';
pattern2 = '(?<=0p)(\d+)(?=um)';

%% read one file
file_name      = files(5).name;
full_file_path = fullfile(full_path, file_name);

% Your code to process each file goes here
numberStr  = regexp(file_name, pattern, 'match');
numberStr2 = regexp(file_name, pattern2, 'match');
if ~isempty(numberStr)
    mesh_size = str2double(numberStr{1});
    disp(mesh_size);
elseif ~isempty(numberStr2)
    mesh_size = str2double(numberStr2{1})/10;
    disp(mesh_size);
else
    disp('No match found!');
end

% read the data
tic;
data = fx_read_inp_file(full_file_path);
toc;

while (mesh_size>=2)
    %% define the basci pogo model
    max_data_node_pos  = data.nodes(:, 1);
    min_data_node_pos  = data.nodes(:, 2);
    diff_data_node_pos = max_data_node_pos - min_data_node_pos;
    nx = round(diff_data_node_pos(1)/mesh_size*1e6 + 1); % adjust accoording to the data
    ny = round(diff_data_node_pos(2)/mesh_size*1e6 + 1);
    nz = 840/mesh_size;
    dx = mesh_size*1e-6;
    dy = mesh_size*1e-6;
    dz = mesh_size*1e-6;
    % put the xy center to the same center as the nodes form texgen
    center = (max_data_node_pos+min_data_node_pos)/2;
    cx     = center(1);
    cy     = center(2);
    model1 = genGrid3D(nx, ny, nz, dx, dy, dz, cx, cy, 0);
    % avoid the random error from float-point number minus
    model1.nodePos     = round(model1.nodePos, 8);
    model1.elTypeRefs  = ones(length(model1.elNodes(1,:)), 1);
    model1.matTypeRefs = ones(length(model1.elNodes(1,:)), 1);

    % calcualte the height for Matrix
    layer_height_dz = (max_data_node_pos(3) - min_data_node_pos(3))/dz;
    % asign the matTypeRefs
    % the element number of one z location layar
    el_num      = size(model1.elTypeRefs, 1);
    z_elnum     = (nx-1) * (ny-1);
    layerup_num = z_elnum * 20/mesh_size; %
    for j = 0:3
        model1.matTypeRefs(data.elsets.Matrix + layerup_num + j*z_elnum*layer_height_dz) = 2;
        model1.matTypeRefs(data.elsets.Yarn0  + layerup_num + j*z_elnum*layer_height_dz) = 3;
        model1.matTypeRefs(data.elsets.Yarn1  + layerup_num + j*z_elnum*layer_height_dz) = 3;
        model1.matTypeRefs(data.elsets.Yarn2  + layerup_num + j*z_elnum*layer_height_dz) = 4;
        model1.matTypeRefs(data.elsets.Yarn3  + layerup_num + j*z_elnum*layer_height_dz) = 4;
    end
    % review
    % plot the elsets with the model basic dataset
    flag_plotelem = 0;
    fx_display_model(model1, flag_plotelem);

    % tof = (19*8 * dz / 2900) / timestep * 2;
    % Settings except for material
    model1.prec    = 8;        % Precision    model1.runName = 'Job';
    model1.nt      = 1e4 /mesh_size*10;
    model1.dt      = 0.25e-9 /10*mesh_size;

    %% Stimulation signal
    frequency = 10e6; % unit: Hz
    cycles    = 3;
    timedelay = 2e-7;
    timestep  = model1.dt;
    endtime   = 3.5e-6;
    phase     = -pi/2;
    filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
    tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

    % Material settings
    % virtual water
    cWater = 1500; rhoWater = 1000;
    model1.matTypes{1,1}.paramsType  = 0;
    E = cWater^2 * rhoWater;
    model1.matTypes{1,1}.paramValues = [E, 0, rhoWater];
    model1.matTypes{2,1}.paramsType  = 0;                  % resin-rich interply
    model1.matTypes{2,1}.paramValues = [3.7e9, 0.4, 1270];
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

    %% Generator
    model1.shots{1, 1}.ntSig = length(tb_signal);
    model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);
    z_loc = min(model1.nodePos(3, :));
    %
    node_index_generator = find(   ...
        round(model1.nodePos(3, :), 8) == round(z_loc + 10/mesh_size*dz, 8) & ...
        model1.nodePos(1, :) == model1.nodePos(1, floor(end/2)) & ...
        model1.nodePos(2, :) == model1.nodePos(2, floor(end/2)));
    % 
    model1.shots{1, 1}.sigs{1, 1}.sigType    = 1; % 0 - force, 1 - displacement
    model1.shots{1, 1}.sigs{1, 1}.isDofGroup = 0;
    model1.shots{1, 1}.sigs{1, 1}.dofSpec    = ones(length(node_index_generator),1)*3;
    model1.shots{1, 1}.sigs{1, 1}.nodeSpec   = node_index_generator';
    model1.shots{1, 1}.sigs{1, 1}.sigAmps    = ones(length(model1.shots{1}.sigs{1}.dofSpec),1)*1e-12 / sqrt(mesh_size);
    model1.shots{1, 1}.sigs{1, 1}.sig        = tb_signal(:,2)';
    % % Receiver
    z_loc = max(model1.nodePos(3, :));
    node_index_receiver = find(   ...
        round(model1.nodePos(3, :), 8) == round(z_loc - 10/mesh_size*dz, 8) & ...
        model1.nodePos(1, :) == model1.nodePos(1, floor(end/2)) & ...
        model1.nodePos(2, :) == model1.nodePos(2, floor(end/2)));
    % node_index_receiver = node_index_generator;
    disp(length(node_index_generator));
    disp(length(node_index_receiver));

    %% Boundary
    xyz_max = max(model1.nodePos, [], 2);
    x_max   = xyz_max(1);
    y_max   = xyz_max(2);
    z_max   = xyz_max(3);

    X_lim_up  = 0;
    X_lim_low = 0;
    Y_lim_up  = y_max;
    Y_lim_low = 0;
    Z_lim_up  = z_max;
    Z_lim_low = 0;
    node_index_generator_yz0 = find(     ...
        (model1.nodePos(1,:)>=X_lim_low)& ...
        (model1.nodePos(1,:)<=X_lim_up)&  ...
        (model1.nodePos(2,:)>=Y_lim_low)& ...
        (model1.nodePos(2,:)<=Y_lim_up)&  ...
        (model1.nodePos(3,:)>=Z_lim_low)& ...
        (model1.nodePos(3,:)<=Z_lim_up));

    X_lim_up  = x_max;
    X_lim_low = x_max;
    node_index_generator_yzx = find(     ...
        (model1.nodePos(1,:)>=X_lim_low)& ...
        (model1.nodePos(1,:)<=X_lim_up)&  ...
        (model1.nodePos(2,:)>=Y_lim_low)& ...
        (model1.nodePos(2,:)<=Y_lim_up)&  ...
        (model1.nodePos(3,:)>=Z_lim_low)& ...
        (model1.nodePos(3,:)<=Z_lim_up));

    X_lim_up  = x_max;
    X_lim_low = 0;
    Y_lim_up  = 0;
    Y_lim_low = 0;
    Z_lim_up  = z_max;
    Z_lim_low = 0;
    node_index_generator_xz0 = find(     ...
        (model1.nodePos(1,:)>=X_lim_low)& ...
        (model1.nodePos(1,:)<=X_lim_up)&  ...
        (model1.nodePos(2,:)>=Y_lim_low)& ...
        (model1.nodePos(2,:)<=Y_lim_up)&  ...
        (model1.nodePos(3,:)>=Z_lim_low)& ...
        (model1.nodePos(3,:)<=Z_lim_up));

    Y_lim_up  = y_max;
    Y_lim_low = y_max;
    node_index_generator_xzy = find(     ...
        (model1.nodePos(1,:)>=X_lim_low)& ...
        (model1.nodePos(1,:)<=X_lim_up)&  ...
        (model1.nodePos(2,:)>=Y_lim_low)& ...
        (model1.nodePos(2,:)<=Y_lim_up)&  ...
        (model1.nodePos(3,:)>=Z_lim_low)& ...
        (model1.nodePos(3,:)<=Z_lim_up));

    model1.fixNodes = [            ...
        node_index_generator_yz0  ...
        node_index_generator_yzx  ...
        node_index_generator_xz0  ...
        node_index_generator_xzy];

    hold on;
    scatter3(...
        model1.nodePos(1, model1.fixNodes), ...
        model1.nodePos(2, model1.fixNodes), ...
        model1.nodePos(3, model1.fixNodes));

    model1.fixDof = [                               ...
        ones(length(node_index_generator_yz0),1)*1; ...
        ones(length(node_index_generator_yzx),1)*1; ...
        ones(length(node_index_generator_xz0),1)*2; ...
        ones(length(node_index_generator_xzy),1)*2];

    %%
    model1.measSets{1, 1}.measDof    = repmat((1:model1.nDims)',length(node_index_receiver),1);
    model1.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model1.nDims,1), ...
        length(node_index_receiver)*model1.nDims,1);
    model1.measSets{1, 1}.name       = 'main';
    model1.measSets{1, 1}.isDofGroup = 0;
    model1.measFreq  = 1;
    model1.measStart = 1;
    %
    node_num               = size(model1.nodePos, 2);
    model1.fieldStoreIncs  = round((1:2:60)/60 * model1.nt)';
    model1.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';
    % save % Save pogo-inp file
    PogoFilename = [currentFileFolder '/' file_name(1:end-4) '_4l_' num2str(mesh_size)];
    savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model1, 1, 15);  % new version POGO
    disp(".pogo-inp saved");

    %% subdivide
    fields = fieldnames(data.elsets); % Get field names for the current struct
    for idx = 2:numel(fields)
        fieldName  = fields{idx};
        fieldValue = data.elsets.(fieldName);
        fprintf('Field: %s, Value numbers: %s\n', fieldName, mat2str(length(fieldValue))); % Use mat2str in case the value is not a string
        newSet = fx_subdivideElementSet(fieldValue);
        fprintf('Field: %s, Value numbers: %s\n', fieldName, mat2str(length(newSet))); % Use mat2str in case the value is not a string
        data.elsets.(fieldName) = newSet;
    end

    mesh_size = mesh_size/2;
    disp(mesh_size);

end


