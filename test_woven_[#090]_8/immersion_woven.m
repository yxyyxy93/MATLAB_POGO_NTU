clear;
clc;
close all;
fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);


%% Open the file
% Create the full file path
full_path = "/home/xiaoyu/TexGen_work/CrossPlyWoven_45n45090/#(090)_10um.inp";
tic;
data_090 = fx_read_inp_file(full_path);
toc;

%%
% define the basci pogo model
mesh_size = 10e-6;

max_data_node_pos = data_090.nodes(:, 1);
min_data_node_pos = data_090.nodes(:, 2);
% calcualte the height for Matrix
diff_data_node_pos = max_data_node_pos - min_data_node_pos;

%% subdivide
fields_090   = fieldnames(data_090.elsets);

% for idx = 2:numel(fields_090)
%     fieldName  = fields_090{idx};
%     fieldValue = data_090.elsets.(fieldName);
%     fprintf('Field: %s, Value numbers: %s\n', fieldName, mat2str(length(fieldValue))); % Use mat2str in case the value is not a string
%     newSet = fx_remeshElement(fieldValue, diff_data_node_pos(1)/mesh_size, ...
%         diff_data_node_pos(2)/mesh_size, diff_data_node_pos(3)/mesh_size);
%     fprintf('Field: %s, Value numbers: %s\n', fieldName, mat2str(length(newSet))); % Use mat2str in case the value is not a string
%     data_090.elsets.(fieldName) = newSet;
% end
%
% mesh_size = mesh_size/2;

%%
nx = round(diff_data_node_pos(1)/mesh_size   + 1); % adjust accoording to the data
ny = round(diff_data_node_pos(2)/mesh_size   + 1);
nz = round(diff_data_node_pos(3)/mesh_size*8*2);
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
layerup_num = z_elnum * 140; %

for seed = 21:30
    rng(seed);
    % reset material
    model1.matTypeRefs = ones(length(model1.elNodes(1,:)), 1);
    numX      = diff_data_node_pos(1)/mesh_size;
    numY      = diff_data_node_pos(2)/mesh_size;
    numZ      = diff_data_node_pos(3)/mesh_size;

    for j = 0:1:7
        % shift the layer
        x_shift = rand();
        y_shift = rand();
        for idx = 2:length(fields_090)
            indices   = data_090.elsets.(fields_090{idx});
            [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, numZ);
            x = mod((x + round(numX*x_shift)), numX);
            y = mod((y + round(numY*y_shift)), numY);
            x(x==0) = numX; % last one should no be zero
            y(y==0) = numY;
            if (idx==2)
                model1.matTypeRefs(coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz) = 2;
            elseif (idx>2 && idx<=length(fields_090)/2+1)
                model1.matTypeRefs(coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz) = 3;
            elseif (idx>length(fields_090)/2+1)
                model1.matTypeRefs(coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz) = 4;
            end
          
        end
    end


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
    model1.nt      = 8e3;
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
    %-----------------------------------------------------%
    % For anisotropic material
    %-----------------------------------------------------%
    % load anisotropic material properties from the text file
    mat_paras = readmatrix([currentFileFolder '/anisotroic_material_prop_0_45_90_135.txt']);

    % x direction is the main direciton
    model1.matTypes{3}.paramsType  = 2; % 2 for anistropic
    model1.matTypes{3}.paramValues = mat_paras(1, :);

    % y direction is the main direciton
    model1.matTypes{4}.paramsType  = 2; % 2 for anistropic
    model1.matTypes{4}.paramValues = mat_paras(3, :);


    %% Generator
    model1.shots{1, 1}.ntSig = length(tb_signal);
    model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);

    % Plane wave
    fd       = 25.6e-3 / 10; % m
    diameter = 12.7e-3 / 10;

    z_loc       = min(model1.nodePos(3, :));

    center = [    ...
        mean(model1.nodePos(1, :)) ...
        mean(model1.nodePos(2, :)) ...
        z_loc];

    dis_to_center = (model1.nodePos(1, :) - center(1)).^2 + ...
        (model1.nodePos(2, :) - center(2)).^2 + ...
        (model1.nodePos(3, :) - center(3)).^2;

    node_index_generator = find(   ...
        model1.nodePos(3, :) == z_loc & ...
        dis_to_center <= (diameter/2).^2);

    [focused_waves, delays] = fx_focused_wave(fd, center, timestep, ...
        model1.nodePos(1:3, node_index_generator), cWater, tb_signal(:,2));

    for i = 1: size(focused_waves, 1)
        model1.shots{1, 1}.sigs{i, 1}.sigType    = 1; % 0 - force, 1 - displacement
        model1.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
        model1.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
        model1.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
        model1.shots{1, 1}.sigs{i, 1}.sigAmps    = ones(length(model1.shots{1}.sigs{i}.dofSpec), 1)*1e-13;
        model1.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
    end

    save([currentFileFolder '/focusing_delays.txt'], 'delays', '-ascii');

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

    model1.measSets{1, 1}.measDof    = repmat((1:model1.nDims)',length(node_index_receiver),1);
    model1.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model1.nDims,1),length(node_index_receiver)*model1.nDims,1);
    model1.measSets{1, 1}.name       = 'main';
    model1.measSets{1, 1}.isDofGroup = 0;
    model1.measSets{1, 1}.measDof    = repmat((1:model1.nDims)',length(node_index_receiver),1);
    model1.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model1.nDims,1),length(node_index_receiver)*model1.nDims,1);
    model1.measFreq        = 1;
    model1.measStart       = 1;

    node_num               = size(model1.nodePos, 2);
    model1.fieldStoreIncs  = round((1:2:60)/60 * model1.nt)';
    model1.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';

    % fill in the "water gap" materials
    indices   = find(model1.matTypeRefs==2);
    first_idx = indices(1);
    last_idx  = indices(end);
    indices   = first_idx + find(model1.matTypeRefs(first_idx:last_idx)==1);
    model1.matTypeRefs(indices) = 2;

    % plot the elsets with the model basic dataset
    flag_plotelem = 0;
    fx_display_model(model1, flag_plotelem);

    %% save
    % Save pogo-inp file

    % PogoFilename = [currentFileFolder '/woven_test_base'];
    % savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model1, 1, 15);  % new version POGO
    % disp(".pogo-inp saved");

    %% save the basic model
    save([currentFileFolder '/base_model_shiftseed_' num2str(seed)], 'model1', 'timestep', 'tb_signal', 'dz', 'center', '-v7.3');
end

% for seed = 123
%     mkdir([currentFileFolder '/base_model_shiftseed_' num2str(seed)]);
% end
