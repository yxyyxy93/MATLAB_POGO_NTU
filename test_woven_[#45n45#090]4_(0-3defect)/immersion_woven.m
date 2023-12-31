clear;
clc;
close all;
fclose all;

addpath(pathdef);

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);


%% Open the file
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

max_data_node_pos = data_090.nodes(:, 1);
min_data_node_pos = data_090.nodes(:, 2);
% calcualte the height for Matrix
diff_data_node_pos = max_data_node_pos - min_data_node_pos;

%% subdivide
fields_090   = fieldnames(data_090.elsets);
fields_45n45 = fieldnames(data_45n45.elsets);

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
% ********** define a center out of function for acceleration
% Preallocate array for faster execution
elementCount = length(model1.elNodes(1,:));
center_elem  = zeros(elementCount, 3); % 3D model

% Calculate all centers at once
for i = 1: size(model1.elNodes, 1)
    center_elem = center_elem + model1.nodePos(:, model1.elNodes(i, :))';
end

center_elem = center_elem / size(model1.elNodes, 1);

%% asign the matTypeRefs
% the element number of one z location layar
z_elnum = (nx-1) * (ny-1);
layerup_num = z_elnum * 140; %

numX    = diff_data_node_pos(1)/mesh_size;
numY    = diff_data_node_pos(2)/mesh_size;
numZ    = diff_data_node_pos(3)/mesh_size;

for seed = 5:8
    rng(seed);
    % Generate a random integer number between 0 and 3
    num_defects = randi([0, 3], 1, 1);
    % Create an instance of the class with the random number of defects
    defects = random_defects(num_defects);
    % reset material
    model1.matTypeRefs = ones(length(model1.elNodes(1,:)), 1);
    for j = 0:2:3
        % Add shift logic
        x_shift = rand();
        y_shift = rand();
        % Loop for data_45n45
        for i = 2:length(fields_45n45)
            indices = data_45n45.elsets.(fields_45n45{i});
            [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, numZ);
            x = mod((x + round(numX*x_shift)), numX)+1;
            y = mod((y + round(numY*y_shift)), numY)+1;
            newIndices = coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz;
            if (i == 2)
                model1.matTypeRefs(newIndices) = 2;
            elseif (i <= length(fields_45n45)/2+1)
                model1.matTypeRefs(newIndices) = 4;
            else
                model1.matTypeRefs(newIndices) = 6;
            end
        end
        x_shift = rand();
        y_shift = rand();
        j = j + 1;
        % Loop for data_090
        for i = 2:length(fields_090)
            indices = data_090.elsets.(fields_090{i});
            [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, numZ);
            x = mod((x + round(numX*x_shift)), numX)+1;
            y = mod((y + round(numY*y_shift)), numY)+1;
            newIndices = coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz;
            if (i == 2)
                model1.matTypeRefs(newIndices) = 2;
            elseif (i <= length(fields_090)/2+1)
                model1.matTypeRefs(newIndices) = 3;
            else
                model1.matTypeRefs(newIndices) = 5;
            end
        end
    end
    for j = 4:2:7
        % Add shift logic
        x_shift = rand();
        y_shift = rand();
        % Loop for data_090
        for i = 2:length(fields_090)
            indices = data_090.elsets.(fields_090{i});
            [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, numZ);
            x = mod((x + round(numX*x_shift)), numX)+1;
            y = mod((y + round(numY*y_shift)), numY)+1;
            newIndices = coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz;
            if (i == 2)
                model1.matTypeRefs(newIndices) = 2;
            elseif (i <= length(fields_090)/2+1)
                model1.matTypeRefs(newIndices) = 3;
            else
                model1.matTypeRefs(newIndices) = 5;
            end
        end
        x_shift = rand();
        y_shift = rand();
        j = j + 1;
        % Loop for data_45n45
        for i = 2:length(fields_45n45)
            indices = data_45n45.elsets.(fields_45n45{i});
            [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, numZ);
            x = mod((x + round(numX*x_shift)), numX)+1;
            y = mod((y + round(numY*y_shift)), numY)+1;
            newIndices = coordinatesToLinearIndex(x, y, z, numX, numY, numZ) + layerup_num + j*z_elnum*layer_height_dz;
            if (i == 2)
                model1.matTypeRefs(newIndices) = 2;
            elseif (i <= length(fields_45n45)/2+1)
                model1.matTypeRefs(newIndices) = 4;
            else
                model1.matTypeRefs(newIndices) = 6;
            end
        end
    end
    % ********** assign defects
    % Generate and plot circles
    defects = defects.generateCircles();
    % defects.plotCircles();
    % Display the generated circles' parameters
    x         = defects.circles(:, 1) * 1e-3; % mm -> m
    y         = defects.circles(:, 2) * 1e-3;
    Radius    = defects.circles(:, 3) * 1e-3; % mm -> m
    Depth     = defects.circles(:, 4) * layer_height_dz * dz; % idx -> m
    Thickness = defects.circles(:, 5) * 1e-6; % um -> m
    disp([x y Radius Depth Thickness]);
    model1 = fx_assign_material_to_cylinder(model1, ...
         Depth-Thickness, Depth+Thickness, Radius, 7, center_elem, [x y]);

    %***********  % Stimulation signal
    frequency = 15e6; % unit: Hz
    cycles    = 3;
    timedelay = 2e-7;
    timestep  = 1e-9;
    endtime   = 1e-6;
    phase     = -pi/2;
    filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
    tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

    tof = (19*8 * dz / 2900) / timestep * 2;

    %***************** Settings except for material
    model1.prec    = 8;        % Precision
    model1.runName = 'Job';
    model1.nt      = 8e3;
    model1.dt      = timestep;

    %***************** Material settings
    % For isotropi49c material
    %-----------------------------------------------------%
    % from paper "Towards using convolutional neural network to locate ..."
    E = 0.04e9;  % Young's modulus in Pascals (Pa)
    nu = 0.497;  % Poisson's ratio (dimensionless)
    rhoWater = 1000;
    % Calculate the speed of longitudinal waves
    c_L = sqrt(E * (1 - nu) / (rhoWater * (1 + nu) * (1 - 2 * nu)));
    % Calculate the speed of shear waves
    c_S = sqrt(E / (2 * rhoWater * (1 + nu)));
    % virtual water
    model1.matTypes{1,1}.paramsType  = 0;
    model1.matTypes{1,1}.paramValues = [0.04e9, 0.497, rhoWater];

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
    % defect --- virtual air
    model1.matTypes{7,1}.paramsType  = 0;
    model1.matTypes{7,1}.paramValues = [1000, 0.497, 1.225];

    %***************** Generator
    model1.shots{1, 1}.ntSig = length(tb_signal);
    model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);
    % Plane wave
    fd       = 25.6e-3 / 10; % m
    diameter = 6.35e-3 / 10;
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
        model1.nodePos(1:3, node_index_generator), c_L, tb_signal(:,2));

    for i = 1: size(focused_waves, 1)
        model1.shots{1, 1}.sigs{i, 1}.sigType    = 1; % 0 - force, 1 - displacement
        model1.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
        model1.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
        model1.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
        model1.shots{1, 1}.sigs{i, 1}.sigAmps    = ones(length(model1.shots{1}.sigs{i}.dofSpec), 1)*1e-13;
        model1.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
    end

    save([currentFileFolder '/focusing_delays.txt'], 'delays', '-ascii');

    %***************** Receiver
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

    %     % plot the elsets with the model basic dataset
    %     flag_plotelem = 1;
    %     fx_display_model(model1, flag_plotelem);

    %***************** save the basic model
    save([currentFileFolder '/base_model_shiftseed_' num2str(seed)], 'model1', 'timestep', 'tb_signal', 'dz', 'center', 'fd', 'diameter', 'frequency', '-v7.3');
    %***************** save
    % Open the text file for appending
    fileID = fopen('defects_log.txt', 'a');
    % Write the set name and defect parameters to the file
    fprintf(fileID, 'Set Name: %s\n', ['shiftseed_' num2str(seed)]);
    % Write the number of defects and their parameters to the file
    fprintf(fileID, 'Num Defects: %d\n', num_defects);
    for i = 1:num_defects
        fprintf(fileID, 'Defect %d: X = %f, Y = %f, Radius = %f, Depth = %f, Thickness = %f\n', ...
            i, x(i), y(i), Radius(i), Depth(i), Thickness(i));
    end
    % Close the file
    fclose(fileID);
    % Save pogo-inp file
    % PogoFilename = [currentFileFolder '/woven_test_base'];
    % savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model1, 1, 15);  % new version POGO
    % disp(".pogo-inp saved");
    fx_display_model(model1, 0);
end
 