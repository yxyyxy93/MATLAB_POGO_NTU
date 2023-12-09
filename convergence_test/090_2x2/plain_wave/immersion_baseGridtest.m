clear;
clc;
close all;
fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

mesh_size = 16;

while (mesh_size>=2)
    %% define the basci pogo model
    nx = 80/mesh_size; % adjust accoording to the data
    ny = 80/mesh_size;
    nz = 6400 /mesh_size;
    dx = mesh_size*1e-6;
    dy = mesh_size*1e-6;
    dz = mesh_size*1e-6;

    model1 = genGrid3D(nx, ny, nz, dx, dy, dz, 0, 0, 0);
    % avoid the random error from float-point number minus
    model1.elTypeRefs  = ones(length(model1.elNodes(1,:)), 1);
    model1.matTypeRefs = ones(length(model1.elNodes(1,:)), 1);

    % tof = (19*8 * dz / 2900) / timestep * 2;
    % Settings except for material
    model1.prec    = 8;        % Precision    model1.runName = 'Job';
    model1.nt      = 30e3 /mesh_size*10;
    model1.dt      = 0.25e-9 /10*mesh_size;

    %% Stimulation signal
    frequency = 10e6; % unit: Hz
    cycles    = 3;
    timedelay = 1e-7;
    timestep  = model1.dt;
    endtime   = 0.5e-6;
    phase     = -pi/2;
    filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
    tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

    % Material settings
    % virtual water
    cWater = 1500; rhoWater = 1000;
    model1.matTypes{1,1}.paramsType  = 0;
    E = cWater^2 * rhoWater;
    model1.matTypes{1,1}.paramValues = [E, 0, rhoWater];
    % Al
    % model1.matTypes{1,1}.paramValues = [69e9, 0.33, 2700];
    v = sqrt(69e9*(1-0.33)/2700/(1+0.33)/(1-2*0.33))

    % Generator
    model1.shots{1, 1}.ntSig = length(tb_signal);
    model1.shots{1, 1}.dtSig = tb_signal(2, 1) - tb_signal(1, 1);
    z_loc = min(model1.nodePos(3, :));
    % Plane wave
    % dis_to_center = (model1.nodePos(1, :) - center(1)).^2 + ...
    %     (model1.nodePos(2, :) - center(2)).^2 + ...
    %     (model1.nodePos(3, :) - center(3)).^2;
    node_index_generator = find(   ...
        round(model1.nodePos(3, :), 8) == round(z_loc + 80/mesh_size*dz, 8));
    % % plain wave
    model1.shots{1, 1}.sigs{1, 1}.sigType    = 1; % 0 - force, 1 - displacement
    model1.shots{1, 1}.sigs{1, 1}.isDofGroup = 0;
    model1.shots{1, 1}.sigs{1, 1}.dofSpec    = ones(length(node_index_generator),1)*3;
    model1.shots{1, 1}.sigs{1, 1}.nodeSpec   = node_index_generator';
    model1.shots{1, 1}.sigs{1, 1}.sigAmps    = ones(length(model1.shots{1}.sigs{1}.dofSpec),1)*1e-13; % / sqrt(mesh_size);
    model1.shots{1, 1}.sigs{1, 1}.sig        = tb_signal(:,2)';
    % % Receiver
    node_index_receiver1 = find(   ...
        round(model1.nodePos(3, :), 8) == round(z_loc + 640/mesh_size*dz, 8));
    z_loc = max(model1.nodePos(3, :));
    node_index_receiver2 = find(   ...
        round(model1.nodePos(3, :), 8) == round(z_loc - 1600/mesh_size*dz, 8));
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

    %% Boundary
    xyz_max = max(model1.nodePos, [], 2);
    x_max   = xyz_max(1);
    y_max   = xyz_max(2);
    z_max   = xyz_max(3);
    xyz_min = min(model1.nodePos, [], 2);
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
    X_lim_low = x_min;
    Y_lim_up  = y_min;
    Y_lim_low = y_min;
    Z_lim_up  = z_max;
    Z_lim_low = z_min;
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
 
    model1.fixDof = [                               ...
        ones(length(node_index_generator_yz0),1)*1; ...
        ones(length(node_index_generator_yzx),1)*1; ...
        ones(length(node_index_generator_xz0),1)*2; ...
        ones(length(node_index_generator_xzy),1)*2];

    %% review
    % plot the elsets with the model basic dataset
    flag_plotelem = 0;
    fx_display_model(model1, flag_plotelem);

    %% Save pogo-inp file
    PogoFilename = [currentFileFolder '/' '3dgrid_' replace(num2str(mesh_size), '.', 'p')];
    save([PogoFilename '_dist.mat'], "dist");
    savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model1, 1, 15);  % new version POGO
    disp(".pogo-inp saved");

    mesh_size = mesh_size/2;
    disp(mesh_size);

end


