clc;
% clear;
close all;
fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

%%
% open the centers file
% centers = dlmread('centers_array.txt');

nx = 20;
ny = 20;
dx = 0.2e-3;
dy = 0.2e-3;
centers = fx_defineloc(dx, dy, nx, ny);
centers = int32(reshape(centers, (nx+1)*(ny+1), 3)*1e6);

figure;
x = centers(:, 1);
y = centers(:, 2);
z = centers(:, 3);

scatter3(x, y, z, 'filled');
% Label each point with its sequence number
for i = 1:length(x)
    text(x(i), y(i), z(i) + 0.2, num2str(i), 'HorizontalAlignment', 'center');
end

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%% read the base model
model_base = load([currentFileFolder '/base_model_20l.mat']);
% model_base = load([currentFileFolder '/base_model_4d_onlyorient_8l.mat']);

timestep  = model_base.timestep;
model1    = model_base.model1;
tb_signal = model_base.tb_signal;
dz        = model_base.dz;

fd       = 25.6e-3 / 5; % m
diameter = 6.35e-3 / 5;

z_loc    = min(model1.nodePos(3, :));

center_ori = int32([    ...
    mean(model1.nodePos(1, :)) ...
    mean(model1.nodePos(2, :)) ...
    z_loc]*1e6);

cWater = 1500;

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

X_unique = unique(X);

X_len = max(X) - min(X);
Y_len = max(Y) - min(Y);
toc;

%% check the centers 
center_list = nan(length(centers), 3);
for step = 1:length(centers)
    center_list(step, :) = center_ori + centers(step, :);
end

scatter3(center_list(:, 1), center_list(:, 2), center_list(:, 3));
hold on;
scatter3([min(X) min(X) max(X) max(X)], [min(Y) max(Y) min(Y) max(Y)], center_list(1:4, 3));

disp(min(center_list(:,1)) - min(X));
disp(min(center_list(:,2)) - min(Y));


model1.prec    = 8;        % Precision
model1.runName = 'Job';
model1.nt      = 15e3;
model1.dt      = timestep;

%% switch to the the scanning point
clc;

chunk_size = 850;
% all unifieid to unit of um

for step = 1:length(centers)
    center = center_ori + centers(step, :);
    disp(step);
    format long e;
    disp(center);
    % cut the model
    elsDelete = find(X<=center(1)-chunk_size | X>=center(1)+chunk_size  ...
        | Y<=center(2)-chunk_size | Y>=center(2)+chunk_size );
    tic;
    [ mDel ] = deleteEls(model1, elsDelete);
    disp(['model cut points:' num2str(length(elsDelete))]);
    toc;

    % update center
    center = round([    ...
    mean(mDel.nodePos(1, :)) ...
    mean(mDel.nodePos(2, :)) ...
    min(mDel.nodePos(3, :))], 5);

    % Plane wave
    dis_to_center = (mDel.nodePos(1, :) - center(1)).^2 + ...
        (mDel.nodePos(2, :) - center(2)).^2 + ...
        (mDel.nodePos(3, :) - center(3)).^2;

    node_index_generator = find(   ...
        mDel.nodePos(3, :) == center(3) & ...
        dis_to_center <= (diameter/2).^2);

    [focused_waves, delays] = fx_focused_wave(fd, center, timestep, ...
        mDel.nodePos(1:3, node_index_generator), cWater, tb_signal(:,2));

    for i = 1: size(focused_waves, 1)
        mDel.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
        mDel.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
        mDel.shots{1, 1}.sigs{i, 1}.sigType    = 1; % 0 - force, 1 - displacement
        mDel.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
        mDel.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
        mDel.shots{1, 1}.sigs{i, 1}.sigAmps    = 1e-13;
    end
    % Receiver
    % node_index_receiver = node_index_generator;
    % Plane wave
    dis_to_center = (mDel.nodePos(1, :) - center(1)).^2 + ...
        (mDel.nodePos(2, :) - center(2)).^2 + ...
        (mDel.nodePos(3, :) - center(3) - dz*2).^2;

    node_index_receiver = find(   ...
        mDel.nodePos(3, :) == center(3)+dz*2 & ...
        dis_to_center <= (diameter/2).^2);

    mDel.measSets{1, 1}.name       = 'main';
    mDel.measSets{1, 1}.isDofGroup = 0;
    mDel.measFreq  = 1;
    mDel.measStart = 1;

    mDel.measSets{1, 1}.measDof   = 3 * ones(length(node_index_receiver),1);
    mDel.measSets{1, 1}.measNodes = node_index_receiver;

    % check point 2
    disp(length(node_index_generator));
    disp(length(node_index_receiver));

    node_num             = size(mDel.nodePos, 2);
    mDel.fieldStoreIncs  = round((1:2:60)/60 * mDel.nt)';
    mDel.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';

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

    % Save pogo-inp file
    PogoFilename = [currentFileFolder '/woven_test_090' num2str(step)];
    savePogoInp(sprintf([PogoFilename,'.pogo-inp']), mDel, 1, 15);  % new version POGO

    disp(".pogo-inp saved");
end


