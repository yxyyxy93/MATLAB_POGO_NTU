function fx_immersion_Cscan(currentFileFolder, model_path, nx, ny, dx, dy, chunk_size)
% open the centers file
% centers = dlmread('centers_array.txt');

centers = fx_defineloc(dx, dy, nx, ny);
centers = int32(reshape(centers, (nx+1)*(ny+1), 3)*1e6);

%% read the base model
try
    model_base = load(model_path);
catch ME
    fprintf('Error loading model: %s\n', ME.message);
    return;  % Exit the function
end
	
% model_base = load([currentFileFolder '/base_model_4d_onlyorient_8l.mat']);
model1    = model_base.model1;
dz        = model_base.dz;

diameter = model_base.diameter;

z_loc    = min(model1.nodePos(3, :));

center_ori = int32([    ...
    mean(model1.nodePos(1, :)) ...
    mean(model1.nodePos(2, :)) ...
    z_loc]*1e6);

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

clear centroids;
clear nodes_pos_all;
clear nodes_pos_reshapes;

% X_unique = unique(X);
% 
% X_len = max(X) - min(X);
% Y_len = max(Y) - min(Y);
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

%% switch to the the scanning point
clc;

% all unifieid to unit of um
for step = 1:length(centers)
    %% skip the existing scanning 
    hist_path = strcat(currentFileFolder,'/woven_test_',num2str(step),'.pogo-hist');
    inp_path  = strcat(currentFileFolder,'/woven_test_090',num2str(step),'.pogo-inp');
    if exist(hist_path, 'file')
        if exist(inp_path, 'file')
            delete(inp_path);
            disp(["deleted file:", inp_path]);
        else
            disp([num2str(step) ": Pogo-hist" "already exist"])
        end
        continue;
    elseif exist(inp_path, 'file')
        disp([num2str(step) ": Pogo-inp" "already exist"])
        continue;
    end

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

    for i = 1: length(node_index_generator)
        mDel.shots{1, 1}.sigs{i, 1}.sigType    = 1; % 0 - force, 1 - displacement
        mDel.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
        mDel.shots{1, 1}.sigs{i, 1}.dofSpec    = 3;
        mDel.shots{1, 1}.sigs{i, 1}.sigAmps    = ones(length(mDel.shots{1}.sigs{i}.dofSpec), 1)*1e-13;
        mDel.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
    end
    % Receiver
    % node_index_receiver = node_index_generator;
    % Plane wave
    dis_to_center = (mDel.nodePos(1, :) - center(1)).^2 + ...
        (mDel.nodePos(2, :) - center(2)).^2 + ...
        (mDel.nodePos(3, :) - center(3)-dz*2).^2;

    node_index_receiver = find(   ...
        mDel.nodePos(3, :) == center(3)+dz*2 & ...
        dis_to_center <= (diameter/2).^2);
		
    mDel.measSets{1, 1}.measDof   = 3 * ones(length(node_index_receiver),1);
    mDel.measSets{1, 1}.measNodes = node_index_receiver;

    % check point 2
    disp(length(node_index_generator));
    disp(length(node_index_receiver));

    node_num             = size(mDel.nodePos, 2);
    mDel.fieldStoreIncs  = round((1:2:60)/60 * mDel.nt)';
    mDel.fieldStoreNodes = round((1:2e6)/2e6 * node_num)';
    % absorbing boundary
    % Absorbing Regions
    % Only isotropic materials supported for SRM (stiffness reduction method)
    % xyz_min = min(mDel.nodePos, [], 2);
    xyz_max = max(mDel.nodePos, [], 2);
    % x_min = xyz_min(1);
    % y_min = xyz_min(2);
    % z_min = xyz_min(3);
    % x_max = xyz_max(1);
    % y_max = xyz_max(2);
    z_max = xyz_max(3);

    % % remove anisotropic boundary
    % sides
    nAbsVals = 60;
    abs_size = 2e-4;
    xLims    = [ ];
    ylims    = [ ];
    zLims    = [-102 -100 z_max-abs_size z_max];
    mDel     = addAbsBound(mDel, xLims, ylims, zLims, nAbsVals, [], [], []);

    % % plot the elsets with the model basic dataset
    % flag_plotelem = 1;
    % fx_display_model(mDel, flag_plotelem);

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
    PogoFilename = strcat(currentFileFolder,'/woven_test_090',num2str(step));
    savePogoInp(sprintf(strcat(PogoFilename,'.pogo-inp')), mDel, 1, 15);  % new version POGO

    disp(".pogo-inp saved");
end

%% saving woven structure
mask_res = X>min(center_list(:, 1)) & X<max(center_list(:, 1)) ...
    & Y>min(center_list(:, 2)) & Y<max(center_list(:, 2)) ...
    & model1.matTypeRefs~=1;

X_del = X(mask_res);
Y_del = Y(mask_res);
Z_del = Z(mask_res);

matetype_del = model1.matTypeRefs(mask_res);

% normalize idx
% there is a ratio to transform the idx to 1, 2, 3, 4 ....
X_del = (X_del - min(X_del))/10 + 1;
Y_del = (Y_del - min(Y_del))/10 + 1;
Z_del = (Z_del - min(Z_del))/10 + 1; % !!!!!!!!!!!!!!! +1 or not there is a doute

x = max(X_del); y = max(Y_del); z = max(Z_del) - 1; % Now I add the -1 term here
% transfer the coordinate to 3d vector
img_pre = zeros(x,y,z);

for i = 1:x*y*z
    img_pre(X_del(i), Y_del(i), Z_del(i)) = matetype_del(i);
end
% reduce the mesh size
img_pre = img_pre(1:10:end, 1:10:end, :);
x = x / 10;
y = y / 10;

% Get the current date
% Convert to a string in the desired format
dateString = datestr(datetime('now'), 'yyyymmdd');

% Save size as the first line, then the flattened 3D array
file_name = [currentFileFolder, '/structure', '_', dateString, '.csv'];
% Open the file for writing
fileID = fopen(file_name, 'w');
% Write the dimensions as a header
fprintf(fileID, '%d,%d,%d\n', x, y, z);
% Iterate through the slices and save them
for i = 1:x
    for j = 1:y
        for k = 1:z
            fprintf(fileID, '%.17f', img_pre(i,j,k)); % Writing each element
            if k < z
                fprintf(fileID, ','); % Separate elements by commas within the same slice
            end
        end
        if j < y
            fprintf(fileID, '\n'); % Newline at the end of each row within a slice
        end
    end
    if i < x
        fprintf(fileID, '\n'); % Extra newline between slices
    end
end
fclose(fileID);

%% save the ultrasound configuration
% Full path for the file
filePath = fullfile(currentFileFolder, 'parameters.txt');
fileID = fopen(filePath, 'w');
% Check if the file was opened successfully
if fileID == -1
    error('File could not be opened');
end
% Write the parameters to the file
fprintf(fileID, 'frequency: %f\n', model_base.frequency);
fprintf(fileID, 'fd: %f\n', model_base.fd);
fprintf(fileID, 'diameter: %f\n', model_base.diameter);
% Close the file
fclose(fileID);

delete(model_path);

end
