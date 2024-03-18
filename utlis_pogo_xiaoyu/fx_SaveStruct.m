function fx_SaveStruct(currentFileFolder, model_path, nx, ny, dx, dy)
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
