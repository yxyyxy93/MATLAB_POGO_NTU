clear;
clc;
close all;
fclose all;
 
%% Open the file of the texgen model
% Create the full file path
full_path = "/home/xiaoyu/TexGen_work/CrossPlyWoven_45n45090/#(45-45)_10um.inp";
tic;
data_10um = fx_read_inp_file(full_path);
toc;
full_path = "/home/xiaoyu/TexGen_work/CrossPlyWoven_45n45090/#(45-45)_20um.inp";
tic;
data_20um= fx_read_inp_file(full_path);
toc;


%%
mesh_size = 10e-6;

max_data_node_pos = data_20um.nodes(:, 1);
min_data_node_pos = data_20um.nodes(:, 2);
diff_data_node_pos = max_data_node_pos - min_data_node_pos;

nx = round(diff_data_node_pos(1)/mesh_size   + 1); % adjust accoording to the data
ny = round(diff_data_node_pos(2)/mesh_size   + 1);
nz = round(diff_data_node_pos(3)/mesh_size*1 + 1);
dx = mesh_size;
dy = mesh_size;
dz = mesh_size;

center = (max_data_node_pos+min_data_node_pos)/2;
cx     = center(1);
cy     = center(2);
cz     = center(3);

model1 = genGrid3D(nx, ny, nz, dx, dy, dz, cx, cy, cz);

%% subdivide
fields = fieldnames(data_20um.elsets); % Get field names for the current struct
for idx = 2:numel(fields)
    fieldName  = fields{idx};
    fieldValue = data_20um.elsets.(fieldName);
    fprintf('Field: %s, Value numbers: %s\n', fieldName, mat2str(length(fieldValue))); % Use mat2str in case the value is not a string
    newSet = fx_remeshElement(fieldValue, diff_data_node_pos(1)/mesh_size/2, ...
        diff_data_node_pos(2)/mesh_size/2, diff_data_node_pos(3)/mesh_size/2);
    fprintf('Field: %s, Value numbers: %s\n', fieldName, mat2str(length(newSet))); % Use mat2str in case the value is not a string
    data_20um_subdivide.elsets.(fieldName) = int32(sort(newSet));
end

%%
close all;
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_10um.elsets.Yarn0, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_10um.elsets.Yarn1, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_10um.elsets.Yarn2, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_10um.elsets.Yarn3, 0.3);

close all;
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn0, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn1, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn2, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn3, 0.3);
 fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn4, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn5, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn6, 0.3);
fx_plot_element_centroids(model1.nodePos.', ...
    model1.elNodes.', data_20um_subdivide.elsets.Yarn7, 0.3);

%%
%
figure, plot(data_20um_subdivide.elsets.Yarn1);
hold on, plot(data_10um.elsets.Yarn1);
