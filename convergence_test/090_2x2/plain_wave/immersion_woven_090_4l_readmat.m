% clear;
clc;
close all;
fclose all;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

%% Stimulation signal
frequency = 10e6; % unit: Hz
cycles    = 3;
timedelay = 2e-7;
timestep  = 0.5e-9;
endtime   = 1e-6;
phase     = -pi/2;
filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

% Plane wave
fd       = 25.6e-3 / 20; % m
diameter = 12.7e-3 / 20;

%% Open the file of the texgen model
% Create the full file path

files   = dir(fullfile(currentFileFolder, '*.mat'));

pattern = '(?<=_)(\d+)(?=um)';

%% loop through all files and directories
for i = 1: length(files)
    %
    file_name = files(i).name;
    full_file_path = fullfile(currentFileFolder, file_name);

    load(full_file_path);
    % Your code to process each file goes here
    numberStr = regexp(file_name, pattern, 'match');
    if ~isempty(numberStr)
        mesh_size = str2double(numberStr{1});
        disp(mesh_size);
    else
        disp('No match found!');
    end

    model1.nt      = 6e3 /mesh_size*10;
    model1.dt      = timestep /10*mesh_size;

    % cut the model
    center = [    ...
        mean(model1.nodePos(1, :)) ...
        mean(model1.nodePos(2, :)) ...
        min(model1.nodePos(3, :))];
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
    elsDelete = find(X<center(1)-10e-5 | X>center(1)+10e-5 ...
        | Y<center(2)-10e-5 | Y>center(2)+10e-5);
    tic;
    [ mDel ] = deleteEls(model1, elsDelete);
    disp('model cut');
    toc;

    % save % Save pogo-inp file
    PogoFilename = [currentFileFolder '/' file_name(1:end-4)];
    savePogoInp(sprintf([PogoFilename,'.pogo-inp']), mDel, 1, 15);  % new version POGO
    disp(".pogo-inp saved");
%
end
