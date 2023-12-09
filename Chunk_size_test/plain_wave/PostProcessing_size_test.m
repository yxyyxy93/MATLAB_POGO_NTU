
clc;
close all;
fclose all;
clear;

%%
% Get the full path to the current script
script_fullpath = mfilename('fullpath');
% Extract the folder from the full path
[script_folder, ~, ~] = fileparts(script_fullpath);
% Define the desired file extension (e.g., '.mat')
file_ext = '*.pogo-hist';
% Get a list of all files with the desired extension in the script's folder
files = dir(fullfile(script_folder, file_ext));

% Extract numbers from the filenames and sort
[~, idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {files.name}));
% Sort the file struct
sortedFiles = files(idx);
% Display the sorted filenames (optional)
disp({sortedFiles.name});

pattern = '(?<=_)(\d+)(?=.p)';
pattern1 = '(?<=0p)(\d+)(?=.p)';
% pattern2 = '(?<=0p)(\d+)(?=um)';
pattern2 = '(?<=0p)(\d+)(?=.p)';

ascan_all = {};
t_all     = {};
idx_cell  = 1;
mesh_all  = [];
% loop through all files and directories
for idx = 1:length(sortedFiles)
    file_name = sortedFiles(idx).name;
    % find the meshsize
    numberStr  = regexp(file_name, pattern, 'match');
    numberStr1 = regexp(file_name, pattern1, 'match');
    numberStr2 = regexp(file_name, pattern2, 'match');
    if ~isempty(numberStr)
        mesh_size = str2double(numberStr{1});
        disp(mesh_size);
    elseif ~isempty(numberStr1)
        mesh_size = str2double(numberStr2{1})/10;
        disp(mesh_size);
    elseif ~isempty(numberStr2)
        mesh_size = str2double(numberStr2{1})/10;
        disp(mesh_size);
    else
        disp('No match found!');
    end
    % 
    % if mesh_size >= 3640 || mesh_size < 1640
    %     continue;
    % end

    mesh_all = [mesh_all mesh_size];
    % read signal
    full_path = fullfile(script_folder, file_name);
    fprintf('Selected file: %s\n', full_path);
    h = loadPogoHist(full_path);
    t = h.startMeas+(0:h.nt-1) * h.dt;

    % only take longitudinal wave signals
    % ascans = h.sets.main.histTraces(:, 3:3:end);
    ascans = h.sets.main.histTraces;
    ascan = mean(ascans, 2);
    inam  = abs(hilbert(ascan));
    inph  = angle(hilbert(ascan));
    
    ascan_all{idx_cell} = ascan;
    t_all{idx_cell}     = t;
    idx_cell = idx_cell + 1;
end

%%
figure,
for idx = 1:size(ascan_all, 2)
    plot(t_all{idx}, ascan_all{idx}, LineWidth=2, DisplayName=['chunk size = ' num2str(mesh_all(idx)/1e3) ' mm']);
    hold on;
end
legend;

xlabel("time (s)");
set(gca, 'linewidth', 2);
set(gca, 'fontname', 'times new roman');
set(gca, 'fontsize', 16);

