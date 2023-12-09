
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

pattern = '(?<=l_).*(?=.p)';
pattern1 = '(?<=0p)(\d+)(?=.p)';
pattern2 = '(?<=0p)(\d+)(?=.p)';

ascan_all = {};
t_all     = {};
idx_cell  = 1;
mesh_all  = {};

% loop through all files and directories
for idx = 1:length(files)
    file_name = files(idx).name;
    % find the meshsize
    numberStr  = regexp(file_name, pattern, 'match');
    numberStr1 = regexp(file_name, pattern1, 'match');
    numberStr2 = regexp(file_name, pattern2, 'match');
    if ~isempty(numberStr)
        mesh_size = numberStr{1};
        disp(mesh_size);
    elseif ~isempty(numberStr1)
        mesh_size = str2double(numberStr1{1})/100;
        disp(mesh_size);
    elseif ~isempty(numberStr2)
        mesh_size = str2double(numberStr2{1})/100;
        disp(mesh_size);
    else
        disp('No match found!');
    end

    mesh_all{idx} = mesh_size;
    % read signal
    full_path = fullfile(script_folder, file_name);
    fprintf('Selected file: %s\n', full_path);
    h = loadPogoHist(full_path);
    t = h.startMeas+(0:h.nt-1) * h.dt;

    % only take longitudinal wave signals
    ascans = h.sets.main.histTraces(:, 3:3:end);
    % ascans = h.sets.main.histTraces;
    ascan = mean(ascans, 2);
    % ascan = ascans(:, round(end/2));
    inam  = abs(hilbert(ascan));
    inph  = angle(hilbert(ascan));
    
    ascan_all{idx_cell} = ascan;
    t_all{idx_cell}     = t;
    idx_cell = idx_cell + 1;
end

%%
figure,
for idx = 1:size(ascan_all, 2)
    plot(t_all{idx}, ascan_all{idx}, LineWidth=2, DisplayName=['mesh size = ' mesh_all{idx} ' \mum']);
    hold on;
end
legend;

xlabel("time (s)");
set(gca, 'linewidth', 2);
set(gca, 'fontname', 'times new roman');
set(gca, 'fontsize', 16);

