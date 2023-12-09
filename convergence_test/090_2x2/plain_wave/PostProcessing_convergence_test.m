
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

pattern  = '(?<=l_).*(?=.p)';
pattern1 = '(?<=_).*(?=.p)';
% pattern2 = '(?<=0p)(\d+)(?=um)';
pattern2 = '(?<=0p)(\d+)(?=.p)';

ascan_all1 = {};
ascan_all2 = {};
t_all      = {};
idx_cell   = 1;
mesh_all   = {};
dists      = [];

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
        mesh_size = numberStr1{1};
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
    ascans = h.sets.main.histTraces;
    ascan1  = mean(ascans(:, 1:floor(end/2)), 2);
    ascan2  = mean(ascans(:, ceil(end/2+1),end), 2);

    inam1   = abs(hilbert(ascan1));
    inam2   = abs(hilbert(ascan2));
    
    ascan_all1{idx_cell} = inam1;
    ascan_all2{idx_cell} = inam2;
    t_all{idx_cell}      = t;
    
    % read distance
    full_path = fullfile(script_folder, [file_name(1:end-10) '_dist.mat']);
    dist_file = load(full_path);
    dist = dist_file.dist;
    dists = [dists dist];

    idx_cell = idx_cell + 1;

end

%%
close all;

% windows to select the peak
% time_win1 = [2000 4000];
t1_vec = [];
t2_vec = [];

figure,
for idx = 1:size(ascan_all1, 2)
    plot(t_all{idx}, ascan_all1{idx}, LineWidth=2, DisplayName=['mesh size = ' mesh_all{idx} ' \mum']);
    hold on;
    [pks, locs1] = findpeaks(ascan_all1{idx}, t_all{idx}, "NPeaks", 1, ...
        "MinPeakWidth", 5e-8, "MinPeakHeight", 5e-14);
    plot(locs1, pks, 'o');
    t1_vec = [t1_vec locs1];
end

legend;
xlabel("time (s)");
set(gca, 'linewidth', 2);
set(gca, 'fontname', 'times new roman');
set(gca, 'fontsize', 16);

figure,
for idx = 1:size(ascan_all2, 2)
    plot(t_all{idx}, ascan_all2{idx}, LineWidth=2, DisplayName=['mesh size = ' mesh_all{idx} ' \mum']);
    hold on;
    [pks, locs2] = findpeaks(ascan_all2{idx}, t_all{idx}, "NPeaks", 1, ...
        "MinPeakWidth", 5e-8, "MinPeakHeight", 1e-14);
    plot(locs2, pks, 'o');
    t2_vec = [t2_vec locs2];
end

% legend;
xlabel("time (s)");
set(gca, 'linewidth', 2);
set(gca, 'fontname', 'times new roman');
set(gca, 'fontsize', 16);

figure,
velocities = abs(dists) ./ (t2_vec - t1_vec);
for idx = 1:size(ascan_all2, 2)
    mesh_size = str2double(replace(mesh_all{idx}, 'p', '.'));
    disp(mesh_size);
    plot(mesh_size/250, velocities(idx), 'o', 'linewidth', 2, ...
        DisplayName=['mesh size = ' mesh_all{idx} ' \mum']);
    hold on;
end






