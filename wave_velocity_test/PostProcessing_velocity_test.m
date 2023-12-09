
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

file_name= files(1).name;

% read distance
full_path = fullfile(script_folder, [file_name(1:end-10) '_dist.mat']);
dist_file = load(full_path);
dist = dist_file.dist;

% read signal
full_path = fullfile(script_folder, files(1).name);
fprintf('Selected file: %s\n', full_path);
h = loadPogoHist(full_path);
t = h.startMeas+(0:h.nt-1) * h.dt;

% only take longitudinal wave signals
ascans = h.sets.main.histTraces;
ascan1  = mean(ascans(:, 1:floor(end/2)), 2);
ascan2  = mean(ascans(:, ceil(end/2+1),end), 2);

inam1   = abs(hilbert(ascan1));
inam2   = abs(hilbert(ascan2));

close all;

figure,
plot(t, inam1, LineWidth=2);
hold on;
[pks, locs1] = findpeaks(inam1, t, "NPeaks", 1, ...
    "MinPeakWidth", 5e-8, "MinPeakHeight", 5e-14);
plot(locs1, pks, 'o');


legend;
xlabel("time (s)");
set(gca, 'linewidth', 2);
set(gca, 'fontname', 'times new roman');
set(gca, 'fontsize', 16);

figure,

plot(t, inam2, LineWidth=2);
hold on;
[pks, locs2] = findpeaks(inam2, t, "NPeaks", 1, ...
    "MinPeakWidth", 5e-8, "MinPeakHeight", 1e-14);
plot(locs2, pks, 'o');

% legend;
xlabel("time (s)");
set(gca, 'linewidth', 2);
set(gca, 'fontname', 'times new roman');
set(gca, 'fontsize', 16);

velocity = abs(dist) ./ (locs2 - locs1)






