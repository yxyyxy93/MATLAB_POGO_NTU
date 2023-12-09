% 
clc;
close all;
fclose all;
clear;

% Get the full path of the current file
currentFileFullPath = mfilename('fullpath');
% Get the path to the folder containing this file
[currentFileFolder, ~, ~] = fileparts(currentFileFullPath);
disp(currentFileFolder);

%%
% Set your directory path here
dir_path = './test_battery/';  % Current directory

% Set your file extension here
file_ext = '.pogo-hist';

% Get a list of all files in the directory with the desired file extension
files = dir([dir_path '*' file_ext]);

filenames   = {files.name};
sorted_file = natsortfiles(filenames); 

Ascans_re = nan(length(files), 10000);
Ascans_im = nan(length(files), 10000);

delays = load([dir_path '/focusing_delays.txt']);
% Loop over all files
for i = 1: length(sorted_file)
    % Display the file name
    fprintf('Working on file: %s\n', sorted_file{i});
    h      = loadPogoHist([dir_path sorted_file{i}]);
    ascans = h.sets.main.histTraces; % only y record
    % ************* add delay
    if length(delays) ~= size(ascans, 2)
        error('delays not correct')
    end
    for si = 1:length(delays)
        ascans(:, si) = circshift(ascans(:, si), -round(delays(si)));
    end
    % ********************
    ascan  = mean(ascans, 2);
    Ascans_re(i, 1:length(ascan)) = ascan;
    Ascans_im(i, 1:length(ascan)) = imag(hilbert(ascan));
end

Ascans_re(:, length(ascan)+1:end) = [];
Ascans_im(:, length(ascan)+1:end) = [];

%%
Ascans_inam = (Ascans_re.^2 + Ascans_im.^2).^1/2;

% % in case there are extra lines
% Ascans_inam = Ascans_inam(2:end-1, :);

plot(Ascans_inam.');

%% save with low occupation

Ascans_reim = [Ascans_re; Ascans_im];
% Get the current date
% Convert to a string in the desired format
dateString = datestr(datetime('now'), 'yyyymmdd');
% % save it to a .mat file
% % save('Ascans_4l_orient_20230717_10m', 'Ascans_reim');
save([dir_path, dir_path(2:end-1), '_4l', '_', dateString, '_10m_rect'], 'Ascans_reim');
