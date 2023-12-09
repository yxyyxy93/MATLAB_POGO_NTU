%
clc;
close all;
fclose all;
clear;

%%
% Set your directory path here
dir_path = './test_woven/';  % Current directory
% dir_path = 'pogo_work/test_woven/';
% dir_path = './test_woven_4l_orient/';  % Current directory
dir_path = './test_woven_4l_shear/';  % Current directory

% Set your file extension here
file_ext = '.pogo-hist';

% Get a list of all files in the directory with the desired file extension
files = dir([dir_path '*' file_ext]);

filenames   = {files.name};
sorted_file = natsortfiles(filenames);

Ascans_re = nan(length(files), 10000);
Ascans_im = nan(length(files), 10000);

delays = load('focusing_delays.txt');
% Loop over all files
for i = 1: length(sorted_file)
    % Display the file name
    fprintf('Working on file: %s\n', sorted_file{i});
    h      = loadPogoHist(sorted_file{i});
    ascans = h.sets.main.histTraces(:, 3:3:end); % full xyz record
    % ascans = h.sets.main.histTraces; % only z record
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

% reshape
Ascans_inam = reshape(Ascans_inam, [9 9 6000]);

% remove the initial signal
Ascans_inam = Ascans_inam(:, :, 2800:6000);

% % display
% orthosliceViewer(Ascans_inam, ...
%     'ScaleFactors', [300 300 1], 'colormap', jet);

% volumeViewer(Ascans_inam);

%% check Bscan

Bscan = squeeze(Ascans_inam(3, :, :)).';
imagesc(Bscan);

%% save with low occupation

Ascans_reim = [Ascans_re; Ascans_im];
% save it to a .mat file
save('Ascans_4l_shear_20230721_10m', 'Ascans_reim');
% save('Ascans_4l_orient_20230721_10m', 'Ascans_reim');


