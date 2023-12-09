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
 
dir_path = './test_woven_[#45n45#090]_5s/';
dir_path = './test_woven_[#45n45#090]_individual/';

dir_path = './test_woven_[#090]_8/';

% Set your file extension here
file_ext = '.pogo-hist';

% Get a list of all files in the directory with the desired file extension
files = dir([dir_path '*' file_ext]);

filenames   = {files.name};
sorted_file = natsortfiles(filenames); 

Ascans_re = nan(length(files), 10000);
Ascans_im = nan(length(files), 10000);

% Loop over all files
for i = 1: length(sorted_file)
    % Display the file name
    fprintf('Working on file: %s\n', sorted_file{i});
    h      = loadPogoHist([dir_path sorted_file{i}]);
    % ascans = h.sets.main.histTraces(:, 3:3:end); % full xyz record
    ascans = h.sets.main.histTraces; % only z record
    % ********************
    ascan  = mean(ascans, 2);
    Ascans_re(i, 1:length(ascan)) = ascan;
    Ascans_im(i, 1:length(ascan)) = imag(hilbert(ascan));
end

dt = h.dt;
fs = 1/dt;

Ascans_re(:, length(ascan)+1:end) = [];
Ascans_im(:, length(ascan)+1:end) = [];

%%
Ascans_inam = (Ascans_re.^2 + Ascans_im.^2).^1/2;

% % in case there are extra lines
% Ascans_inam = Ascans_inam(2:end-1, :);

% downsample
Ascans_im_ds = Ascans_im(:, 1:int32(fs/2000e6):end);

plot(Ascans_im_ds.');

%% save with low occupation

row = 27;
col = size(Ascans_im_ds, 1) / row;
img_pre = reshape(Ascans_im_ds, [row, col, size(Ascans_im_ds, 2)]);

% Get the current date
% Convert to a string in the desired format
dateString = datestr(datetime('now'), 'yyyymmdd');

% Save size as the first line, then the flattened 3D array
file_name = [dir_path, dir_path(3:end-1), '_15m', '_', dateString, '.csv'];
[x, y, z] = size(img_pre);

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



