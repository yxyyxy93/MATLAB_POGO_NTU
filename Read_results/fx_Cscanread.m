function fx_Cscanread(dir_path)
%%
% Set your file extension here
file_ext = '.pogo-hist';

dir_path = strcat(dir_path, '/');
% Get a list of all files in the directory with the desired file extension
files = dir(strcat(dir_path, '*', file_ext));

disp(dir_path);
display(files);

filenames   = {files.name};
sorted_file = natsortfiles(filenames);

Ascans_re = nan(length(files), 10000);
Ascans_im = nan(length(files), 10000);

delays = load(strcat(dir_path, 'focusing_delays.txt'));
% Loop over all files
for i = 1: length(sorted_file)
    % Display the file name
    fprintf('Working on file: %s\n', sorted_file{i});
    h      = loadPogoHist(strcat(dir_path, sorted_file{i}));
    % ascans = h.sets.main.histTraces(:, 3:3:end); % full xyz record
    ascans = h.sets.main.histTraces; % only z record
    % % ************* add delay
    if length(delays) ~= size(ascans, 2)
        error('delays not correct')
    end
    for si = 1:length(delays)
        ascans(:, si) = circshift(ascans(:, si), -round(delays(si))); % 20231211, I compared 3 cases (no dealy, -delay, +delay), -delay show best resutls
    end
    % ********************
    ascan  = mean(ascans, 2);
    Ascans_re(i, 1:length(ascan)) = ascan;
    Ascans_im(i, 1:length(ascan)) = imag(hilbert(ascan));
    % 
end

dt = h.dt;
fs = 1/dt;

Ascans_re(:, length(ascan)+1:end) = [];
% Ascans_im(:, length(ascan)+1:end) = [];
% delete after saving
for i = 1: length(sorted_file)
    delete(strcat(dir_path, sorted_file{i}));
end

% in case of bug 
% save Ascans_re;
% save(strcat(dir_path, "Ascans.mat"), "Ascans_re");

% Define the full path to the Ascans.mat file
filePath = fullfile(dir_path, "Ascans.mat");
% Check if the file exists
if exist(filePath, 'file') == 2
    % File exists, so load Ascans_re from the .mat file
    load(filePath, 'Ascans_re');
    disp('Ascans.mat file found and data loaded.');
else
    % File does not exist
    disp('Ascans.mat file does not exist.');
end

% downsample
Ascans_im_ds = Ascans_re(:, 1:int32(fs/2000e6):end);

%% save with low occupation
row = 17;
col = size(Ascans_im_ds, 1) / row;
img_pre = reshape(Ascans_im_ds, [row, col, size(Ascans_im_ds, 2)]);

% Get the current date
% Convert to a string in the desired format
dateString = datestr(datetime('now'), 'yyyymmdd');

% Save size as the first line, then the flattened 3D array
file_name = strcat(dir_path, '_20m', '_', dateString, '.csv');
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

end
