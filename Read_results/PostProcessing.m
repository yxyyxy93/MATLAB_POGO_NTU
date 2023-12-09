% 
% clc;
% close all;
% fclose all;
% clear;

% Display a file selection dialog
[file_name, file_path] = uigetfile( ...
    '*.pogo-hist; *.pogo-field', 'Select a  file');

% Check if a file was selected
if isequal(file_name, 0)
    disp('No file was selected.');
else
    % Create the full file path
    full_path = fullfile(file_path, file_name);

    % Get the file parts (path, name, and extension)
    [~, ~, file_ext] = fileparts(full_path);

    % Display the selected file and its extension
    fprintf('Selected file: %s\n', full_path);
    fprintf('File extension: %s\n', file_ext);

    % Perform actions based on the file extension
    switch file_ext
        case '.pogo-field'
            f = loadPogoField(full_path);
        case '.pogo-hist'
            h = loadPogoHist(full_path);
        otherwise
            disp('Unknown file type.');
            % Add code to handle unknown file types
    end

end

%%
t = h.startMeas+(0:h.nt-1) * h.dt;

% only take longitudinal wave signals
% ascans = h.sets.main.histTraces(:, 3:3:end);
ascans = h.sets.main.histTraces;

% ************* add delay
if exist(fullfile(file_path, 'focusing_delays.txt'), 'file')
    delays = load(fullfile(file_path, 'focusing_delays.txt'));
    if length(delays) ~= size(ascans, 2)
        error('delays not correct')
    end

    ascans_defoc = ascans;
    for i = 1:length(delays)
        ascans_defoc(:, i) = circshift(ascans(:, i), -round(delays(i)));
    end
    % ascans = ascans_defoc;
    ascan = mean(ascans_defoc, 2);
else
    ascan = mean(ascans, 2);
end

inam  = abs(hilbert(ascan));

inph  = angle(hilbert(ascan));

figure,
subplot(2, 1, 1);
plot(t, ascan);
hold on;
plot(t, inam);
hold on;
% xlim([2.5e-5 4.4e-5]);

subplot(2, 1, 2);
plot(t, inph, "Color", 'magenta');
hold on;

% scatter(t(I+5000), inph(I+5000), 'red', LineWidth=2);
% xlim([2.5e-5 4.4e-5]);

%%
figure,
plot(t, mean(ascans, 2));
% hold on;
% plot(t, mean(ascans_defoc, 2));



