function [focused_waves, delays] = fx_focused_wave(fd, center, dt, nodepos, c, signal)
% ************* 
% nodepos: derived from node data from model in pogo
% 1st row x axis, 2nd row y axis
% center, (x, y）

% Calculate distances from each point on the plane to the focal point
distances = sqrt((nodepos(1,:)-center(1)).^2 + (nodepos(2, :)-center(2)).^2 ...
    + (nodepos(3,:) - center(3)).^2 + fd^2);
distances = distances - min(distances);

% % Virtual focusing simulation
focused_waves = zeros(numel(distances), length(signal));  % Initialize focused wave
delays        = zeros(numel(distances), 1);

for i = 1:numel(distances)
    delay = distances(i) / c;  % Calculate time delay for each point on the plane
    delay = round(delay/dt);
    %     disp(delay);
    focused_waves(i, :) = circshift(signal, -delay); % notice: plus or reduce
    delays(i) = delay;
end

% % for debug
% figure,
% for i = 1:20
%     subplot(20, 1, i);
%     plot(focused_waves(i*10 - 9, :));
%     xlim([1200 1800]);
%     hold on;
% end
% 

end