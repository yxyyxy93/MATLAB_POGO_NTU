function fx_plot_element_centroids(nodes, elements, elsets, ~)
% nodes is the entire nodes data structure
% elements is the entire elements data structure

figure(1);

dims = size(nodes, 2);

% check if 2D
if dims==2
    centroids = zeros(size(elsets, 2), 2);  % Preallocate centroid array

    for i = 1:1:size(elsets, 2) % skip by ? to save time
        element_nodes   = elements(elsets(:, i), 1:end); % the first colume includeds
        nodes_pos       = nodes(element_nodes, 1:end);
        centroids(i, :) = mean(nodes_pos, 1);
    end
    % skip to save memory
    X = centroids(1:1:end, 1);
    Y = centroids(1:1:end, 2);
    % scatter3(X, Y, Z, 'filled', 'MarkerFaceAlpha', trans);
    plot(X, Y, '.', 'MarkerSize', 1);
    hold on;
    % check 3D or higher
else
    centroids = zeros(size(elsets, 2), 3);  % Preallocate centroid array

    for i = 1:10:size(elsets, 2) % skip by 10 to save time
        element_nodes   = elements(elsets(:, i), 1:end); % the first colume includeds
        nodes_pos       = nodes(element_nodes, 1:end);
        centroids(i, :) = mean(nodes_pos, 1);
    end

    % skip to save memory
    X = centroids(1:10:end, 1);
    Y = centroids(1:10:end, 2);
    Z = centroids(1:10:end, 3);

    % scatter3(X, Y, Z, 'filled', 'MarkerFaceAlpha', trans);
    plot3(X, Y, Z, '.', 'MarkerSize', 1);
    hold on;
end

axis equal;  % make the axes scale equal
% grid on;
end
