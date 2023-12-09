function fx_plot_slices(nodes, elements, elsets, x_slice, fignum)

centroids = zeros(size(elsets, 2), 3);  % Preallocate centroid array

for i = 1:10:size(elsets, 2) % skip by 10 to save time
    element_nodes   = elements(elsets(:, i), 1:end); % the first colume includeds
    nodes_pos       = nodes(element_nodes, 1:end);
    centroids(i, :) = mean(nodes_pos, 1);
end

% skip to save memory
X = centroids(1:1:end, 1);
Y = centroids(1:1:end, 2);
Z = centroids(1:1:end, 3);

% Find all points where Z is close to z_slice within a tolerance
tolerance = 5e-5;
indices = abs(X - x_slice) < tolerance;

% Plot the 2D slice
figure(fignum);
scatter(Y(indices), Z(indices), '.');
title(['2D slice at x = ' num2str(x_slice)]);
xlabel('Y');
ylabel('Z');
hold on;

end