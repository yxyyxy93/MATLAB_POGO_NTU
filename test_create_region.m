n = 10; % Define the diamond size; from -n to n
spacing = 1; % Define the spacing between the points

% Initialize an empty matrix to store points
points = [];

% Loop over the range
for x = -n:spacing:n
    for y = -n:spacing:n
        if abs(x) + abs(y) <= n
            points = [points; x, y];
        end
    end
end

% Display the diamond region
scatter(points(:, 1), points(:, 2), 'filled');
axis equal;
