
% Example usage:
center  = [5, 5];
n       = 14;
spacing = 2.5;
points  = generateDiamond(center, n, spacing);

x = points(:, 1);
y = points(:, 2);

% Display the diamond region
scatter(x, y, 'filled');
axis equal;

% Label each point with its sequence number
for i = 1:length(x)
    text(x(i), y(i), num2str(i), 'HorizontalAlignment', 'center');
end
xlabel('X');
ylabel('Y');
grid on;

rectangular_points = fx_diamondToRectangle(points);

% Display the diamond region
figure,
x = rectangular_points(:, 1);
y = rectangular_points(:, 2);
scatter(x, y, 'filled');
axis equal;

% Label each point with its sequence number
for i = 1:length(x)
    text(x(i), y(i), num2str(i), 'HorizontalAlignment', 'center');
end
xlabel('X');
ylabel('Y');
grid on;

function points = generateDiamond(center, n, spacing)
    % center: [x_center, y_center]
    % n: half of the diamond's diagonal length
    % spacing: diagonal spacing between the points

    % Calculate the horizontal and vertical increments for diagonal spacing
    delta = spacing / sqrt(2);

    % Generate points for the top half of the diamond, starting from the top point
    top_points = [];
    
    % Iterate over the main diagonal
    for d = 0:1:n
        y = center(2) + (n - d)*delta;
        x_start = center(1) - d*delta;

        % For each y, create a row of points
        row_points = x_start:2*delta:center(1) + d*delta;
        for x = row_points
            top_points = [top_points; x, y];
        end
    end
    
    % Reflect the top half to get the bottom half
    bottom_points = [top_points(:,1), 2*center(2) - top_points(:,2)];
    % Exclude the points on the primary diagonal to avoid duplicates
    bottom_points = bottom_points(bottom_points(:,2) ~= center(2), :);
    % Combine the top and bottom points
    points = [top_points; bottom_points];
end

function rectangular_points = fx_diamondToRectangle(points)
    % Find the bottommost point of the diamond
    [~, idx] = min(points(:,2));
    x_bottom = points(idx,1);
    y_min = points(idx,2);
    
    % Constants for the left diagonal
    c1 = y_min - x_bottom;
    
    % Calculate the distance from each point to the left diagonal
    distance_to_diag1 = abs(points(:,1) - points(:,2) + c1) / sqrt(2);
    
    % Sort primarily by distance to the left diagonal and then by y-coordinate for ties
    [~, sorted_indices] = sortrows(distance_to_diag1);
    
    num_row = sqrt(size(points, 1));
    num_col = size(points, 1) / num_row;
    for i = 0:num_row-1
        [~, sorted_indices_row] = sortrows(points(sorted_indices(i*num_col+1:(i+1)*num_col), 1));
        sorted_indices(i*num_col+1:(i+1)*num_col) = sorted_indices(i*num_col+sorted_indices_row);
    end

    rectangular_points = points(sorted_indices, :);
end






