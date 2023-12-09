function points = fx_generateDiamond_region(center, n, spacing)
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

