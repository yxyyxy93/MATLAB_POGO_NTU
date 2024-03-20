function hex_centers = fx_calculate_honeycomb_centers(rows, side_length, center_offset, centerZ)
    hex_centers = [];
    vertical_spacing = side_length * 1.1; % Spacing between centers of hexagons vertically
    
    for row = 1:rows
        % The x-coordinate is constant because there's only one column
        x = center_offset(1);
        % Calculate the y-coordinate
        y = (row - 1) * vertical_spacing + center_offset(2);
        % Store the center coordinates
        hex_centers(end + 1, :) = [x, y, centerZ];
    end
end

