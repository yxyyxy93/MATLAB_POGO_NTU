function gm_voided = fx_HoneyCombwithPlate(zTopOuter, rows, cols, rMedian, wallThickness)
    % Initialize array to hold all points
    % Calculate the outer and inner radii
    rInner = rMedian - wallThickness / 2;
    t = (pi/6:pi/3:2*pi)' + pi/6;
    % Inner prism vertices (scaled)
    xInner = rInner*cos(t);
    yInner = rInner*sin(t);

    % Updated for an expanded grid (adding one extra row and column on each side)
    dx = 1.5 * rMedian; % Horizontal distance between centers
    dy = sqrt(3) * rMedian / 2; % Vertical distance for the shift
    centers = zeros(rows*cols, 2);  % Resize the centers array

    % Populate the centers array
    index = 1;
    for row = -1:rows  % Start from -1 to add an extra row at the beginning
        for col = -1:cols  % Start from -1 to add an extra column at the beginning
            cx_temp =  col * dx;
            cy_temp =  - 2 * row * dy - mod(col + 1, 2) * dy;  % Adjust for extra row and staggered columns
            centers(index, :) = [cx_temp, cy_temp];
            index = index + 1;
        end
    end
    % Calculate the min and max of centers directly
    xMin = min(centers(:, 1)) - wallThickness / 2;
    xMax = max(centers(:, 1)) + wallThickness / 2;
    yMin = min(centers(:, 2)) - wallThickness / 2;
    yMax = max(centers(:, 2)) + wallThickness / 2;
    base_surf = polyshape([xMin xMax xMax xMin], [yMin yMin yMax yMax]);

    % Calculate the min and max of centers directly
    xMin = min(centers(:, 1)) - wallThickness;
    xMax = max(centers(:, 1)) + wallThickness;
    yMin = min(centers(:, 2)) - wallThickness;
    yMax = max(centers(:, 2)) + wallThickness;
    base_surf_outer = polyshape([xMin xMax xMax xMin], [yMin yMin yMax yMax]);

    tr = triangulation(base_surf_outer);
    model = createpde;
    tnodes = tr.Points';
    telements = tr.ConnectivityList';
    geometryFromMesh(model,tnodes,telements);
    gm_2d = fegeometry(model.Geometry);
    gm_3d_base = extrude(gm_2d, 1.2*zTopOuter);
    gm_3d_base = translate(gm_3d_base, [0 0 -0.1*zTopOuter]);

    % Create bottom surface hexagons, ensuring they are all within the boundary
    bottom_surf = polyshape();
    for idx = 1:size(centers, 1)
        cx = centers(idx, 1);
        cy = centers(idx, 2);
        hex_shape = polyshape(xInner + cx, yInner + cy);
        hex_shape = intersect(hex_shape, base_surf);
        if idx == 1
            bottom_surf = hex_shape;
        else
            bottom_surf = union(bottom_surf, hex_shape);
        end
    end

    % % Create a new figure
    % figure;
    % hold on;
    % axis equal;
    % % Plot the base surface in blue
    % plot(base_surf, 'FaceColor', 'blue', 'FaceAlpha', 0.5);
    % title('Base and Hexagonal Surfaces');
    % xlabel('X coordinate');
    % ylabel('Y coordinate');
    % % Plot the bottom surface (hexagons) in red
    % plot(bottom_surf, 'FaceColor', 'red', 'FaceAlpha', 0.5);
    % % Add grid and legend
    % grid on;
    % legend('Base Surface', 'Hexagonal Surface', 'Location', 'best');
    % hold off;

    tr = triangulation(bottom_surf);
    model = createpde;
    tnodes = tr.Points';
    telements = tr.ConnectivityList';
    geometryFromMesh(model,tnodes,telements);
    gm_2d = fegeometry(model.Geometry);
    gm_3d = extrude(gm_2d, zTopOuter);

    % Attempt to add voids
    try
        gm_voided = addVoid(gm_3d_base, gm_3d);
    catch ME
        disp('Error encountered while adding voids:');
        disp(ME.message);
        % Optional: Return the base geometry for visualization
        gm_voided = gm_3d_base;
    end
end
