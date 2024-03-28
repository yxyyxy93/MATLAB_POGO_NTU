% Function to plot a single hexagon
function fx_plotHexagon(cx, cy, r)
    % Angles for the hexagon vertices (60 degree steps)
    theta = (0:6)*2*pi/6; % 0:6 gives us the 7 points needed to return to the start
    % Hexagon vertices
    x = r * cos(theta) + cx;
    y = r * sin(theta) + cy;
    plot(x, y, 'LineWidth', 2);
end