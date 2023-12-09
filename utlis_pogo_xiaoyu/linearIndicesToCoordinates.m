% Define the conversion functions for arrays
function [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, ~)
    indices = double(indices);
    z       = ceil(indices / (numX * numY));
    indices = indices - (z - 1) * (numX * numY);
    y       = ceil(indices / numX);
    x       = indices - (y - 1) * numX;
end
