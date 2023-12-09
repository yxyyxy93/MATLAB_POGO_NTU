function [index] = coordinatesToLinearIndex(x, y, z, numX, numY, ~)
    index = (z - 1) * numX * numY + (y - 1) * numX + x;
end