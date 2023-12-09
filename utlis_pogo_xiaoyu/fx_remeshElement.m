function newElementIndices = fx_remeshElement(indices, numX, numY, numZ)
    % Convert the given element indices to X, Y, and Z coordinates
    [x, y, z] = linearIndicesToCoordinates(indices, numX, numY, numZ);

    % Initialize an array for the new element indices
    newElementIndices = zeros(1, 8*numel(indices));

    % Loop to generate the new element indices
    for i = 1:numel(indices)
        % Calculate the X, Y, and Z coordinates for the sub-element
        subX = (x(i) - 1) * 2 + [1 1 1 1 2 2 2 2];
        subY = (y(i) - 1) * 2 + [1 1 2 2 1 1 2 2];
        subZ = (z(i) - 1) * 2 + [1 2 1 2 1 2 1 2];

        % Convert the X, Y, and Z coordinates to a linear index
        newIndex = coordinatesToLinearIndex(subX, subY, subZ, numX*2, numY*2, numZ*2);

        % Store the new element index in the new array
        newElementIndices(i*8-7: i*8) = newIndex;
    end
end
