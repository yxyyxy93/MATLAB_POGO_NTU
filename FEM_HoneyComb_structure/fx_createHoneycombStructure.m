function gm = fx_createHoneycombStructure(zBottomOuter, zTopOuter, centers, rMedian, wallThickness, nPoints_inner, nPoints_outer)
    % Initialize array to hold all points
    allPrismPoints = [];

    % Calculate outer and inner radius
    rOuter = rMedian + wallThickness / 2;
    rInner = rMedian - wallThickness / 2;
    t = (pi/6:pi/3:2*pi)' + pi/6;
    % Outer prism vertices
    xOuter = rOuter*cos(t);
    yOuter = rOuter*sin(t);
    % Inner prism vertices (scaled)
    xInner = rInner*cos(t);
    yInner = rInner*sin(t);
    % Median vertices
    xMedian = rMedian * cos(t);
    yMedian = rMedian * sin(t);
    % Interpolate vertices for outer and inner prisms, and median vertices
    interpBottomOuter = interpolateHexagon(xOuter, yOuter, zBottomOuter, nPoints_outer);
    interpTopOuter = interpolateHexagon(xOuter, yOuter, zTopOuter, nPoints_outer);
    interpBottomInner = interpolateHexagon(xInner, yInner, zBottomOuter, nPoints_inner);
    interpTopInner = interpolateHexagon(xInner, yInner, zTopOuter, nPoints_inner);
    interpBottomMedian = interpolateHexagon(xMedian, yMedian, zBottomOuter, (nPoints_inner + nPoints_outer) / 2);
    interpTopMedian = interpolateHexagon(xMedian, yMedian, zTopOuter, (nPoints_inner + nPoints_outer) / 2);

    % Loop through each center to generate hexagonal prisms
    for idx = 1:size(centers, 1)
        % Extract center coordinates for the current hexagon
        cx = centers(idx, 1);
        cy = centers(idx, 2);
        % Combine all points
        fullPrismPoints = [interpBottomOuter; interpTopOuter; interpBottomInner; interpTopInner; interpBottomMedian; interpTopMedian];
        fullPrismPoints(:, 1) = fullPrismPoints(:, 1) + cx;
        fullPrismPoints(:, 2) = fullPrismPoints(:, 2) + cy;
        allPrismPoints = [allPrismPoints; fullPrismPoints]; % Accumulate all points
    end
    
    % Create alphaShape and generate mesh for the entire structure
    shp = alphaShape(allPrismPoints, 'HoleThreshold', max(allPrismPoints(:))*2);
    [elements, nodes] = boundaryFacets(shp);
    gm = fegeometry(nodes, elements);
end

function interpVertices = interpolateHexagon(x, y, z, nPoints)
    totalInterpPoints = (nPoints + 1) * (length(x) - 1);
    interpVertices = zeros(totalInterpPoints, 3);
    zInterp = z * ones(1, nPoints + 2);
    for i = 1:length(x)
        xNext = circshift(x, -1);
        yNext = circshift(y, -1);
        xInterp = linspace(x(i), xNext(i), nPoints + 2);
        yInterp = linspace(y(i), yNext(i), nPoints + 2);
        idxStart = (i-1) * (nPoints + 1) + 1;
        idxEnd = idxStart + nPoints;
        interpVertices(idxStart:idxEnd, :) = [xInterp(1:end-1)' yInterp(1:end-1)' zInterp(1:end-1)'];
    end
    interpVertices = unique(interpVertices, 'rows', 'stable');
end
