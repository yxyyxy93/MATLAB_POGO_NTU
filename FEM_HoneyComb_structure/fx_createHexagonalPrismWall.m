function gm = fx_createHexagonalPrismWall(zBottomOuter, zTopOuter, cx, cy, rMedian, wallThickness, nPoints_inner, nPoints_outer)
    % Calculate outer and inner radius
    rOuter = rMedian + wallThickness / 2;
    rInner = rMedian - wallThickness / 2;

    t = (pi/6:pi/3:2*pi)' + pi/6;

    % Outer prism vertices
    xOuter = rOuter*cos(t) + cx;
    yOuter = rOuter*sin(t) + cy;

    % Inner prism vertices (scaled)
    xInner = rInner*cos(t) + cx;
    yInner = rInner*sin(t) + cy;

    % Median vertices
    xMedian = rMedian * cos(t) + cx;
    yMedian = rMedian * sin(t) + cy;

    % Interpolate outer and inner prisms
    interpBottomOuter = interpolateHexagon(xOuter, yOuter, zBottomOuter, nPoints_outer);
    interpTopOuter = interpolateHexagon(xOuter, yOuter, zTopOuter, nPoints_outer);
    interpBottomInner = interpolateHexagon(xInner, yInner, zBottomOuter, nPoints_inner);
    interpTopInner = interpolateHexagon(xInner, yInner, zTopOuter, nPoints_inner);

    % Interpolate median vertices
    interpBottomMedian = interpolateHexagon(xMedian, yMedian, zBottomOuter, (nPoints_inner+nPoints_outer)/2);
    interpTopMedian = interpolateHexagon(xMedian, yMedian, zTopOuter, (nPoints_inner+nPoints_outer)/2);

    % Combine all points
    fullPrismPoints = [interpBottomOuter; interpTopOuter; interpBottomInner; interpTopInner; interpBottomMedian; interpTopMedian];

    % Create alphaShape and generate mesh
    shp = alphaShape(fullPrismPoints, 'HoleThreshold', max(fullPrismPoints(:))*2);
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

