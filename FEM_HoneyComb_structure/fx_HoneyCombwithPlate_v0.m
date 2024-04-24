function gm_voided = fx_HoneyCombwithPlate_v0(zTopOuter, centers, rMedian, wallThickness)
% Initialize array to hold all points
% Calculate outer and inner radius
rInner = rMedian - wallThickness / 2;
rOuter = rMedian + wallThickness / 2;
t = (pi/6:pi/3:2*pi)' + pi/6;
% Inner prism vertices (scaled)
xInner = rInner*cos(t);
yInner = rInner*sin(t);
% Outer prism vertices (scaled)
xmed = rMedian *cos(t);
ymed = rMedian *sin(t);
centers = centers - mean(centers);

% % Loop through each center to generate base
% for idx = 1:size(centers, 1)
%     % Extract center coordinates for the current hexagon
%     cx = centers(idx, 1);
%     cy = centers(idx, 2);
%     % Use the polyshape function
%     if (idx==1)
%         base_surf = polyshape(xmed+cx, ymed+cy);
%     else
%         base_surf = addboundary(base_surf, xmed+cx, ymed+cy);
%     end
% end

% Calculate the min and max of centers directly
xMin = min(centers(:, 1)) - rMedian;
xMax = max(centers(:, 1)) + rMedian;
yMin = min(centers(:, 2)) - rMedian;
yMax = max(centers(:, 2)) + rMedian;
base_surf = polyshape([xMin xMax xMax xMin], [yMin yMin yMax yMax]);

tr = triangulation(base_surf);
gm_2d = fegeometry(tr);
gm_3d_base = extrude(gm_2d, 1.2*zTopOuter);
% figure;
% pdegplot(gm_3d_base,"VertexLabels","on","FaceAlpha",0.5);
gm_3d_base = translate(gm_3d_base, [0 0 -0.1*zTopOuter]);

for idx = 1:size(centers, 1)
    % Extract center coordinates for the current hexagon
    cx = centers(idx, 1);
    cy = centers(idx, 2);
    % Use the polyshape function
    if (idx==1)
        bottom_surf = polyshape(xInner+cx, yInner+cy);
    else
        bottom_surf = addboundary(bottom_surf, xInner+cx, yInner+cy);
    end
end
tr = triangulation(bottom_surf);
gm_2d = fegeometry(tr);
gm_3d = extrude(gm_2d, zTopOuter);
% figure;
% pdegplot(gm_3d,"VertexLabels","on","FaceAlpha",0.5)
gm_voided= addVoid(gm_3d_base, gm_3d);
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

function allPoints = interpolateLayers(interpBottomOuter, interpTopOuter, nLayers)
    % Interpolates layers between two sets of points and generates multiple
    % layers in between.
    %
    % Parameters:
    % interpBottomOuter - Points for the bottom surface, as an Nx3 matrix.
    % interpTopOuter - Points for the top surface, as an Nx3 matrix.
    % nLayers - Number of layers to interpolate between the bottom and top.
    %
    % Returns:
    % allPoints - A combined set of points including original and interpolated layers.

    % Pre-allocate array to hold all interpolated points, including top and bottom
    allPoints = zeros(size(interpBottomOuter, 1) * (nLayers + 2), 3);
    % Include the bottom layer
    allPoints(1:size(interpBottomOuter, 1), :) = interpBottomOuter;
    % Loop to interpolate and fill in each layer
    for i = 1:nLayers
        fraction = i / (nLayers + 1); % Fractional distance between bottom and top
        % Linear interpolation for current layer
        interpLayer = interpBottomOuter + fraction * (interpTopOuter - interpBottomOuter);
        % Insert the interpolated layer into the allPoints array
        allPoints(i * size(interpBottomOuter, 1) + (1:size(interpBottomOuter, 1)), :) = interpLayer;
    end
    % Include the top layer
    allPoints(end-size(interpTopOuter, 1)+1:end, :) = interpTopOuter;
end
