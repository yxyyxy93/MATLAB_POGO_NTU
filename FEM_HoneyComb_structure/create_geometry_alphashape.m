% Outer prism dimensions
zBottomOuter = 0;
zTopOuter = 1;
cx = 0.5; % Center X coordinate
cy = 1; % Center Y coordinate
r = 2; % radius

tOuter = (pi/6:pi/3:2*pi)' + pi/6;
xOuter = r*cos(tOuter) + cx;
yOuter = r*sin(tOuter) + cy;

% Scale factor for the inner prism relative to the outer prism
scaleFactor = 0.8; % Example: inner prism is half the size of the outer

% Inner prism dimensions (scaled)
xInner = scaleFactor * r*cos(tOuter) + cx;
yInner = scaleFactor * r*sin(tOuter) + cy;

xmedian = (0.5+scaleFactor/2) * r*cos(tOuter) + cx;
ymedian = (0.5+scaleFactor/2) * r*sin(tOuter) + cy;

% % Combine into a single set of points (including both prisms)
% vertices = [xBottomOuter, yBottomOuter, zBottomOuter*ones(length(xBottomOuter),1); 
%             xTopOuter, yTopOuter, zTopOuter*ones(length(xTopOuter),1);
%             xBottomInner, yBottomInner, zBottomOuter*ones(length(xBottomInner),1);
%             xTopInner, yTopInner, zTopOuter*ones(length(xTopInner),1)];
% 
% % Scatter plot for bottom and top frames
% scatter3(vertices(:, 1), vertices(:, 2), vertices(:, 3), 'LineWidth', 2);
% % Annotating each point with its index number
% numVertices = size(vertices, 1); % Total number of vertices
% for idx = 1:numVertices
%     % Adjust the position of the text for better visibility
%     textOffset = 0.1; % Offset for text annotation to avoid overlapping with the point
%     text(vertices(idx, 1) + textOffset, vertices(idx, 2) + textOffset, vertices(idx, 3),...
%         num2str(idx), 'FontSize', 8, 'HorizontalAlignment', 'center');
% end

% **************** interpolation ******************
nPoints_inner = 8; % Adjust as needed for density
nPoints_outer = 10; % Adjust as needed for density
% nPoints_v_edge = 1; % Adjust as needed for density
% nPoints_h_edge = 1; % Adjust as needed for density
% Interpolate the outer prism
interpBottomOuter = interpolateHexagon(xOuter, yOuter, zBottomOuter, nPoints_outer);
interpTopOuter = interpolateHexagon(xOuter, yOuter, zTopOuter, nPoints_outer);
% wallPointsOuter = interpolatePrismWalls(interpBottomOuter(:,1), interpBottomOuter(:,2), zBottomOuter, ...
%     interpTopOuter(:,1), interpTopOuter(:,2), zTopOuter, nPoints_v_edge, nPoints_h_edge);
% Interpolate the inner prism
interpBottomInner = interpolateHexagon(xInner, yInner, zBottomOuter, nPoints_inner);
interpTopInner = interpolateHexagon(xInner, yInner, zTopOuter, nPoints_inner);
% wallPointsInner = interpolatePrismWalls(interpBottomInner(:,1), interpBottomInner(:,2), zBottomOuter, ...
%     interpTopInner(:,1), interpTopInner(:,2), zTopOuter, nPoints_v_edge, nPoints_h_edge);
% Combine all points
fullPrismPointsOuter = [interpBottomOuter; interpTopOuter];
fullPrismPointsInner = [interpBottomInner; interpTopInner];
interpbotmed = interpolateHexagon(xmedian, ymedian, zBottomOuter, (nPoints_inner+nPoints_outer)/2);
interptopmed = interpolateHexagon(xmedian, ymedian, zTopOuter, (nPoints_inner+nPoints_outer)/2);
% Create an alphaShape object of the hexagonal prism wall
fullPrismPoints = [fullPrismPointsOuter; fullPrismPointsInner; interpbotmed; interptopmed];
% Visualize the combined outer and inner prisms
figure, hold on
scatter3(fullPrismPoints(:, 1), fullPrismPoints(:, 2), fullPrismPoints(:, 3), 'LineWidth', 2);
hold off
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Interpolated Hexagonal Prism with Wall');
grid on
% *******************

% Create an alphaShape object of the hexagonal prism wall
shp = alphaShape(fullPrismPoints, 'HoleThreshold', max(fullPrismPoints(:))*2); 

% Adjust the 'HoleThreshold' to ensure the inner space is recognized as a hole.
% You may need to tweak this value based on the size of your prisms.

% Generate a surface mesh of the alphaShape object
[elements, nodes] = boundaryFacets(shp);

% Convert to fegeometry format
gm = fegeometry(nodes, elements);

% Generate a mesh for the geometry
gm = generateMesh(gm, 'GeometricOrder', 'linear', 'Hmax', 0.1);

% Plot the geometry with the face labels
figure;
pdegplot(gm, 'FaceLabels', "on", 'FaceAlpha', 0.5);

% Plot the mesh
figure;
pdemesh(gm);

function interpVertices = interpolateHexagon(x, y, z, nPoints)
    % Initialize the output matrix with the correct size for efficiency
    totalInterpPoints = (nPoints + 1) * (length(x) - 1);
    interpVertices = zeros(totalInterpPoints, 3);
    % Precompute the repeated z values to avoid recalculating them
    zInterp = z * ones(1, nPoints + 2);
    % Vectorized loop to populate interpVertices
    for i = 1:length(x)
        % Determine the current and next vertices
        xNext = circshift(x, -1);
        yNext = circshift(y, -1);
        % Compute the interpolated values between vertices
        xInterp = linspace(x(i), xNext(i), nPoints + 2);
        yInterp = linspace(y(i), yNext(i), nPoints + 2);
        % Calculate the index range in the output array for these points
        idxStart = (i-1) * (nPoints + 1) + 1;
        idxEnd = idxStart + nPoints;
        % Fill the corresponding section of the output array
        interpVertices(idxStart:idxEnd, :) = [xInterp(1:end-1)' yInterp(1:end-1)' zInterp(1:end-1)'];
    end
    % Remove potential duplicates at the end due to circshift overlap
    interpVertices = unique(interpVertices, 'rows', 'stable');
end


% % Function to interpolate the walls of a hexagonal prism
% function wallPoints = interpolatePrismWalls(xBottom, yBottom, zBottom, xTop, yTop, zTop, nVerticalPoints, nHorizontalPoints)
%     wallPoints = [];
%     for i = 1:length(xBottom)
%         % Bottom and top vertex of current edge
%         bottomVertex = [xBottom(i), yBottom(i), zBottom];
%         topVertex = [xTop(i), yTop(i), zTop];
%         % Interpolate points vertically along the current edge
%         xLine = linspace(bottomVertex(1), topVertex(1), nVerticalPoints);
%         yLine = linspace(bottomVertex(2), topVertex(2), nVerticalPoints);
%         zLine = linspace(bottomVertex(3), topVertex(3), nVerticalPoints);
%         for j = 1:nVerticalPoints - 1
%             % Interpolate points horizontally to form the wall
%             if i == length(xBottom)
%                 nextBottomVertex = [xBottom(1), yBottom(1), zBottom];
%                 nextTopVertex = [xTop(1), yTop(1), zTop];
%             else
%                 nextBottomVertex = [xBottom(i+1), yBottom(i+1), zBottom];
%                 nextTopVertex = [xTop(i+1), yTop(i+1), zTop];
%             end
%             % Calculate points for bottom and top edge of the current horizontal segment
%             bottomEdgePoint = [xLine(j), yLine(j), zLine(j)];
%             topEdgePoint = [xLine(j+1), yLine(j+1), zLine(j+1)];
%             xEdge = linspace(bottomEdgePoint(1), topEdgePoint(1), nHorizontalPoints);
%             yEdge = linspace(bottomEdgePoint(2), topEdgePoint(2), nHorizontalPoints);
%             zEdge = linspace(bottomEdgePoint(3), topEdgePoint(3), nHorizontalPoints);
%             % Combine the interpolated points
%             wallPoints = [wallPoints; [xEdge' yEdge' zEdge']];
%         end
%     end
% end