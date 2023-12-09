clear
% define element size
dx = 0.2e-3;
% lambda = 4mm so this is 20 els per wavelength at 1 MHz

% define nodes along the bottom of the free meshed domain
nodex = (0.1: dx :0.2) .';

% put into a points array , along with another point at the top to
% make a triangle
points = [ 
    [ nodex zeros(length(nodex), 1)] 
    0.15 0.06].';

% add points inside for a rectangle to be deleted
rectPoints = [ 
    0.14 0.01
    0.16 0.01
    0.16 0.025
    0.14 0.025].';
% combine :
points = [ points rectPoints ];

nPointsOuter = length ( nodex ) +1;

% define the segments to joint the points around the outside together
segsOut = [ (1: nPointsOuter ).' [2: nPointsOuter 1].'].';

% segments for the internal rectangle
segsRect = [(nPointsOuter+1: nPointsOuter+4).' ...
    [(nPointsOuter+2: nPointsOuter+4) nPointsOuter+1].'].';

% combined segments
segs = [ segsOut segsRect ];

% one hole which will be deleted
holes = [0.148 0.015].';

% save into a poly format
savePoly('temp. poly ', points, segs, holes);

% IMPORTANT
% This might need to be changed depending on your system setup
% Essentially you want to run pogoMesh with the file just generated
% and maximum length set to dx *1.5
system (sprintf('pogoMesh temp . poly -l %f', dx *1.5))

% having done the meshing can delete the file
delete ('temp . poly ')

% load in mesh
m = loadVtuMesh ('temp . vtu ');
% can delete mesh file now loaded
delete ('temp . vtu ')

% generate a 2D grid
% centre locations for grid
% - note that we want to make the top of this line up with the nodes from
% the free mesh
cx = 0.15;
cy = -0.05;
% points in x and y
nx = 0.1/ dx +1;
ny = 0.1/ dx +1;
% actually make grid
m2 = genGrid2D (nx ,ny ,dx ,dx ,cx ,cy);

% get the overlapping node numbers from each model ; the top line from 2, the
% bottom line from 1
n2 = (1: nx) + (ny -1) *nx;
n1 = 1: nx;



















