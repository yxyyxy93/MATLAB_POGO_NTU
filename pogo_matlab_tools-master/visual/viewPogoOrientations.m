function  [ ] = viewPogoOrientations( m, npix, plane, offset )
% viewPogoOrientations - view material data from Pogo in 2D or 3D
%   viewPogoOrientations( m, npix, plane, offset )
%
%NB - in 3D particular, the speed is often low due to the interpolation
%routine
%
%m - model structure, as loaded from loadPogoInp()
%all other inputs optional (set to [] if unused, or just omit if at end of list)
%matNum - which of the material parameters to plot (default 1 - see manual for more details)
%npix - number of pixels along longest axis (default 100)
%plane - (3D only) which plane to plot - x-y = 1, x-z = 2, y-z = 3, 
%        auto = 0 (default). Auto plots the two dimensions with the largest range.
%offset - (3D only) position along the other axis (default = 0)
%
%Written Jul 2020, Peter Huthwaite

nDims = m.nDims;

%if nargin < 2 || isempty(matNum)
    matNum = 1;
%end
if nargin < 2 || isempty(npix)
    npix = 100;
end
if npix < 1 || npix ~= round(npix)
    error('npix must be an integer and greater than 0')
end


if nargin < 3 || isempty(plane)
    plane = 0;
    %1 is x-y
    %2 is x-z
    %3 is y-z
    %0 is auto (biggest dims)
end
if nDims == 2
    plane = 0;
end
if plane < 0 || plane ~= round(plane) 
    error('plotDir not in range.')
end

if nargin < 4 || isempty(offset)
    offset = 0;
    %position of plane along other axis
end

px = m.nodePos(1,:);
py = m.nodePos(2,:);
if nDims == 3
    pz = m.nodePos(3,:);
end

xMax = max(px);
xMin = min(px);
xRange = xMax - xMin;

yMax = max(py);
yMin = min(py);
yRange = yMax - yMin;

if nDims == 3
    zMax = max(pz);
    zMin = min(pz);
    zRange = zMax - zMin;


    if plane == 0
        if zRange < xRange && zRange < yRange
           plane = 1; 
        elseif yRange < xRange && yRange < zRange
           plane = 2; 
        elseif xRange < yRange && xRange < zRange
           plane = 3; 
        end
    end
else 
    zRange = 0;
end

maxRange = max([xRange yRange zRange]);

dpx = maxRange/npix; %so npix pixels in maximum direction

npx = round((xMax-xMin)/dpx);
npy = round((yMax-yMin)/dpx);
cpx = (xMax+xMin)/2;
cpy = (yMax+yMin)/2;
gpx = ((1:npx)-(npx+1)/2)*dpx+cpx;
gpy = ((1:npy)-(npy+1)/2)*dpx+cpy;

if nDims == 3
    npz = round((zMax-zMin)/dpx);
    cpz = (zMax+zMin)/2;
    gpz = ((1:npz)-(npz+1)/2)*dpx+cpz;
end


if plane == 1 || nDims == 2
    %x-y plane
    [GPX, GPY] = ndgrid(gpx,gpy);
elseif plane == 2
    %x-z plane
    [GPX, GPZ] = ndgrid(gpx,gpz);
else
    %y-z plane
    [GPY, GPZ] = ndgrid(gpy,gpz);
end

nEls = length(m.matTypeRefs);
v = zeros(nEls,1);
for oCnt = 1:length(m.or)
%     if matNum > length(m.matTypes{mCnt}.paramValues)
%         v(m.matTypeRefs == mCnt) = 0;
%     else
%         v(m.matTypeRefs == mCnt) = m.matTypes{mCnt}.paramValues(matNum);
%     end
    v(m.orientRefs == oCnt) = m.or{oCnt}.paramValues(8);
end

%v = m.matTypes{m.matTypeRefs}.paramValues(matNum);

if nDims == 2
    [ex, ey] = getElCents(m);
else
    [ex, ey, ez] = getElCents(m);
end
% figure
% plot3(ex,ey,ez,'b.')

if nDims == 3 && (max(ez)-min(ez) == 0)
    nDims = 2;
end
if (max(ey)-min(ey) == 0)
    nDims = 2;
    ey = ez;
end
if (max(ex)-min(ex) == 0)
    nDims = 2;
    ex = ey;
    ey = ez;
end
    
    %do the interpolation onto a uniform grid:
    if nDims == 2
        %vq = griddata(px,py,v,GPX,GPY); 
        vq = griddata(ex,ey,v,GPX,GPY); 
    else
        if plane == 1
            %x-y plane
            vq = griddata(ex,ey,ez,v,GPX,GPY,zeros(npx,npy)+offset); 
        elseif plane == 2
            %x-z plane
            vq = griddata(ex,ey,ez,v,GPX,zeros(npx,npz),GPZ+offset); 
        else
            %y-z plane
            vq = griddata(ex,ey,ez,v,zeros(npy,npz),GPY,GPZ+offset); 
        end
    end
    
    if plane == 1 || nDims == 2
        %x-y plane
        %cartimagesc(gpx,gpy,vq) 
        imagesc(gpx,gpy,vq.') 
        xlabel('x')
        ylabel('y')
    elseif plane == 2
        %x-z plane
        %cartimagesc(gpx,gpz,vq) 
        imagesc(gpx,gpz,vq.') 
        xlabel('x')
        ylabel('z')
    else
        %y-z plane
        %cartimagesc(gpy,gpz,vq) 
        imagesc(gpy,gpz,vq.') 
        xlabel('y')
        ylabel('z')
    end
    axis xy
    axis image
    colorbar
    title('Orientation')

    drawnow %draw the image now rather than waiting
end
