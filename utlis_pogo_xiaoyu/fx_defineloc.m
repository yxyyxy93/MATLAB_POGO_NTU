function B = fx_defineloc(dx, dy, nx, ny)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[X, Y] = meshgrid(-ny/2:ny/2, -nx/2:nx/2);
centers = nan(nx+1, ny+1, 3);

centers(:,:,1) = X*dx;
centers(:,:,2) = Y*dy;
centers(:,:,3) = 0;

% Reshape A to a 2D array B of size mn x p
B = reshape(centers, (nx+1)*(ny+1), 3);

% Save B to a text file
dlmwrite('centers_array.txt', B);

end
