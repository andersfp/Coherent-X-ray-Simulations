function A = imshift3(A,shft)
% Use the Fourier shifting theorem to shift the 3D image stack A by any 
% amount (also subpixel) in the x-, y-, and z-direction.
%
% Usage:
% B = imshift2(A,[dx dy dz]);
%

% Get the number of pixels
nx = size(A,2);
ny = size(A,1);
nz = size(A,3);

% Generate the axes
if mod(nx,2) == 0
    x = (-nx/2):(nx/2 - 1);
else
    x = (-(nx - 1)/2):((nx - 1)/2);
end
if mod(ny,2) == 0
    y = (-ny/2):(ny/2 - 1);
else
    y = (-(ny - 1)/2):((ny - 1)/2);
end
y = y.';
if mod(nz,2) == 0
    z = (-nz/2):(nz/2 - 1);
else
    z = (-(nz - 1)/2):((nz - 1)/2);
end
z = permute(z,[3 1 2]);

% Shift the axes
x = ifftshift(x);
y = ifftshift(y);
z = ifftshift(z);

% Separate the shifts
dx = -shft(1);
dy = -shft(2);
dz = -shft(3);

% Perform the shifting
A = ifftn(fftn(A).*exp((1i.*2.*pi.*dx./nx).*x + (1i.*2.*pi.*dy./ny).*y + (1i.*2.*pi.*dz./nz).*z));

