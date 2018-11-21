function A = imshift2(A,shft)
% Use the Fourier shifting theorem to shift the 2D image A by any amount
% (also subpixel) in the x- and y-direction.
%
% Usage:
% B = imshift2(A,[dx dy]);
%

% Get the number of pixels
nx = size(A,2);
ny = size(A,1);

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

% Shift the axes
x = ifftshift(x);
y = ifftshift(y);

% Separate the shifts
dx = -shft(1);
dy = -shft(2);

% Perform the shifting
A = ifft2(fft2(A).*exp((1i.*2.*pi.*dx./nx).*x + (1i.*2.*pi.*dy./ny).*y));


