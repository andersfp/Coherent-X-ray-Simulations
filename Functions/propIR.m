function u2 = propIR(u1,L,lambda,z)
% Propagation - impulse response approach
% Assumes same x and y side lengths and uniform sampling.
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field

% Get size of input array
[M,N] = size(u1);

% Get sample interval
dx = L/M;

% Get wavenumber
k = 2*pi/lambda;

% Generate spatial coordinates
x = (-L/2):(dx):(L/2-dx);

% Make 2D coordinate grid
[X,Y] = meshgrid(x,x);

% Generate the impulse function
h = 1./(1i.*lambda.*z).*exp(1i.*k./(2.*z).*(X.^2 + Y.^2));

% Get the transfer function
H = fft2(fftshift(h)).*dx.^2;

% Shift and transform the source field
U1 = fft2(fftshift(u1));

% Multiply transfer and transformed source field
U2 = U1.*H;

% Inverse transform to get obervation field
u2 = ifftshift(ifft2(U2));

end
