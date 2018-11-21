function u2 = propTF(u1,L,lambda,z)
% Propagation - transfer function approach
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

% Generate frequency coordinates
fx = (-1/(2*dx)):(1/L):(1/(2*dx)-1/L);

% Make 2D frequency grid
[Fx,Fy] = meshgrid(fx,fx);

% Generate the transfer function
H = exp(-1i*pi*lambda*z*(Fx.^2 + Fy.^2));

% Shift the transfer function
H = fftshift(H);

% Shift and transform the source field
U1 = fft2(fftshift(u1));

% Multiply transfer and transformed source field
U2 = U1.*H;

% Inverse transform to get obervation field
u2 = ifftshift(ifft2(U2));

end
