function u2 = propTFopt(u1,L,lambda,z)
% Propagation - transfer function approach
% Assumes same x and y side lengths and uniform sampling.
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field

% Get size of input array
M = size(u1,1);

% Get sample interval
dx = L/M;

% Generate frequency coordinates
fx = (-1/(2*dx)):(1/L):(1/(2*dx)-1/L);

% Generate the transfer function
H = exp(-1i*pi*lambda*z*(fx.^2 + fx.'.^2));

% Shift the transfer function
H = fftshift(H);

% Shift and transform the source field
u1 = fft2(fftshift(u1));

% Multiply transfer and transformed source field
u2 = u1.*H;

% Inverse transform to get obervation field
u2 = ifftshift(ifft2(u2));

end
