function [u2,L2] = propFF(u1,L1,lambda,z)
% Propagation - Fraunhofer pattern
% Assumes uniform sampling.
% u1 - source plane field
% L1 - source plane side length
% lambda - wavelength
% z - propagation distance
% L2 - observation plane side length
% u2 - observation plane field

% Get size of input array
[M,N] = size(u1);

% Get sample interval
dx1 = L1/M;

% Get wavenumber
k = 2*pi/lambda;

% Get observation plane side length
L2 = lambda*z/dx1;

% Get the observation plane sampling interval
dx2 = lambda*z/L1;

% Generate obervation plane coordinates
x2 = (-L2/2):dx2:(L2/2 - dx2);

% Generate 2D observation plane coordinates
[X2,Y2] = meshgrid(x2,x2);

% Generate chirp
c = 1./(1i.*lambda.*z).*exp(1i.*k./(2*z).*(X2.^2 + Y2.^2));

% Generate observation plane field
u2 = c.*ifftshift(fft2(fftshift(u1))).*dx1.^2;

