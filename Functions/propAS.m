function [u2,x2,y2] = propAS(u1,lambda,d1,d2,z)
% Propagate - using the angular spectrum method

% Get the size of the input
N = size(u1,1); % assume square grid

% Get the wavenumber
k = 2*pi/lambda; % optical wavevector

% Generate source plane coordinates
[x1,y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
r1sq = x1.^2 + y1.^2;

% Generate spatial frequency (of source plane) coordinates
df1 = 1/(N*d1);
[fX,fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
fsq = fX.^2 + fY.^2;

% Calculate the scaling parameter
m = d2/d1;

% Generate observation plane coordinates
[x2,y2] = meshgrid((-N/2 : 1 : N/2 - 1) * d2);
r2sq = x2.^2 + y2.^2;

% Calculate quadratic phase factors
Q1 = exp(1i*k/2*(1-m)/z*r1sq);
Q2 = exp(-1i*pi^2*2*z/m/k*fsq);
Q3 = exp(1i*k/2*(m-1)/(m*z)*r2sq);

% Compute the propagated field
%u2 = Q3.* ift2(Q2 .* ft2(Q1 .* u1 / m, d1), df1);
%u2 = Q3.*ifftshift(ifft2(fftshift(Q2).*fft2(fftshift(Q1.*u1/m))*d1^2)*(N*df1)^2);
u2 = Q3.*ifftshift(ifft2(fftshift(Q2).*fft2(fftshift(Q1.*u1/m))));
