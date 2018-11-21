function u2 = propAS1gpuss(u1,lambda,d1,z)
% Propagate - using the angular spectrum method. Special case using m = 1.

% Get the size of the input
N = size(u1,1); % assume square grid

% Stop if the input is too large
if N^2 > 110250000
    error('Input array cannot be larger than 10500x10500.');
end

% Get the wavenumber
k = 2*pi/lambda; % optical wavevector

% Generate spatial frequency (of source plane) coordinates
df1 = 1/(N*d1);
[fX,fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
fsq = fX.^2 + fY.^2;
clear fX fY;

% Calculate quadratic phase factors
Q2 = exp(-1i*pi^2*2*z/k*fsq);
clear fsq;

% Compute the propagated field
ug = gpuArray(fftshift(u1));
ug = ifft2(fftshift(Q2).*fft2(ug));
u2 = gather(ug);
u2 = ifftshift(u2);
