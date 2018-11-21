function C = fconv2(A,B)
% Perform a 2D convolution between A and B by using Fourier transforms

% Make the convolution
C = fftshift(ifft2(fft2(ifftshift(A)).*fft2(ifftshift(B))));

