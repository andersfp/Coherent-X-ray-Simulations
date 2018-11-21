function C = fconv1(A,B)
% Perform a 1D convolution between A and B by using Fourier transforms

% Make the convolution
C = fftshift(ifft(fft(ifftshift(A)).*fft(ifftshift(B))));

