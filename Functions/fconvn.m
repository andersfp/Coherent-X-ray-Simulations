function C = fconvn(A,B)
% Perform an N-D convolution between A and B by using Fourier transforms

% Make the convolution
C = fftshift(ifftn(fftn(ifftshift(A)).*fftn(ifftshift(B))));

