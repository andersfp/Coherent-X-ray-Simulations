function F = xcorr2fft(A,B)
% Calculate the 2D cross-correlation by using FFTs.
%
% Example of usage:
% F = xcorr2fft(A,B)
%

% Fourier transform and complex conjugate the first input
a = fft2(ifftshift(A));
a = conj(a);

% Fourier transform the second input
b = fft2(ifftshift(B));

% Calculate the cross-correlation
F = fftshift(ifft2(a.*b));

