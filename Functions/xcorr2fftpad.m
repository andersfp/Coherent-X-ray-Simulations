function F = xcorr2fftpad(A,B,f)
% Calculate the 2D cross-correlation by using FFTs. It pads the inputs to
% double their size. The output is cropped if f (full) is 0.
%
% Example of usage:
% F = xcorr2fftpad(A,B)
%

% Get the size of the inputs
s1 = size(A,1);
s2 = size(A,2);

% Pad the inputs
A = padarray(A,[s1/2 s2/2],0,'both');
B = padarray(B,[s1/2 s2/2],0,'both');

% Fourier transform and complex conjugate the first input
a = fft2(ifftshift(A));
a = conj(a);

% Fourier transform the second input
b = fft2(ifftshift(B));

% Calculate the cross-correlation
F = fftshift(ifft2(a.*b));

% Only return the center
if ~f
    F = F(s1/2+1:end-s1/2,s2/2+1:end-s2/2);
end
