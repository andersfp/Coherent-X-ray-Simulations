function f = MTF_Fitting(r,x)
% Model the Fourier domain signal from a reconstruction data. Convolute an
% ideal rectangle function (true object) with a Gaussian (PSF), and
% multiply by a second rectangle function (support). Finally Fourier
% transform, take the amplitude, and normalize.

% Separate the variables
w0 = x(1);
ws = x(2);
sig = x(3);

% Generate the ideal object
a = recta(r./w0);

% Generate the Gaussian PSF
g = exp(-r.^2./(2.*sig.^2));

% Convolute the two real space function
b = fconv1(a,g);

% Generate the support
s = recta(r./ws);

% Multiply the support and lower resolution object
c = b.*s;

% Perform Fourier transform
f = abs(fftshift(fft(c)));

% Normalize the Fourier transform
f = f./max(f);

