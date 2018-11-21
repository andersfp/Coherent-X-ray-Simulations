function y = gaussRMS(x,rms)
% Calculate an unnormalized Gaussian function centered around 0 with
% RMS width of rms.
%
% Example of usage:
% y = gaussRMS(x,rms)

% Calculate the gaussian
y = exp(-x.^2./(2.*rms.^2));


