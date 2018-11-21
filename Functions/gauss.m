function y = gauss(x,FWHM)
% Calculate an area normalized Gaussian function centered around 0 with
% full width at half max of FWHM.
%
% Example of usage:
% y = gauss(x,FWHM)

% Calculate the c-parameter from the FWHM
c = FWHM./(2.*sqrt(2.*log(2)));

% Calculate the gaussian
y = 1./(c.*sqrt(2.*pi)).*exp(-x.^2./(2.*c.^2));


