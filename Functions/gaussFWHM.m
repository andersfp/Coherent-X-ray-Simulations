function y = gaussFWHM(x,FWHM,norm)
% Calculate a Gaussian function centered around 0 with full width at half 
% max of FWHM. If norm is 1, then it will be area normalized.
%
% Example of usage:
% y = gauss(x,FWHM,norm)

% Calculate the c-parameter from the FWHM
c = FWHM./(2.*sqrt(2.*log(2)));

% Calculate the gaussian
if norm
    y = 1./(c.*sqrt(2.*pi)).*exp(-x.^2./(2.*c.^2));
else
    y = exp(-x.^2./(2.*c.^2));
end

