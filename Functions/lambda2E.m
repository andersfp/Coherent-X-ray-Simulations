function E = lambda2E(lambda)
% Calculate the photon energy (eV) from the wavelength (m).
%
% Example of usage:
% E = lambda2E(lambda);
%

% Set the conversion factor
c = 12398.42.*1e-10;

% Calculate the wavelength
E = c./lambda;

