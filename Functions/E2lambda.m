function lambda = E2lambda(E)
% Calculate the wavelength (m) from the photon energy (eV).
%
% Example of usage:
% lambda = E2lambda(E);
%

% Set the conversion factor
c = 12398.42.*1e-10;

% Calculate the wavelength
lambda = c./E;

