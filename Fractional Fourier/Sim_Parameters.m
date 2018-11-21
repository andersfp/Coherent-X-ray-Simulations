% Initialization
clear;
close all;
clc;


%% Set parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42/E;

% Set CRL parameters
R_CRL = 50e-6;
delta = 1.17843774e-6;
f = R_CRL/(2*delta);
T = 1.6e-3;
N = 69;
mu = 1/24590.6e-6;
Tweb = 2e-6;

% Calculate focal length
phi = sqrt(T/f);
fN = f*phi*cot(N*phi);

% Set object and image positions (see CRL_Focal_Length_Sim)
M = 10;
d0 = fN*(1 + 1/(M*cos(N*phi)));
dd = fN*(1 + M/cos(N*phi));

% Set the effective aperture size at the end of the CRL (Gaussian RMS)
aper = 8.246e-5;

% Calculate sigma_a
sigA = sqrt(R_CRL / (mu * N * (d0^2 + (f * phi)^2)))/sqrt(1 + 1 / N - 1 / (N * phi) * sin((N + 1) * phi) * cos((N - 1) * phi + 2 * atan(d0 / (f * phi))));

% Large sim resul saving location
p = 'C:\Users\anfils\Documents\Simulation_Results\CRL_Paper_2\';

% Save the parameters
save('Sim_Parameters.mat');
save([p 'Sim_Parameters.mat']);



