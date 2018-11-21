% Initialization
clear;
close all;
clc;


%% Set up the experimental parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Get material properties
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6;
T = 1.6e-3;
Tweb = 2e-6;
N = 70;

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set the magnification
M = 10;
d2 = f.*phi.*(M + cos(N.*phi))./sin(N.*phi);
d1 = fN.*(d2 + f.*phi.*tan(N.*phi))./(d2 - fN);

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);


%% Generate the object
% Set the number of pixels
m = 256;

% Field-of-view
dx = 10e-6;

% Generate axis
x0 = ((-m/2):(m/2 - 1)).'./m.*dx;

% Make object
w = 8e-6;
E0 = recta(x0/w).*recta(x0.'/w);
%E0 = E0.*exp(10.*1i.*2.*pi.*abs(x0)./dx);


%% Partial coherence parameters
% Coherence length: lambda*L/s, ID06: L = 18 m, s = 1 mm
L = 18;
s = 1e-3;
l = lambda*L/s;

% Set the energy bandwidth
dE = 1e-4*E;


%% Save the settings
% Save parameters
save('Example_Image_Settings.mat');



