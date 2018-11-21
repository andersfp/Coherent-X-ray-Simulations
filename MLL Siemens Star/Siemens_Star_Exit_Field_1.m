% Initialization
clear;
close all;
clc;


%% Generate the physical star pattern
% Set the scaling (amount simulated beyond the detector range)
s = 2;

% Set the binning (oversampling compared to detector)
b = 2;

% Set the number of pixels
mx = b*s*2048;
my = b*s*2048;

% Set the field of view
fovx = s*23.46e-6;
fovy = s*35.16e-6;

% Generate the axes
x0 = fovx./mx.*((-mx/2):(mx/2-1)).';
y0 = fovy./my.*((-my/2):(my/2-1)).';

% Generate 2D array of axes
[X,Y] = meshgrid(x0,y0);

% Get the phase and radius
P = atan2(Y,X);
R = hypot(X,Y);

% Make the spoke pattern
sp = 24;
S = sin(sp*P);
S = double(S > 0);

% Include rings
r1 = 1.25e-6; % 1.16e-6, 1.25e-6
r2 = 2.50e-6; % 2.30e-6, 2.56e-6
r3 = 6.00e-6; % 5.66e-6, 6.42e-6
rw = 100e-9;
S(R > r1 - rw/2 & R < r1 + rw/2) = 0;
S(R > r2 - rw/2 & R < r2 + rw/2) = 0;
S(R > r3 - rw/2 & R < r3 + rw/2) = 0;
S(R < 50e-9) = 0;

% Plot the base pattern
figure;
imagesc(x0,y0,S);
axis equal tight;


%% Calculate the complex wavefront
% Photon energy
E = 17.3e3;
lambda = 1e-10*12398.42/E;
k = 2*pi/lambda;

% The refractive index decrement of Au
delta = 1.04822539e-5;

% Absorption length of Au
mu = 4.66180e-6;

% Thickness of the star
t = 600e-9;

% Calculate the complex exit field
E0 = 1e4*sqrt(exp(-t*S/mu)).*exp(-1i*k*delta*t*S);

% Make a gentle cut-off at the edges
E0 = E0.*sigmf(-X,[1e7 -22e-6]).*sigmf(X,[1e7 -22e-6]).*sigmf(-Y,[1e7 -25e-6]).*sigmf(Y,[1e7 -25e-6]);

% Plot the wave field
figure;
subplot(1,2,1);
imagesc(x0,y0,abs(E0));
axis equal tight;
colorbar;
subplot(1,2,2);
imagesc(x0,y0,angle(E0));
axis equal tight;
colorbar;


%% Save the data
% Save the exit field
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';
save([dp 'Exit_Field_1.mat'],'E0','fovx','fovy','x0','y0','s','b','-v7.3');


