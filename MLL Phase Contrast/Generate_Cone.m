% Initialization
clear;
close all;
clc;


%% Generate cone object
% Set up physical parameters
FOV = 40e-6;

% Set up cone parameters
x0 = -7e-6;
y0 = 7e-6;
theta = 20;
n = 8192;

% Generate coordinates
x = (-n/2:n/2-1)'*FOV/n;
[X,Y] = meshgrid(x,x);

% Calculate line parameters
at = tand(-45 + theta/2);
bt = y0 - at*x0;
ab = tand(-45 - theta/2);
bb = y0 - ab*x0;
b0 = y0 + x0;

% Calculate a mask for the cone area
mask = Y > ab*X + bb & Y < at*X + bt;

% Calculate the distance from the center axis
D1 = abs(X + Y - (y0 + x0))/sqrt(2);
D2 = abs(-X + Y - (y0 + -x0))/sqrt(2);

% Calculate maximum thickness
Tmax = 2*tand(theta/2)*D2;

% Calculate thickness
T = zeros(n);
T(mask) = sqrt(Tmax(mask).^2 - D1(mask).^2);

% Plot the thickness
figure;
imagesc(x,x,T);
axis equal tight;
set(gca,'YDir','normal');
colorbar;

% Save the cone
C = T;
save('Cone.mat','C','FOV');

