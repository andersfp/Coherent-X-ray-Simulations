% Initialization
clear;
close all;
clc;

% Load the experimental parameters
load('Sim_Parameters.mat');

% Set incidence angles and locations
y1 = sigA*d0;
a1 = sigA;
y2 = 2*sigA*d0;
a2 = 2*sigA;

% Calculate the trajectories
yn1 = y1*(cos((1:N)*phi) + phi/2*sin((1:N)*phi)) + a1*(f*phi*sin((1:N)*phi) - T/2*cos((1:N)*phi));
yn2 = y2*(cos((1:N)*phi) + phi/2*sin((1:N)*phi)) + a2*(f*phi*sin((1:N)*phi) - T/2*cos((1:N)*phi));

% Plot the two trajectories
figure;
plot([0.5 1:N],[y1 yn1],[0.5 1:N],[y2 yn2]);
xlabel('Lens number');
ylabel('Ray location');

% Calculate the path length inside Be
l1 = sum(yn1.^2/R_CRL);
l2 = sum(yn2.^2/R_CRL);

% Calculate the total propagation distance
L = d0 + dd + N*T;

% Calculate Be fractions
r1 = l1/L*100;
r2 = l2/L*100;
disp([r1 r2]);

