% Initialization
clear;
close all;
clc;


%% Load data
% Set the full file locations
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\';
f = 'Plotting_Data.mat';

% Load the data
load([p f]);


%% Extract the phase aberrations
% Get the phases
ii = 129;
pobj = angle(obj(:,:,ii));
prec = angle(rec(:,:,ii));
pna4 = angle(na4(:,:,ii));

% Get the phase differences
drec = prec - pobj;
dna4 = pna4 - pobj;

% Unwrap the phase differences
drec(drec < -pi) = drec(drec < -pi) + 2*pi;
drec(drec > pi) = drec(drec > pi) - 2*pi;
dna4(dna4 < -pi) = dna4(dna4 < -pi) + 2*pi;
dna4(dna4 > pi) = dna4(dna4 > pi) - 2*pi;

% Plot the phase difference
figure;
subplot(1,2,1);
imagesc(drec,[-pi pi]);
axis equal tight;
subplot(1,2,2);
imagesc(dna4,[-pi pi]);
axis equal tight;

% Plot a single phase difference
drec(abs(rec(:,:,ii)) < 0.2) = -pi;
figure;
imagesc(drec,[-pi pi]);
axis equal tight off;
colorbar;



