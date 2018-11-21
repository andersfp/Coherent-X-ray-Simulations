% Initialization
clear;
close all;
clc;


%% Load data
% Set the full file locations
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\';
fo = 'Lens_Aberrations_Dislocations_Phase\Phantom.mat';
fp = 'Partial_Coherence\Reconstruction.mat';

% Load the data
load([p fo],'obj','rx','ry','rz','x','y','z');
par = load([p fp],'rho');
par = double(par.rho);


%% Process the data
% Normalize the data
par = 0.9*par.*(sum(abs(obj(:)))./sum(abs(par(:))));

% Adjust the phase
par = par.*exp(-1i.*angle(par(134,129,129)));

% Correct the phase parabola
par = par.*exp(-1i.*1.29922e-4.*(x.^2 + y.^2));

% Save the data
clear obj;
save([p 'Plotting_Data_Partial_Coherence.mat']);

