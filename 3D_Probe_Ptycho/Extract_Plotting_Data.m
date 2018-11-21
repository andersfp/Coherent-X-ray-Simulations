% Initialization
clear;
close all;
clc;


%% Load data
% Set the full file locations
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\';
fo = 'Lens_Aberrations_Dislocations_Phase\Phantom.mat';
fr = 'Lens_Aberrations_Dislocations_Phase\Reconstruction.mat';
fn = 'NA_4x\Reconstruction.mat';

% Load the data
load([p fo],'obj','rx','ry','rz','x','y','z');
rec = load([p fr],'rho');
rec = double(rec.rho);
na4 = load([p fn],'rho');
na4 = double(na4.rho);


%% Process the data
% Normalize the data
rec = 0.9*rec.*(sum(abs(obj(:)))./sum(abs(rec(:))));
na4 = 0.9*na4.*(sum(abs(obj(:)))./sum(abs(na4(:))));
obj = 0.9*obj;

% Adjust the phase
obj = obj.*exp(-1i.*angle(obj(134,129,129)));
rec = rec.*exp(-1i.*angle(rec(134,129,129)));
na4 = na4.*exp(-1i.*angle(na4(134,129,129)));

% Correct the phase parabola
rec = rec.*exp(-1i.*1.29922e-4.*(x.^2 + y.^2));
na4 = na4.*exp(-1i.*1.29922e-4.*(x.^2 + y.^2));

% Save the data
save([p 'Plotting_Data.mat']);

