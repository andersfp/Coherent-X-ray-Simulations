% Initialization
clear;
close all;
clc;


%% Set up the object
% Data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the Siemens star exit field
tic;
load([dp 'Exit_Field.mat']);
toc;


%% Set up propagation parameters
% Set the object parameters
param.wx = fovx;
param.wy = fovy;

% Set the detector parameters
param.det_dx = 0.74e-6/b;
param.det_dy = 0.74e-6/b;
param.det_nx = b*2048;
param.det_ny = b*2048;

% Set the distance and lens parameters
param.D = [9.6e-3 4.7e-3 615.7e-3];
param.F = [9.4537e-3 Inf;Inf 13.9754e-3]; % Mx = 64.6, My = 43.1
%param.ap = [40e-6 100e-6;100e-6 46e-6];
%param.ap_shape = [1 1];
%param.ap_shape = 0;
ew = 100e-9;
apx = @(x,y) sigmf(-x,[1/ew -20e-6]).*sigmf(x,[1/ew -20e-6]).*sigmf(-y,[1/ew -100e-6]).*sigmf(y,[1/ew -100e-6]);
apy = @(x,y) sigmf(-x,[1/ew -100e-6]).*sigmf(x,[1/ew -100e-6]).*sigmf(-y,[1/ew -23e-6]).*sigmf(y,[1/ew -23e-6]);
param.ap_shape = [5 5];
param.ap = {apx apy};

% Propagation parameters
param.phase = 0;
param.energy = 17.3e3;
param.method = 'gpu';

% Set the image struct
img.obj = E0;


%% Perform the propagation
% Calculate and check experiment parameters
param = parameters(img,param);

% Calculate FrFT parameters
param = propagation_parameters(param);

% Perform the propagation
tic;
[img,param] = propagation(img,param);
toc;

% Cut the pixels out for the detector
[img,param] = dezeropad(img,param);


%% Extract the result
% Get the axes
x = param.x;
y = param.y;
x = x(1:b:end);
y = y(1:b:end);

% Get the intensity image
I = abs(img.det).^2;
I = bin2(I,b);

% Plot the intensity
figure;
imagesc(x,y,I,[3.5e4 3.8e4]);
axis equal tight;

% Plot a line
figure;
plot(y,I(:,1));


