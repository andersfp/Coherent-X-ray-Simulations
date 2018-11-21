% Initialization
clear;
close all;
clc;


%% Set up the object
% Data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the Siemens star exit field
tic;
load([dp 'Exit_Field_1.mat']);
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
%param.D = [9.6e-3 4.7e-3 615.7e-3];
%param.F = [9.4537e-3 Inf;Inf 13.9754e-3]; % Mx = 64.6, My = 43.1
%param.F = [9.5e-3 Inf;Inf 13.5e-3]; % Mx = 64.3, My = 44.6
param.F = [9.9e-3 Inf;Inf 13.1e-3]; % Mx = 61.4, My = 46.0
%param.ap = [40e-6 100e-6;100e-6 46e-6];
%param.ap_shape = [1 1];
aph = 37.2e-6;
apv = 43.5e-6;
ew = 100e-9;
apx = @(x,y) sigmf(-x,[1/ew -aph/2]).*sigmf(x,[1/ew -aph/2]).*sigmf(-y,[1/ew -100e-6]).*sigmf(y,[1/ew -100e-6]);
apy = @(x,y) sigmf(-x,[1/ew -100e-6]).*sigmf(x,[1/ew -100e-6]).*sigmf(-y,[1/ew -apv/2]).*sigmf(y,[1/ew -apv/2]);
param.ap_shape = [5 5];
param.ap = {apx apy};

% Propagation parameters
param.phase = 0;
param.energy = 17.3e3;
param.method = 'gpu';

% Set the image struct
img.obj = E0;

% Set the distances and distance variations
%D = [9.6e-3 4.7e-3 615.7e-3];
D = [10.1e-3 3.3e-3 616.6e-3];
d = linspace(-2e-3,2e-3,41);


%% Perform the propagation
% Calculate the stepwise propagation
n = length(d);
m = b*2048;
Id = zeros(m,m,n);
xd = zeros(m,n);
yd = xd;
for i = 1:n
    if isfield(param,'lambda')
        param = rmfield(param,'lambda');
    end
    param.D = D + d(i)*[1 0 -1];
    param = parameters(img,param);
    param = propagation_parameters(param);
    [img,param] = propagation(img,param);
    [img,param] = dezeropad(img,param);
    xd(:,i) = param.x.';
    yd(:,i) = param.y;
    Id(:,:,i) = abs(img.det).^2;
    fprintf('.');
end
fprintf('\n');

% Interpolate the simulated results
m = 2048;
det_dx = 0.74e-6;
x = ((-m/2):(m/2 - 1)).*det_dx;
y = x.';
I = zeros(m,m,n);
for i = 1:n
    I(:,:,i) = interp2(xd(:,i).',yd(:,i),Id(:,:,i),x,y,'linear');
    fprintf('.');
end
fprintf('\n');


%% Plot and save the result
% Plot the image stack
Slicer(I);

% Save the image stack
save([dp 'Star_Simulation_Series_1.mat'],'I','x','y','d','-v7.3');



