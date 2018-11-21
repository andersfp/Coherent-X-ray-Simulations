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

% Set the distances and distance variations
%D = [9.6e-3 4.7e-3 615.7e-3];
D = [10.1e-3 3.3e-3 616.6e-3];
d = linspace(-0.9e-3,0.9e-3,7);


%% Perform the propagation
% Coherence length: lambda*L/s, ID06: L = 18 m, s = 1 mm
lambda = E2lambda(param.energy);
L = 18;
sx = 2e-3;
sy = 1e-3;
lx = lambda*L/sx;
ly = lambda*L/sy;
lx = lx/10;
ly = ly/10;

% Partial coherence parameters
sigma_fx = 2.5*lx;
sigma_fy = 2.5*ly;
sigma_rx = sqrt(4*pi*sigma_fx^4/lx^2);
sigma_ry = sqrt(4*pi*sigma_fy^4/ly^2);

% Partial coherence sampling spectrum
fx = ((-size(E0,2)/2):(size(E0,2)/2 - 1))./fovx;
fy = ((-size(E0,1)/2):(size(E0,1)/2 - 1)).'./fovy;
FF = exp(-pi.^2.*(sigma_fx.^2.*fx.^2 + sigma_fy.^2.*fy.^2));

% Number of partial coherence repetitions
p = 200;

% Calculate the stepwise propagation
n = length(d);
m = b*2048;
Id = zeros(m,m,n);
xd = zeros(m,n);
yd = xd;
img.obj = E0;
for i = 1:n
    param.D = D + d(i)*[1 0 -1];
    for j = 1:p
        if isfield(param,'lambda')
            param = rmfield(param,'lambda');
        end
        param = parameters(img,param);
        param = propagation_parameters(param);
        phi = abs(ifftshift(ifft2(FF.*randn(param.ny,param.nx))).*sqrt(sigma_rx.*param.nx.^2./fovx).*sqrt(sigma_ry.*param.ny.^2./fovy));
        img.obj = E0.*exp(1i.*phi);
        [img,param] = propagation(img,param);
        [img,param] = dezeropad(img,param);
        Id(:,:,i) = Id(:,:,i) + abs(img.det).^2;
        fprintf('.');
    end
    xd(:,i) = param.x.';
    yd(:,i) = param.y;
    fprintf('\n');
end
fprintf('\n');

% Bin the images
Id = bin3(Id,[b b 1])/p;
xd = xd(1:2:end,:);
yd = yd(1:2:end,:);

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
save([dp 'Star_Simulation_Series_1_Partial_9_200.mat'],'I','x','y','d','lx','ly','-v7.3');

% Shutdown the computer
%system('shutdown -s');

% Combine _4_1000, _5_200, _6_200, _7_200, _8_800, _9_200

