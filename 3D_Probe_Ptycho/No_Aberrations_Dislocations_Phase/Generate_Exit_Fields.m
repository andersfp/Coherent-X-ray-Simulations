% Initialization
clear;
close all;
clc;


%% Load the phantom
% Set data path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\No_Aberrations_Dislocations_Phase\';

% Load .mat file
load([p 'Phantom.mat']);


%% Generate the probe
% Set the photon energy and wavelength
E = 17e3;
lambda = E2lambda(E);

% Set the bragg angle
th = asin(sqrt(dot(G,G)).*lambda./(4.*pi));

% Calculate the shift amount
shft = shift_amount(mean(diff(ry)),mean(diff(rz)),th);

% Define a probe function
prob = @(x0,y0,wp,x,y,z) sqrt(gaussRMS(x - x0,wp).*gaussRMS(y - tan(2.*th).*z - y0,wp./cos(2.*th)));

% Define a 3D probe position grid
npx = 50;
nx = length(x);
ny = length(y);
nz = length(z);
xp = linspace(-nx./2,nx./2,npx);
yp = linspace(-ny./2 - nz./2.*tan(2.*th),ny./2 + nz./2.*tan(2.*th),round(npx.*(ny.*cos(2.*th) + nz.*sin(2.*th))./nx));
[xp,yp] = meshgrid(xp,yp);
xp = xp(:);
yp = yp(:);

% Determine the optimum probe size
wp = 4;
g1 = gaussRMS(x,wp);
g2 = gaussRMS(x - max(diff(xp)),wp);
gc = conv(g1,g1,'same');
gc = gc./max(gc);
ol = interp1(x,gc,max(diff(xp)));
figure;
plot(x,[g1.' g2.']);


%% Generate exit fields
% Get the number of probe positions
np2d = length(xp);

% Pre-allocate the exit fields
E0 = zeros(ny,nx,np2d);

% Send data to the GPU
obj = gpuArray(obj);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

% Generate the exit fields
for i = 1:np2d
    P = prob(xp(i),yp(i),wp,x,y,z);
    A = P.*obj;
    B = sum(A,3);
    E0(:,:,i) = gather(B);
    fprintf('.');
end
fprintf('\n');

% Collect data from the GPU
obj = gather(obj);
x = gather(x);
y = gather(y);
z = gather(z);
P = gather(P);
A = gather(A);
B = gather(B);

% Remove the empty and low intensity exit fields
I = abs(E0).^2;
ii = find(squeeze(sum(sum(I)) > 1e-3.*max(sum(sum(I)))));
E0 = E0(:,:,ii);
xp = xp(ii);
yp = yp(ii);
np2d = length(ii);

% Save the data
save([p 'Exit_Fields.mat'],'E','E0','fov','G','lambda','np2d','npx','nx','ny','nz','ol','rx','ry','rz','shft','th','wp','xp','yp','x','y','z');




