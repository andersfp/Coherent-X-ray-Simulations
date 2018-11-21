% Initialization
clear;
close all;
clc;


%% Load the phantom
% Set data path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Partial_Coherence\';

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


%% Partial coherence parameters
% Set the monochromator relative bandwidth
dEr = 0.5e-4; % 0.5e-4

% Calculate the energy bandwidth
dE = dEr*E;

% Generate an energy spectrum
nE = 31; % 41
iE = ceil(nE/2);
Es = linspace(E - 2*dE,E + 2*dE,nE);
lambdas = E2lambda(Es);

% Calculate the energy weigths
wE = gaussFWHM(Es - E,dE,0);

% Calculate the shift in reciprocal space due to energy shift
dqx = 0;
dqy = -sin(2.*th)./lambdas + sin(2.*th)./lambda;
dqz = (1-cos(2.*th))./lambdas - (1-cos(2.*th))./lambda;

% Calculate the relative pixel shift in reciprocal space
dqxr = dqx.*fov;
dqyr = dqy.*fov;
dqzr = dqz.*fov;
dqxr = permute(dqxr,[3 4 5 1 2]);
dqyr = permute(dqyr,[3 4 5 1 2]);
dqzr = permute(dqzr,[3 4 5 1 2]);

% Use the Fourier shift theorem to apply a phase ramp in real space
obj = obj.*exp((1i.*2.*pi.*dqxr./nx).*x + (1i.*2.*pi.*dqyr./ny).*y + (1i.*2.*pi.*dqzr./nz).*z);


%% Generate exit fields
% Get the number of probe positions
np2d = length(xp);

% Pre-allocate the exit fields
E0 = zeros(ny,nx,np2d,1,nE);

% Send data to the GPU
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

% Generate the exit fields
for j = 1:nE
    tmp = gpuArray(obj(:,:,:,:,j));
    for i = 1:np2d
        P = prob(xp(i),yp(i),wp,x,y,z);
        A = P.*tmp;
        B = sum(A,3);
        E0(:,:,i,:,j) = gather(B);
        fprintf('.');
    end
    fprintf('\n');
end

% Collect data from the GPU
tmp = gather(tmp);
x = gather(x);
y = gather(y);
z = gather(z);
P = gather(P);
A = gather(A);
B = gather(B);

% Remove the empty and low intensity exit fields
I = abs(E0(:,:,:,:,iE)).^2;
ii = find(squeeze(sum(sum(I)) > 1e-3.*max(sum(sum(I)))));
E0 = E0(:,:,ii,:,:);
xp = xp(ii);
yp = yp(ii);
np2d = length(ii);

% Save the data
save([p 'Exit_Fields.mat'],'E','E0','fov','G','lambda','np2d','npx','nx','ny','nz','ol','rx','ry','rz','shft','th','wp','xp','yp','x','y','z','dEr','dE','nE','iE','Es','lambdas','wE');


