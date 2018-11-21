% Initialization
clear;
close all;
clc;


%% Load data
% Set the data path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\No_Aberrations_No_Phase\';

% Load data
load([p 'Phantom.mat'],'obj','x','y','z');
load([p 'Detector_Field.mat'],'rp','ws','np','xp','yp');
load([p 'Exit_Fields.mat'],'E','lambda','nx','ny','nz','th','wp');


%% Perform projections only
% Generate the probe function
prob = @(x0,y0,wp,x,y,z) sqrt(gaussRMS(x - x0,wp).*gaussRMS(y - tan(2.*th).*z - y0,wp./cos(2.*th)));

% Generate the slit function
slit = @(r0,ws,r) rect((r - r0)./ws);
S = slit(permute(rp,[2 3 1]),ws,y);
S = flip(S,1);

% Transfer to GPU
obj = gpuArray(obj);
xp = gpuArray(xp);
yp = gpuArray(yp);
wp = gpuArray(wp);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
nx = gpuArray(nx);
ny = gpuArray(ny);
nz = gpuArray(nz);
S = gpuArray(S);

% Projections
rho = zeros(ny,nx,nz,'gpuArray');
W = rho;
for j = 1:np
    Pj = prob(xp(j),yp(j),wp,x,y,z);
    psij = obj.*Pj;
    psij = sum(psij,3)./nz;
    psij = S(:,:,j).*psij;
    psij = repmat(psij,1,1,nz);
    psij = conj(Pj).*psij;
    rho = rho + psij;
    W = W + Pj.*S(:,:,j);
end

% Transfer to CPU
obj = gather(obj);
xp = gather(xp);
yp = gather(yp);
wp = gather(wp);
x = gather(x);
y = gather(y);
z = gather(z);
nx = gather(nx);
ny = gather(ny);
nz = gather(nz);
S = gather(S);
rho = gather(rho);
W = gather(W);
Pj = gather(Pj);
psij = gather(psij);

% Normalize the result
rho = rho./W;

% Remove NaNs
rho(isnan(rho)) = 0;

% Plot the results
Slicer(obj);
Slicer(rho);

% Save the results
save('C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Projection_Probe_Only.mat','rho');



