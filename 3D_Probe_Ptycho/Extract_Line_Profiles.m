% Initialization
clear;
close all;
clc;


%% Load data
% Set the paths
p = {...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Lens_Aberrations_Dislocations_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Lens_Aberrations_No_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\No_Aberrations_Dislocations_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\No_Aberrations_No_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_2x\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_4x\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_8x\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Partial_Coherence\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Partial_Coherence_2\',...
    };

% Select data sets
p = p([1 6 8 9]);

% Load the true object
load([p{1} 'Phantom.mat']);
nx = length(x);
ny = length(y);
nz = length(z);

% Load the reconstructions
n = length(p);
R(n) = struct('rho',[]);
for i = 1:n
    R(i) = load([p{i} 'Reconstruction.mat'],'rho');
end

% Load the probe/projection only
P = load('C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Projection_Probe_Only.mat','rho');
R(n+1).rho = P.rho;
n = n + 1;


%% Get the line profiles
% Get the x-, y-, and z-lines
lx = zeros(nx,n);
ly = zeros(ny,n);
lz = zeros(nz,n);
for i = 1:n
    lx(:,i) = R(i).rho(ny/2+1,:,nz/2+1).';
    ly(:,i) = R(i).rho(:,nx/2+1+10,nz/2+1);
    lz(:,i) = squeeze(R(i).rho(ny/2+1,nx/2+1,:));
end

% Get indices for diagonal
ipp = sub2ind([ny nx nz],1:ny,1:nx,1:nz);
ipm = sub2ind([ny nx nz],ny:-1:1,1:nx,1:nz);
imm = sub2ind([ny nx nz],ny:-1:1,nx:-1:1,1:nz);
imp = sub2ind([ny nx nz],1:ny,nx:-1:1,1:nz);

% Get the diagonal lines
dpp = zeros(nx,n);
dpm = zeros(nx,n);
dmm = zeros(nx,n);
dmp = zeros(nx,n);
for i = 1:n
    dpp(:,i) = R(i).rho(ipp);
    dpm(:,i) = R(i).rho(ipm);
    dmm(:,i) = R(i).rho(imm);
    dmp(:,i) = R(i).rho(imp);
end

% Generate the real space axes
rl = ry;
rd = sqrt(3).*ry;

% Save the line profiels
save('Line_Profiles.mat','lx','ly','lz','dpp','dpm','dmm','dmp','rl','rd');


