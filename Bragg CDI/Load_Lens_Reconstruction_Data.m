% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the original object
tic;
O = load_binary([p 'Oxyz.bin'],[ny nx n_omega]);
toc;

% Get the object reconstructed from 1000000 counts
tic;
OL1 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_object_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 1000000 counts
tic;
OLA1 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Get the object reconstructed from 50000 counts
tic;
OL2 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_2_object_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 50000 counts
tic;
OLA2 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_2_object_avg_xyz.bin'],[ny nx n_omega]);
toc;


%% Process data
% Scale the objects to plot nicely in the range [0 1]
s0 = 6500;
tic;
OL1 = s0*OL1;
OLA1 = (s0/50)*OLA1;
OL2 = (s0*sqrt(20))*OL2;
OLA2 = (s0*sqrt(20)/50)*OLA2;
toc;

% Center the mean phase at 0
tic;
mp1 = mean(angle(OL1(OL1 ~= 0)));
mpa1 = mean(angle(OLA1(OLA1 ~= 0)));
mp2 = mean(angle(OL2(OL2 ~= 0)));
mpa2 = mean(angle(OLA2(OLA2 ~= 0)));
OL1 = OL1.*exp(-1i*mp1);
OLA1 = OLA1.*exp(-1i*mpa1);
OL2 = OL2.*exp(-1i*mp2);
OLA2 = OLA2.*exp(-1i*mpa2);
toc;

% Separate amplitude and phase
tic;
OL1a = abs(OL1);
OL1p = angle(OL1);
OLA1a = abs(OLA1);
OLA1p = angle(OLA1);
OL2a = abs(OL2);
OL2p = angle(OL2);
OLA2a = abs(OLA2);
OLA2p = angle(OLA2);
clear OL1 OLA1 OL2 OLA2;
toc;


%% Make plots
% Plot the original object
sO = Slicer(O,'slice',n_omega/2+1,'name','True object','displayRange',[0 1]);

% Plot the instantaneous 1000000 object
sOL1a = Slicer(OL1a,'slice',n_omega/2+1,'name','Object 1000000 amplitude','displayRange',[0 1]);
sOL1p = Slicer(OL1p,'slice',n_omega/2+1,'name','Object 1000000 phase','displayRange',[-0.1 0.1]);

% Plot the averaged 1000000 object
sOLA1a = Slicer(OLA1a,'slice',n_omega/2+1,'name','Object 1000000 averaged amplitude','displayRange',[0 1]);
sOLA1p = Slicer(OLA1p,'slice',n_omega/2+1,'name','Object 1000000 averaged phase','displayRange',[-0.1 0.1]);

% Plot the instantaneous 50000 object
sOL2a = Slicer(OL2a,'slice',n_omega/2+1,'name','Object 50000 amplitude','displayRange',[0 1]);
sOL2p = Slicer(OL2p,'slice',n_omega/2+1,'name','Object 50000 phase','displayRange',[-0.1 0.1]);

% Plot the averaged 50000 object
sOLA2a = Slicer(OLA2a,'slice',n_omega/2+1,'name','Object 50000 averaged amplitude','displayRange',[0 1]);
sOLA2p = Slicer(OLA2p,'slice',n_omega/2+1,'name','Object 50000 averaged phase','displayRange',[-0.1 0.1]);


