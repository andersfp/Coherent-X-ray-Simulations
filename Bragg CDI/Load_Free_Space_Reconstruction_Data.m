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
OF1 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_object_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 1000000 counts
tic;
OFA1 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Get the object reconstructed from 50000 counts
tic;
OF2 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 50000 counts
tic;
OFA2 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_avg_xyz.bin'],[ny nx n_omega]);
toc;


%% Process data
% Scale the objects to plot nicely in the range [0 1]
s0 = 6500;
tic;
OF1 = s0*OF1;
OFA1 = (s0/50)*OFA1;
OF2 = (s0*sqrt(20))*OF2;
OFA2 = (s0*sqrt(20)/50)*OFA2;
toc;

% Center the mean phase at 0
tic;
OF2 = OF2.*exp(1i*2.0);
mp1 = mean(angle(OF1(OF1 ~= 0)));
mpa1 = mean(angle(OFA1(OFA1 ~= 0)));
mp2 = mean(angle(OF2(OF2 ~= 0)));
mpa2 = mean(angle(OFA2(OFA2 ~= 0)));
OF1 = OF1.*exp(-1i*mp1);
OFA1 = OFA1.*exp(-1i*mpa1);
OF2 = OF2.*exp(-1i*mp2);
OFA2 = OFA2.*exp(-1i*mpa2);
toc;

% Separate amplitude and phase
tic;
OF1a = abs(OF1);
OF1p = angle(OF1);
OFA1a = abs(OFA1);
OFA1p = angle(OFA1);
OF2a = abs(OF2);
OF2p = angle(OF2);
OFA2a = abs(OFA2);
OFA2p = angle(OFA2);
clear OF1 OFA1 OF2 OFA2;
toc;


%% Make plots
% Plot the original object
sO = Slicer(O,'slice',n_omega/2+1,'name','True object','displayRange',[0 1]);

% Plot the instantaneous 1000000 object
sOF1a = Slicer(OF1a,'slice',n_omega/2+1,'name','Object 1000000 amplitude','displayRange',[0 1]);
sOF1p = Slicer(OF1p,'slice',n_omega/2+1,'name','Object 1000000 phase','displayRange',[-0.1 0.1]);

% Plot the averaged 1000000 object
sOFA1a = Slicer(OFA1a,'slice',n_omega/2+1,'name','Object 1000000 averaged amplitude','displayRange',[0 1]);
sOFA1p = Slicer(OFA1p,'slice',n_omega/2+1,'name','Object 1000000 averaged phase','displayRange',[-0.1 0.1]);

% Plot the instantaneous 50000 object
sOF2a = Slicer(OF2a,'slice',n_omega/2+1,'name','Object 50000 amplitude','displayRange',[0 1]);
sOF2p = Slicer(OF2p,'slice',n_omega/2+1,'name','Object 50000 phase','displayRange',[-0.1 0.1]);

% Plot the averaged 50000 object
sOFA2a = Slicer(OFA2a,'slice',n_omega/2+1,'name','Object 50000 averaged amplitude','displayRange',[0 1]);
sOFA2p = Slicer(OFA2p,'slice',n_omega/2+1,'name','Object 50000 averaged phase','displayRange',[-0.1 0.1]);


