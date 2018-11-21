% Initialization
clear;
close all;
clc;


%% Load the object and setup settings
% Load the settings .mat file
load('Example_Image_Settings.mat');


%% Calculate FrFT parameters
% Focal lengths
F = [f.*ones(N,1);Inf];

% Propagation distance
D = [d1 + T/2;T.*ones(N-1,1);T/2 + fN;d2 - fN];

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);

% Calculate the propagation parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);


%% Perform coherent propagation
% Calculate the vignetting function
v = gaussRMS(x0,sigma_v);
V = v.*v.';

% Propagate to the pupil plane
[EP,~,xP] = propFrFT2(E0.*V,x0.',x0,Inf,Inf,s0,sp(end-2),sum(a(1:end-2)),lambda,0,'gpu');

% Calculate the pupil function
p = gaussRMS(xP,sigma_p);
P = p.*p.';

% Propagate to the BFP
[EB,~,xB] = propFrFT2(EP.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda,0,'gpu');

% Propagate to the image plane
[EI,~,x] = propFrFT2(EP.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda,0,'gpu');

% Calculate intensities
IP = abs(EP).^2;
B = abs(EB).^2;
I = abs(EI).^2;


%% Save the results
% Save intensity maps and coordinates
save('Example_Image_Coherent.mat','I','B','x','xB');


