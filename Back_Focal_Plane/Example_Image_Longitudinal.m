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


%% Perform coherent propagation to get pupil function
% Calculate the vignetting function
v = gaussRMS(x0,sigma_v);
V = v.*v.';

% Generate axes
xP = x0.*(sp(end-2)./s0);
xB = x0.*(sp(end-1)./s0);
x = x0.*(sp(end)./s0);

% Calculate the pupil function
p = gaussRMS(xP,sigma_p);
P = p.*p.';


%% Partial coherence propagation
% Set the number of longitudinal repetitions
nl = 35;

% Set the energy bandwidth
E2 = linspace(-2*dE,2*dE,nl).' + E;
lambda2 = 1e-10*12398.42./E2;

% Get the intensity distribution
S = gaussFWHM(E2 - E,dE,1);

% Send data to the GPU
method = 'vec';
IP = zeros(m,m,'gpuArray');
B = IP;
I = IP;
E0 = gpuArray(E0);
x0 = gpuArray(x0);
xP = gpuArray(xP);
V = gpuArray(V);
P = gpuArray(P);

% Perform the propagation
for i = 1:nl
    [delta,~] = Be_Prop(E2(i));
    [f,~,~] = CRL_Parameters_1(R,T,N,delta);
    F = [f.*ones(N,1);Inf];
    [a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda2(i),R0,s0);
    [EP,~,xPt] = propFrFT2(E0.*V,x0.',x0,Inf,Inf,s0,sp(end-2),sum(a(1:end-2)),lambda2(i),0,method);
    [EB,~,xBt] = propFrFT2(EP.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda2(i),0,method);
    [EI,~,xIt] = propFrFT2(EP.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda2(i),0,method);
    IPt = abs(EP).^2;
    IPt = interp2(xPt.',xPt,IPt,xP.',xP,'linear',NaN);
    IP = IP + S(i).*IPt;
    Bt = abs(EB).^2;
    Bt = interp2(xBt.',xBt,Bt,xB.',xB,'linear',NaN);
    B = B + S(i).*Bt;
    It = abs(EI).^2;
    It = interp2(xIt.',xIt,It,x.',x,'linear',NaN);
    I = I + S(i).*It;
    fprintf('.');
end
fprintf('\n');

% Collect data from the GPU
EP = gather(EP);
EB = gather(EB);
EI = gather(EI);
IP = gather(IP);
B = gather(B);
I = gather(I);
E0 = gather(E0);
x0 = gather(x0);
xP = gather(xP);
V = gather(V);
P = gather(P);
IPt = gather(IPt);
Bt = gather(Bt);
It = gather(It);

% Scale intensities
IP = IP/sum(S);
B = B/sum(S);
I = I/sum(S);


%% Save the results
% Save intensity maps
save('Example_Image_Longitudinal.mat','I','B');

