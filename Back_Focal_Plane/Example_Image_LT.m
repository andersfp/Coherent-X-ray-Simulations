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
% Partial transversal coherence parameters
sigma_f = 2.5*l;
sigma_r = sqrt(4*pi*sigma_f^4/l^2);

% Partial transversal coherence sampling spectrum
fx = ((-m/2):(m/2 - 1)).'./dx;
FF = exp(-pi^2*sigma_f^2*(fx.'.^2 + fx.^2));

% Set number of transversal repetitions
nt = 100;

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
FF = gpuArray(FF);

% Perform the propagation
for i = 1:nl
    [delta,~] = Be_Prop(E2(i));
    [f,~,~] = CRL_Parameters_1(R,T,N,delta);
    F = [f.*ones(N,1);Inf];
    [a,Rm,Rp,sm,sp,~,~] = FrFT_parameters(D,F,lambda2(i),R0,s0);
    IPt = zeros(m,m,'gpuArray');
    Bt = IPt;
    It = IPt;
    for j = 1:nt
        phi = abs(ifftshift(ifft2(FF.*randn(m,m,'gpuArray'))).*(sigma_r.*m.^2./dx));
        [EP,~,xPt] = propFrFT2(E0.*V.*exp(1i.*phi),x0.',x0,Inf,Inf,s0,sp(end-2),sum(a(1:end-2)),lambda2(i),0,method);
        [EB,~,xBt] = propFrFT2(EP.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda2(i),0,method);
        [EI,~,xIt] = propFrFT2(EP.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda2(i),0,method);
        IPt = IPt + abs(EP).^2;
        Bt = Bt + abs(EB).^2;
        It = It + abs(EI).^2;
        fprintf('.');
    end
    IPt = interp2(xPt.',xPt,IPt,xP.',xP,'linear',NaN);
    IP = IP + S(i).*IPt;
    Bt = interp2(xBt.',xBt,Bt,xB.',xB,'linear',NaN);
    B = B + S(i).*Bt;
    It = interp2(xIt.',xIt,It,x.',x,'linear',NaN);
    I = I + S(i).*It;
    fprintf('\n');
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
FF = gather(FF);
phi = gather(phi);
IPt = gather(IPt);
Bt = gather(Bt);
It = gather(It);

% Scale intensities
IP = IP/nt/sum(S);
B = B/nt/sum(S);
I = I/nt/sum(S);


%% Save the results
% Save intensity maps
save('Example_Image_LT.mat','I','B');


