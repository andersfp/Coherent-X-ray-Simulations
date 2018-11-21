% Initialization
clear;
close all;
clc;


%% Set parameters
% Load experimental parameters
load('Sim_Parameters.mat');

% Simulation name extension
ext = '_10um';

% Sampling parameters
dx = 10e-6;


%% Make the object
% Load the image
E0 = imread('USAF_980.png');

% Process the image
E0 = double(E0)/255;
E0 = 1 - E0(:,:,1);

% Get the size of the image
m = size(E0,1);

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Set the energy RMS bandwidth
dEmax = 1e-2*E;

% Make the attenuation function
att = @(x,sig) sqrt(exp(-(x.^2 + x.'.^2)/(2*sig.^2)));

% Make the propagation distance and focal length arrays
D = [d0+T/2;repmat(T,N-1,1);dd+T/2];

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);


%% Make the propagation
% Generate energy and wavelength samples
Es = ((-3*max(dEmax)):2:(3*max(dEmax))) + E;
lambda = 1e-10*12398.42./Es;

% Get the delta and mu parameter values
[delta,mu] = Be_Prop(Es);

% Generate lens parameters
f = R_CRL./(2*delta);
phi = sqrt(T./f);
fN = f.*phi.*cot(N.*phi);

% Get aperture
sigA = sqrt(R_CRL./(mu.*N.*(d0.^2 + (f.*phi).^2)))./sqrt(1 + 1./N - 1./(N.*phi).*sin((N + 1).*phi).*cos((N - 1).*phi + 2.*atan(d0./(f.*phi))));
y0 = sigA.*d0;
a0 = sigA;
yN = y0.*(cos(N.*phi) + phi./2.*sin(N.*phi)) + a0.*(f.*phi.*sin(N.*phi) - T./2.*cos(N.*phi));

% Propagation
t = initGPU();
nE = length(lambda);
E = zeros(m,m,nE);
x = zeros(m,nE);
sd = zeros(1,nE);
h = waitbar(0,'Progress');
tic;
for i = 1:nE
    waitbar(i/nE,h);
    F = repmat(f(i),N,1);
    [a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda(i),R0,s0);
    sd(i) = sp(end);
    [Et,~,xt] = propFrFFT(E0,x0,x0,Inf,Inf,sm(1),sp(N),sum(a(1:N)),lambda(i),'gpu');
    [E(:,:,i),~,x(:,i)] = propFrFFT(Et.*att(xt,yN(i)),xt,xt,Inf,Inf,sm(N+1),sp(N+1),a(N+1),lambda(i),'gpu');
end
toc;
close(h);

% Generate the intensity
I = abs(E).^2;

% Clear variables
clear Et xt E

% Save the results
tic;
save([p 'Pink_Beam_USAF_980' ext '.mat'],'-v7.3');
toc;


