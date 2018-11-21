% Initialization
clear;
close all;
clc;


%% Set parameters
% Load experimental parameters
load('Sim_Parameters.mat');

% Simulation name extension
ext = '';

% Sampling parameters
m = 1000;
dx = 10e-6;


%% Make the object
% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
E0 = zeros(m,m);
E0(x0==0,:) = 1;

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
Es = ((-3*max(dEmax)):1:(3*max(dEmax))) + E;
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
nE = length(lambda);
L = zeros(m,nE);
x = L;
sd = zeros(1,nE);
h = waitbar(0,'Progress');
tic;
for i = 1:nE
    waitbar(i/nE,h);
    F = repmat(f(i),N,1);
    [a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda(i),R0,s0);
    sd(i) = sp(end);
    [Et,~,xt] = propFrFFT(E0,x0,x0,Inf,Inf,sm(1),sp(N),sum(a(1:N)),lambda(i),'gpu');
    [Et,~,xt] = propFrFFT(Et.*att(xt,yN(i)),xt,xt,Inf,Inf,sm(N+1),sp(N+1),a(N+1),lambda(i),'gpu');
    L(:,i) = Et(:,xt == 0);
    x(:,i) = xt;
end
toc;
close(h);

% Save the results
save(['Pink_Beam_Profiles' ext '.mat']);



