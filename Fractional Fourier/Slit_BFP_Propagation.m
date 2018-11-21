% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters
load('Sim_Parameters.mat');

% Sampling parameters
m = 1000;
dx = 10e-6;

% Get aperture
sigA = sqrt(R_CRL / (mu * N * (d0^2 + (f * phi)^2)))/sqrt(1 + 1 / N - 1 / (N * phi) * sin((N + 1) * phi) * cos((N - 1) * phi + 2 * atan(d0 / (f * phi))));
y0 = sigA*d0;
a0 = sigA;
yN = y0*(cos(N*phi) + phi/2*sin(N*phi)) + a0*(f*phi*sin(N*phi) - T/2*cos(N*phi));


%% Make the object
% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
E0 = zeros(m,m);
E0(x0==0,x0==0) = 1;
%E0(x0==0,:) = 1;

% Define the slit size
n = 40;
sl = logspace(-4.5,-2.65,n).';


%% Calculate the wave propagation
% Make the propagation distance and focal length arrays
D = [d0+T/2;repmat(T,N-1,1);fN+T/2;dd-fN];
F = [repmat(f,N,1);Inf];

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);
[a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda,R0,s0);

% Propagate to the end of the CRL
t = initGPU();
tic;
[E3,~,x3] = propFrFFT(E0,x0,x0,Rm(1),Inf,sm(1),sp(N),sum(a(1:N)),lambda,'gpu');
toc;

% Apply effective aperture
E3 = E3.*sqrt(exp(-(x3.^2 + x3.'.^2)/(2*yN^2)));

% Propagate to the back focal plane
tic;
[E4,~,x4] = propFrFFT(E3,x3,x3,Inf,Inf,sm(N+1),sp(N+1),a(N+1),lambda,'gpu');
toc;

% Apply the slit
L = zeros(m,n);
h = waitbar(0,'Progress');
tic;
for i = 1:n
    waitbar(i/n,h);
    Ed = E4.*rect(x4/sl(i)).*rect(x4.'/sl(i));
    [Ed,~,xd] = propFrFFT(Ed,x4,x4,Inf,Rp(end),sm(end),sp(end),a(end),lambda,'gpu');
    L(:,i) = Ed(:,m/2+1);
end
toc;
close(h);

% Save the electric field interesting points
save('Slit_BFP_Propagation_Profiles.mat');

% Plot the object
figure;
imagesc(x0,x0,E0);
axis equal tight;

% Plot the CRL exit after aperture
figure;
imagesc(x3,x3,abs(E3).^2);
axis equal tight;

% Plot the back focal plane
figure;
imagesc(x4,x4,abs(E4).^2);
axis equal tight;

% Plot the last image on the detector plane
figure;
imagesc(xd,xd,abs(Ed).^2);
axis equal tight;


