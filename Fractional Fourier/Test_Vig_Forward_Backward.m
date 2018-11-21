% Initialization
clear;
close all;
clc;


%% Set parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42/E;

% Set CRL parameters
R_CRL = 50e-6;
delta = 1.17843774e-6;
f = R_CRL/(2*delta);
T = 1.6e-3;
N = 69;
mu = 1/24590.6e-6;
Tweb = 2e-6;

% Calculate focal length
phi = sqrt(T/f);
fN = f*phi*cot(N*phi);

% Set object and image positions (see CRL_Focal_Length_Sim)
%M0 = (1:1:100)*0.1;
M0 = logspace(-1.3,3,50);
%M0 = 0.04;
n = length(M0);
SA = zeros(n,1);
SV = SA;
FF = SA;
D1 = SA;
D2 = SA;
YN = SA;
GM = SA;
YM = SA;
h = waitbar(0,'Progress');
for j = 1:n
waitbar(j/n,h);
M = M0(j);
d1 = fN*(1 + 1/(M*cos(N*phi)));
d2 = fN*(1 + M/cos(N*phi));

% Calculate sigma_a
sigA = sqrt(R_CRL / (mu * N * (d1^2 + (f * phi)^2)))/sqrt(1 + 1 / N - 1 / (N * phi) * sin((N + 1) * phi) * cos((N - 1) * phi + 2 * atan(d1 / (f * phi))));

% Set PSF width
sigpsf = 9.861e-5;

% Set vignetting width
%sigvig = 8.293e-5;
%sigvig = Inf;
%sigvig = 1.8507e-4;
%sigvig = 1.8497e-4;
sigvig = 187.57e-6; % M = 10
%sigvig = 242.43e-6; % M = 1
%sigvig = 703.32e-6; % M = 0.1

% Vignetting from geometrical optics
sigV = delta/(mu*sigA*sqrt((N*phi)^2 - sin(N*phi)^2));

% Gamma parameter
gamma = ((f + 2*d1*N)*sqrt(d1^2 + (f*phi)^2)*cos(2*N*phi + atan(d1/(f*phi))))/(2*(d1*f + N*(d1^2 + f^2)) + (d1^2 + (f*phi)^2)*sin(2*N*phi)/phi);

% Theoretical PSF width
y0 = sigA*d1;
a0 = sigA;
yn = y0*(cos((1:N)*phi) + phi/2*sin((1:N)*phi)) + a0*(f*phi*sin((1:N)*phi) - T/2*cos((1:N)*phi));
yN = yn(end);

% Effective aperture
sigD = sqrt(R_CRL/(mu*N))/sqrt(1 + sinc(2*N*phi));


%% Make the object
% Set up object parameters
m = 1000;
dx = 2000e-6;

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
w = 1000e-6;
E0 = rect(x0/w).*rect(x0.'/w);

% Make line and intensity from object
L0 = E0(:,m/2+1);
I0 = abs(E0).^2;


%% Calculate the wave propagation
% Make the propagation distance and focal length arrays
D = [d1+T/2;T*ones(N-1,1);T/2+d2];
F = f*ones(N,1);

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);
[a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda,R0,s0);

% Initialize GPU
t = initGPU();

% Propagate in N steps with attenuation
att = @(x) sqrt(exp(-mu*(x.^2 + x.'.^2)./R_CRL));
EN = E0;
xN = x0;
%tic;
for i = 1:N
    [EN,~,xN] = propFrFFT(EN,xN,xN,Inf,Inf,sm(i),sp(i),a(i),lambda,'gpu');
    EN = EN.*att(xN);
end
%toc;
LN = EN(:,m/2+1);
IN = abs(EN).^2;

% Propagate backward in one step
attpsf = sqrt(exp(-(xN.^2 + xN.'.^2)./(2*(yN)^2)));
mask = rect(xN./w).*rect(xN.'./w);
attpsf(mask < 0.1) = 1;
%tic;
[E1,~,x1] = propFrFFT(EN./attpsf,xN,xN,Inf,Inf,sp(N),sm(1),-sum(a(1:N)),lambda,'gpu');
%toc;
L1 = E1(:,m/2+1);
I1 = abs(E1).^2;

% Fit the CRL plane and image plane
x = x1;
y = abs(L1).^2./abs(L0).^2;
fun = @(sig,x) exp(-x.^2./(2*sig.^2));
ff = fit(x,y,fun,'StartPoint',100e-6,'Lower',0,'Exclude',isinf(y),'Robust','bisquare');

figure;
plot(ff,x,y);

SA(j) = sigA;
SV(j) = sigV;
FF(j) = ff.sig;
D1(j) = d1;
D2(j) = d2;
YN(j) = yN;
YM(j) = max(yn);
GM(j) = gamma;
end
close(h);


%% Make plots
% Object
figure;
imagesc(1e6*x0,1e6*x0,I0);
axis equal tight;
title('Object intensity');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (lens-by-lens)
figure;
imagesc(1e6*xN,1e6*xN,IN);
axis equal tight;
title('Exit intensity (lens by lens attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Object field intensity (back-propagation)
figure;
imagesc(1e6*x1,1e6*x1,I1);
axis equal tight;
title('Object intensity (back-propagation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');



