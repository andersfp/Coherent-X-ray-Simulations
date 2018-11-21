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
M0 = logspace(-2,2,100);
n = length(M0);
SA = zeros(n,1);
SV = SA;
FF = SA;
D0 = SA;
DD = SA;
YN = SA;
GM = SA;
YM = SA;
for j = 1:n
M = M0(j);
d0 = fN*(1 + 1/(M*cos(N*phi)));
dd = fN*(1 + M/cos(N*phi));

% Calculate sigma_a
sigA = sqrt(R_CRL / (mu * N * (d0^2 + (f * phi)^2)))/sqrt(1 + 1 / N - 1 / (N * phi) * sin((N + 1) * phi) * cos((N - 1) * phi + 2 * atan(d0 / (f * phi))));

% Set PSF width
sigpsf = 9.861e-5;

% Set vignetting width
%sigvig = 8.293e-5;
%sigvig = Inf;
%sigvig = 1.8507e-4;
%sigvig = 1.8497e-4;
%sigvig = 187.57e-6; % M = 10
%sigvig = 242.43e-6; % M = 1
sigvig = 703.32e-6; % M = 0.1

% Vignetting from geometrical optics
sigV = delta/(mu*sigA*sqrt((N*phi)^2 - sin(N*phi)^2));

% Gamma parameter
gamma = ((f + 2*d0*N)*sqrt(d0^2 + (f*phi)^2)*cos(2*N*phi + atan(d0/(f*phi))))/(2*(d0*f + N*(d0^2 + f^2)) + (d0^2 + (f*phi)^2)*sin(2*N*phi)/phi);

% Theoretical PSF width
y0 = sigA*d0;
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
D = [d0+T/2;T*ones(N-1,1);T/2+dd];
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
tic;
for i = 1:N
    [EN,~,xN] = propFrFFT(EN,xN,xN,Inf,Inf,sm(i),sp(i),a(i),lambda,'gpu');
    EN = EN.*att(xN);
end
toc;
LN = EN(:,m/2+1);
IN = abs(EN).^2;

% Propagate to the image plane
tic;
[EM,~,xM] = propFrFFT(EN,xN,xN,Inf,Inf,sm(end),sp(end),a(end),lambda,'gpu');
toc;
LM = EM(:,m/2+1);
IM = abs(EM).^2;

% Propagate in 1 step
attvig = sqrt(exp(-(x0.^2 + x0.'.^2)./(2*(sigvig)^2)));
tic;
[E1,~,x1] = propFrFFT(E0.*attvig,x0,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,'gpu');
toc;
%attpsf = sqrt(exp(-(x1.^2 + x1.'.^2)./(2*(sigpsf)^2)));
attpsf = sqrt(exp(-(x1.^2 + x1.'.^2)./(2*(yN)^2)));
E1 = E1.*attpsf;
L1 = E1(:,m/2+1);
I1 = abs(E1).^2;

% % Propagate in 1 step
% tic;
% [E1,~,x1] = propFrFFT(E0,x0,x0,Inf,Inf,sm(1),sp(2),a(1),lambda,'gpu');
% toc;
% attvig = sqrt(exp(-(x1.^2 + x1.'.^2)./(2*(sigvig)^2)));
% tic;
% [E1,~,x1] = propFrFFT(E1.*attvig,x1,x1,Inf,Inf,sm(2),sp(end-1),sum(a(2:end-1)),lambda,'gpu');
% toc;
% attpsf = sqrt(exp(-(x1.^2 + x1.'.^2)./(2*(sigpsf)^2)));
% E1 = E1.*attpsf;
% L1 = E1(:,m/2+1);
% I1 = abs(E1).^2;

% Propagate to the image plane
tic;
[E2,~,x2] = propFrFFT(E1,x1,x1,Inf,Inf,sm(end),sp(end),a(end),lambda,'gpu');
toc;
L2 = E2(:,m/2+1);
I2 = abs(E2).^2;

% Fit the CRL plane and image plane
%yN = abs(LN).^2;
yM = abs(LM).^2;
%fun = @(a,sig,x) a.*exp(-x.^2./(2*sig.^2));
%ffN = fit(xN,yN,fun,'StartPoint',[max(yN) sqrt(sum(yN.*xN.^2)/sum(yN))]);
%ffM = fit(xM,yM,fun,'StartPoint',[max(yM) sqrt(sum(yM.*xM.^2)/sum(yM))]);

% Fit an object plane vignetting function
%fun = @(sig) Single_Prop_Fit(E0,x0,sm,sp,a,lambda,m,yM,sigpsf,sig);
fun = @(sig) Single_Prop_Fit(E0,x0,sm,sp,a,lambda,m,yM,yN,sig);
%opt = optimset('TolX',1e-9);
opt = optimset('TolX',1e-8);
tic;
%[ff,fval,exitflag,output] = fminsearch(fun,sigvig,opt);
%[ff,fval,exitflag,output] = fminsearch(fun,sigV/3,opt);
[ff,fval,exitflag,output] = fminsearch(fun,2e-7*sigA^(-0.75) + sqrt(R_CRL/(2*N*mu)),opt);
toc;

SA(j) = sigA;
SV(j) = sigV;
FF(j) = ff;
D0(j) = d0;
DD(j) = dd;
YN(j) = yN;
YM(j) = max(yn);
GM(j) = gamma;
end

%% Make plots
% Object
figure;
imagesc(1e6*x0,1e6*x0,I0);
axis equal tight;
title('Object');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (no attenuation)
figure;
imagesc(1e6*x1,1e6*x1,I1);
axis equal tight;
title('Exit field (effective attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (with attenuation)
figure;
imagesc(1e6*xN,1e6*xN,IN);
axis equal tight;
title('Exit field (lens by lens attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Image field intensity (effective attenuation)
figure;
imagesc(1e6*x2,1e6*x2,I2);
axis equal tight;
title('Image (effective attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Image field intensity (lens by lens attenuation)
figure;
imagesc(1e6*xM,1e6*xM,IM);
axis equal tight;
title('Image (lens by lens attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Compare exit fields
figure;
plot(1e6*x1,abs(L1).^2,1e6*xN,abs(LN).^2);
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
figure;
plot(1e6*x1,angle(L1*exp(-1i*angle(L1(m/2+1)))),1e6*xN,angle(LN*exp(-1i*angle(LN(m/2+1)))));
xlabel('Position (\mum)');
ylabel('Phase (rad)');

% Compare image fields
figure;
plot(1e6*x2,abs(L2).^2,1e6*xM,abs(LM).^2);
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
figure;
plot(1e6*x2,angle(L2*exp(-1i*angle(L2(m/2+1)))),1e6*xM,angle(LM*exp(-1i*angle(LM(m/2+1)))));
xlabel('Position (\mum)');
ylabel('Phase (rad)');

% Image intensity differences
figure;
plot(1e6*x2,(abs(L2).^2 - abs(LM).^2)*M^2);
xlabel('Position (\mum)');
ylabel('Intensity difference (a.u.)');


