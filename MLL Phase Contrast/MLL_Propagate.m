% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the cone thickness
load('Cone.mat');

% Set physical parameters
E = 30e3;
lambda = 1e-10*12398.42/E;
k = 2*pi/lambda;

% Set lens parameters
fh = 38e-3;
fv = 43e-3;

% Set distances
dll = 5e-3;
dod = 2;

% Detector and pixel size
pix = 1e-6;
siz = 2e-3;

% Zeropadding factor
zp = 2;


%% Generate object field
% Material parameters
Au.delta.E17 = 1.08325503e-5;
Au.delta.E30 = 3.55550856e-6;
Au.mu.E17 = 4.45817e-6;
Au.mu.E30 = 19.9816e-6;
GaAs.delta.E17 = 3.393101e-6;
GaAs.delta.E30 = 1.09541213e-6;
GaAs.mu.E17 = 28.0344e-6;
GaAs.mu.E30 = 136.143e-6;

% Pick the correct parameter
mu = Au.mu.E30;
delta = Au.delta.E30;

% Generate coordinates
n = size(C,1);
x0 = (-n/2:n/2-1)'*FOV/n;

% Plot the thickness
figure;
imagesc(x0,x0,C);
axis equal tight;
set(gca,'YDir','normal');
title('Thickness');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Calculate the complex object field
u0 = sqrt(exp(-C/mu)).*exp(-1i*k*delta*C);

% Zeropad the object
N = zp*n;
u0 = padarray(u0,[(N-n)/2 (N-n)/2],1,'both');

% Make updated coordinates
x = (-N/2:N/2-1)'*FOV/n;
[X,Y] = meshgrid(x,x);

% Make fuzzy support boundary
p = rect(X/FOV).*rect(Y/FOV);
g = exp(-(X.^2 + Y.^2)/(2*(1e-6)^2));
p = fconv(p,g);
u0 = u0.*p;

% Plot the object amplitude and phase
figure;
subplot(1,2,1);
imagesc(1e6*x,1e6*x,abs(u0).^2);
title('Intensity');
xlabel('x (\mum)');
ylabel('y (\mum)');
axis equal tight;
set(gca,'YDir','normal');
subplot(1,2,2);
imagesc(1e6*x,1e6*x,angle(u0));
title('Phase');
xlabel('x (\mum)');
ylabel('y (\mum)');
axis equal tight;
set(gca,'YDir','normal');


%% Perform propagation
% Make lens pupils
ph = rect(X/30e-6).*rect(Y/100e-6);
pv = rect(X/100e-6).*rect(Y/30e-6);
ph = fconv(ph,g);
pv = fconv(pv,g);

% Make distances
d1 = 38.7508e-3 - 2e-3;
d2 = dll;
d3 = dod - d1 - d2;

% Calculate propagation parameters
fn1 = FOV^2/(lambda*d1);
fn2 = FOV^2/(lambda*d2);
fn3 = FOV^2/(lambda*d3);
disp(['Fresnel numbers: ' num2str(fn1) ' ' num2str(fn2) ' ' num2str(fn3)]);
dx = FOV/(n-1);
disp(['dx and lambda*z/L: ' num2str(dx) ' ' num2str(lambda*d1/FOV) ' ' num2str(lambda*d2/FOV) ' ' num2str(lambda*d3/FOV)]);
fny = 1/(2*dx);
disp(['Nyquist frequency: ' num2str(fny)]);

% Propagate to first lens
u1 = propAS(u0,lambda,dx,dx,d1);

% Apply the lens phase shift and the pupil
u1b = u1.*exp(-1i*k/(2*fh)*X.^2).*ph;

% Propagate to the second lens
u2 = propAS(u1b,lambda,dx,dx,d2);

% Apply the lens phase shift and the pupil
u2b = u2.*exp(-1i*k/(2*fv)*Y.^2).*pv;

% Propagate to the detector
[u3,Xd,Yd] = propAS(u2b,lambda,dx,pix,d3);

% Plot the intermediate results
figure;
subplot(2,2,1);
imagesc(1e6*x,1e6*x,abs(u1));
axis equal tight;
set(gca,'YDir','normal');
title('Before horizontal lens');
xlabel('x (\mum)');
ylabel('y (\mum)');
subplot(2,2,2);
imagesc(1e6*x,1e6*x,abs(u1b));
axis equal tight;
set(gca,'YDir','normal');
title('After horizontal lens');
xlabel('x (\mum)');
ylabel('y (\mum)');
subplot(2,2,3);
imagesc(1e6*x,1e6*x,abs(u2));
axis equal tight;
set(gca,'YDir','normal');
title('Before vertical lens');
xlabel('x (\mum)');
ylabel('y (\mum)');
subplot(2,2,4);
imagesc(1e6*x,1e6*x,abs(u2b));
axis equal tight;
set(gca,'YDir','normal');
title('After vertical lens');
xlabel('x (\mum)');
ylabel('y (\mum)');

% Calculate the intensity
I = abs(u3).^2;

% Flip the image 180 degrees
I = rot90(I,2);

% Normalize the image
I = I/max(max(I));

% Extract the detector coordinates
xd = Xd(1,:);
yd = Yd(:,1);

% Extract the active detector area
ix = abs(xd) <= siz/2;
iy = abs(yd) <= siz/2;
xd = xd(ix);
yd = yd(iy);
I = I(iy,ix);

% Plot the intensity
figure;
imagesc(1e3*xd,1e3*yd,I);
axis equal tight;
set(gca,'YDir','normal');
title('Detector image');
xlabel('x (mm)');
ylabel('y (mm)');


