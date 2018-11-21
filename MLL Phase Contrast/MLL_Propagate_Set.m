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
title('Amplitude');
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

% Save the object
%save('M:\PostDoc Data\MLL Phase Contrast\MLL_Object.mat','u0','x','-v7.3');


%% Perform propagation
% Make lens pupils
ph = rect(X/30e-6).*rect(Y/100e-6);
pv = rect(X/100e-6).*rect(Y/30e-6);
ph = fconv(ph,g);
pv = fconv(pv,g);

% Calculate parameters
dx = FOV/(n-1);

% Set displacements
nd = 9;
d = linspace(-4e-3,4e-3,nd);

% Preallocate image
m = round(siz/pix) + 1;
I = zeros(m,m,nd);

for i = 1:nd
    d1 = 38.7508e-3 + d(i);
    d2 = dll;
    d3 = dod - d1 - d2;
    u1 = propAS(u0,lambda,dx,dx,d1);
    u1b = u1.*exp(-1i*k/(2*fh)*X.^2).*ph;
    u2 = propAS(u1b,lambda,dx,dx,d2);
    u2b = u2.*exp(-1i*k/(2*fv)*Y.^2).*pv;
    [u3,Xd,~] = propAS(u2b,lambda,dx,pix,d3);
    A = abs(u3).^2;
    xd = Xd(1,:)';
    ii = abs(xd) <= siz/2;
    xd = xd(ii);
    A = A(ii,ii);
    A = rot90(A,2);
    I(:,:,i) = A/max(max(A));
    figure;
    imagesc(1e3*xd,1e3*xd,I(:,:,i));
    axis equal tight;
    set(gca,'YDir','normal');
    title(['Detector image: ' num2str(d(i)*1e3,'%+1.1f') ' mm']);
    xlabel('x (mm)');
    ylabel('y (mm)');
end

% Save the detector images
%save('MLL_Displace_XX.mat','I','d','xd');




