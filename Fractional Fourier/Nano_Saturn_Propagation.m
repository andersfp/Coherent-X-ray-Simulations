% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters
load('Sim_Parameters.mat');

% Save results?
s = 1;


%% Make the object
% Load saturn
sat = imread('saturn.png');

% Convert to grayscale
I0 = double(rgb2gray(sat));

% Set up object parameters
mx = size(I0,2);
my = size(I0,1);
dx = 5e-9*mx;
dy = 5e-9*my;

% Make coordinate
x0 = ((-mx/2):(mx/2 - 1));
x0 = x0/mx*dx;
y0 = ((-my/2):(my/2 - 1))';
y0 = y0/my*dy;

% Make amplitude of object
E0 = sqrt(I0);


%% Calculate the wave propagation
% Make the propagation distance and focal length arrays
D = [d0+T/2;T*ones(N-1,1);T/2+dd];
F = f*ones(N,1);

% Calculate the propagation parameters
R0 = Inf;
s0x = dx/sqrt(mx);
s0y = dy/sqrt(my);
[ax,Rmx,Rpx,smx,spx,gmx,gpx] = Lens_Stack(D,F,lambda,R0,s0x);
[ay,Rmy,Rpy,smy,spy,gmy,gpy] = Lens_Stack(D,F,lambda,R0,s0y);

% Initialize GPU
t = initGPU();

% Propagate in N steps with attenuation
att = @(x,y) sqrt(exp(-mu*(x.^2 + y.^2)./R_CRL));
EN = E0;
xN = x0;
yN = y0;
tic;
for i = 1:N
    [EN,xN,yN] = propFrFFT(EN,xN,yN,Inf,Inf,[smx(i) smy(i)],[spx(i) spy(i)],[ax(i) ay(i)],lambda,'gpu');
    EN = EN.*att(xN,yN);
end
toc;
IN = abs(EN).^2;

% Propagate to the image plane
tic;
[EM,xM,yM] = propFrFFT(EN,xN,yN,Inf,Inf,[smx(end) smy(end)],[spx(end) spy(end)],[ax(end) ay(end)],lambda,'gpu');
toc;
IM = abs(EM).^2;

% Propagate in 1 step
tic;
[E1,x1,y1] = propFrFFT(E0,x0,y0,Inf,Inf,[smx(1) smy(1)],[spx(end-1) spy(end-1)],[sum(ax(1:end-1)) sum(ay(1:end-1))],lambda,'gpu');
toc;
I1 = abs(E1).^2;

% Fit the effective aperture
[X,Y] = meshgrid(x1,y1);
Z = IN./I1;
fun = @(sig,x,y) exp(-(x.^2 + y.^2)./(2*(1e-4*sig).^2));
sp = 1e4*mean([sqrt(sum(sum(X.^2.*Z))./sum(sum(Z))) sqrt(sum(sum(Y.^2.*Z))./sum(sum(Z)))]);
ff = fit([X(:),Y(:)],Z(:),fun,'StartPoint',sp);

% Propagate to the image plane
att = sqrt(exp(-(x1.^2 + y1.^2)./(2*(1e-4*ff.sig)^2)));
tic;
[E2,x2,y2] = propFrFFT(E1.*att,x1,y1,Inf,Inf,[smx(end) smy(end)],[spx(end) spy(end)],[ax(end) ay(end)],lambda,'gpu');
toc;
I2 = abs(E2).^2;


%% Make plots
% Object
figure;
imagesc(1e6*x0,1e6*y0,I0);
axis equal tight;
title('Object');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (no attenuation)
figure;
imagesc(1e6*x1,1e6*y1,I1);
axis equal tight;
title('Exit field (no attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (with attenuation)
figure;
imagesc(1e6*xN,1e6*yN,IN);
axis equal tight;
title('Exit field (with attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Image field intensity (effective attenuation)
figure;
imagesc(1e6*x2,1e6*y2,I2);
axis equal tight;
title('Image (effective attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Image field intensity (lens by lens attenuation)
figure;
imagesc(1e6*xM,1e6*yM,IM);
axis equal tight;
title('Image (lens by lens attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');


%% Save data
% Save the electric field on the CRL entrance
if s == 1
    tic;
    %save([p 'Nano_Saturn_Propagation.mat'],'-v7.3');
    clear E0 E1 E2 EN EM X Y Z sat att
    save('Nano_Saturn_Propagation.mat','-v7.3');
    toc;
end



