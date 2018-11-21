% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters
load('Sim_Parameters.mat');

% Save results?
s = 0;


%% Make the object
% Set up object parameters
m = 1000;
dx = 1e-3;
w = 400e-6;

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
E0 = rect(x0./w).*rect(x0.'./w);

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

% Propagate in 1 step
tic;
[E1,~,x1] = propFrFFT(E0,x0,x0,Rm(1),Rm(end),sm(1),sp(end-1),sum(a(1:end-1)),lambda,'gpu');
toc;
L1 = E1(:,m/2+1);
I1 = abs(E1).^2;
if s == 1
    Save_Tiff([p 'Effective_Aperture_No_Attenuation_I1.tiff'],I1,max(max(I1)));
end

% Propagate to the image plane
tic;
[E2,~,x2] = propFrFFT(E1,x1,x1,Rm(end),Rp(end),sm(end),sp(end),a(end),lambda,'gpu');
toc;
L2 = E2(:,m/2+1);
I2 = abs(E2).^2;
if s == 1
    Save_Tiff([p 'Effective_Aperture_No_Attenuation_I2.tiff'],I2,max(max(I2)));
end

% Propagate in N steps without attenuation
tic;
[EN,~,xN] = propFrFFT(E0,x0,x0,Inf,Inf,sm(1),sp(1),a(1),lambda,'gpu');
for i = 2:N
    [EN,~,xN] = propFrFFT(EN,xN,xN,Inf,Inf,sm(i),sp(i),a(i),lambda,'gpu');
end
[EN,~,xN] = propFrFFT(EN,xN,xN,Inf,Rm(end),sp(end-1),sp(end-1),0,lambda,'gpu');
toc;
LN = EN(:,m/2+1);
IN = abs(EN).^2;
if s == 1
    Save_Tiff([p 'Effective_Aperture_No_Attenuation_IN.tiff'],IN,max(max(IN)));
end

% Propagate to the image plane
tic;
[EM,~,xM] = propFrFFT(EN,xN,xN,Rm(end),Rp(end),sm(end),sp(end),a(end),lambda,'gpu');
toc;
LM = EM(:,m/2+1);
IM = abs(EM).^2;
if s == 1
    Save_Tiff([p 'Effective_Aperture_No_Attenuation_IM.tiff'],IM,max(max(IM)));
end

% Save the electric field on the CRL entrance
if s == 1
    tic;
    save([p 'Effective_Aperture_No_Attenuation.mat'],'E0','x0','L0','E1','x1','L1','E2','x2','L2','EN','xN','LN','EM','xM','LM','-v7.3');
    save('Effective_Aperture_No_Attenuation_Profiles.mat','L0','x0','L1','x1','L2','x2','LN','xN','LM','xM','-v7.3');
    toc;
end


%% Make plots
% Object
figure;
imagesc(1e6*x0,1e6*x0,I0);
axis equal tight;
title('Object');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (single propagation)
figure;
imagesc(1e6*x1,1e6*x1,I1);
axis equal tight;
title('Exit field (single propagation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (N propagations)
figure;
imagesc(1e6*xN,1e6*xN,IN);
axis equal tight;
title('Exit field (N propagations)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Absolute difference map
figure;
imagesc(1e6*xN,1e6*xN,IN-I1);
axis equal tight;
title('Absolute difference');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Relative difference map
figure;
imagesc(1e6*xN,1e6*xN,(IN-I1)./I1);
axis equal tight;
title('Relative difference');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Image intensity (2 propagations)
figure;
imagesc(1e6*x2,1e6*x2,I2);
axis equal tight;
title('Image (single propagation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Image intensity (N+1 propagations)
figure;
imagesc(1e6*xM,1e6*xM,IM);
axis equal tight;
title('Image (N propagations)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Absolute difference map
figure;
imagesc(1e6*xM,1e6*xM,IM-I2);
axis equal tight;
title('Absolute difference');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Relative difference map
figure;
imagesc(1e6*xM,1e6*xM,(IM-I2)./I2);
axis equal tight;
title('Relative difference');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

