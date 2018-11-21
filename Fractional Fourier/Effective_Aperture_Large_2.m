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
% Set up object parameters
m = 4096;
dx = 1e-3;
w = Inf; % dx*2/5; % 600e-6

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
E0 = rect(x0./w).*rect(x0.'./w);

% Make line and intensity from object
L0 = E0(:,m/2+1);
I0 = abs(E0).^2;

% Calculate analytical apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R_CRL,mu,f,phi,d0);

% Calculate the effective vignetting width
sigma_v = Vignetting(R_CRL,N,mu,d0,T,f,lambda,sigma_p);


%% Calculate the wave propagation
% Make the propagation distance and focal length arrays
D = [d0+T/2;T*ones(N-1,1);T/2+dd];
F = f*ones(N,1);

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);

% Optimize the d2 distance
%d2c = optimImagePos(L0,x0,D,F,lambda,R0,s0,sigma_v,sigma_p);
d2c = -0.0241;
%d2c = 0;
D(end) = D(end) + d2c;

% Calculate the FrFT parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);

% Propagate in 1 step
tic;
[E1,~,x1] = propFrFT2(E0,x0,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,0,'gpu');
toc;
L1 = E1(:,m/2+1);
I1 = abs(E1).^2;

% Propagate to the image plane
att = sqrt(exp(-(x1.^2 + x1.'.^2)/(2*(8.293e-5)^2)));
tic;
[E2,~,x2] = propFrFT2(E1.*att,x1,x1,Inf,Rp(end),sm(end),sp(end),a(end),lambda,0,'gpu');
toc;
L2 = E2(:,m/2+1);
I2 = abs(E2).^2;

% Propagate in N steps with attenuation
att = @(x) sqrt(exp(-mu*(x.^2 + x.'.^2)./R_CRL));
tic;
[EN,~,xN] = propFrFT2(E0,x0,x0,Inf,Inf,sm(1),sp(1),a(1),lambda,0,'gpu');
EN = EN.*att(xN);
for i = 2:N
    [EN,~,xN] = propFrFT2(EN,xN,xN,Inf,Inf,sm(i),sp(i),a(i),lambda,0,'gpu');
    EN = EN.*att(xN);
end
[EN,~,xN] = propFrFT2(EN,xN,xN,Inf,Inf,sp(end-1),sp(end-1),0,lambda,0,'gpu');
toc;
LN = EN(:,m/2+1);
IN = abs(EN).^2;

% Propagate to the image plane
tic;
[EM,~,xM] = propFrFT2(EN,xN,xN,Inf,Rp(end),sm(end),sp(end),a(end),lambda,0,'gpu');
toc;
LM = EM(:,m/2+1);
IM = abs(EM).^2;


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
title('Exit field (no attenuation)');
xlabel('Position (\mum)');
ylabel('Position (\mum)');

% Exit field intensity (with attenuation)
figure;
imagesc(1e6*xN,1e6*xN,IN);
axis equal tight;
title('Exit field (with attenuation)');
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


%% Save data
% Save the electric field on the CRL entrance
if s == 1
    tic;
    %save([p 'Effective_Aperture_Large_2.mat'],'E0','x0','L0','E1','x1','L1','E2','x2','L2','EN','xN','LN','EM','xM','LM','-v7.3');
    clear E0 E1 E2 EN EM I0 I1 I2 IN IM;
    save('Effective_Aperture_Large_Profiles_2.mat','-v7.3');
    toc;
end



