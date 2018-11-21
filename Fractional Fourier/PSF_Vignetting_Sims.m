% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters
load('Sim_Parameters.mat');

% Set PSF width
sigpsf = 9.861e-5;

% Set vignetting width
sigvig = 1.8507e-4;


%% Make the object
% Set up object parameters
m = 1000;
dx = 3e-6;

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
w = 1e-6;
E0 = rect(x0/w).*rect(x0.'/w);


%% Calculate the wave propagation
% Set up side lengths
DX = [5 10 20 40 80 160 320 640 1280 2560]*1e-6;

% Make the propagation distance and focal length arrays
D = [d0+T/2;T*ones(N-1,1);T/2+dd];
F = f*ones(N,1);

% Initialize GPU
t = initGPU();

% Run the simulation for each object size
n = length(DX);
L2 = zeros(m,n);
L70 = L2;
for i = 1:n
    dx = DX(i);
    
    % Make coordinate
    x0 = ((-m/2):(m/2 - 1))';
    x0 = x0/m*dx;
    
    % Calculate the propagation parameters
    R0 = Inf;
    s0 = dx/sqrt(m);
    [a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda,R0,s0);
    
    % Propagate in N steps with attenuation
    att = @(x) sqrt(exp(-mu*(x.^2 + x.'.^2)./R_CRL));
    EN = E0;
    xN = x0;
    tic;
    for j = 1:N
        [EN,~,xN] = propFrFFT(EN,xN,xN,Inf,Inf,sm(j),sp(j),a(j),lambda,'gpu');
        EN = EN.*att(xN);
    end
    toc;
    
    % Propagate to the image plane
    tic;
    [EM,~,xM] = propFrFFT(EN,xN,xN,Inf,Inf,sm(end),sp(end),a(end),lambda,'gpu');
    toc;
    L70(:,i) = EM(:,m/2+1);
    
    % Propagate in 1 step
    attvig = sqrt(exp(-(x0.^2 + x0.'.^2)./(2*(sigvig)^2)));
    tic;
    [E1,~,x1] = propFrFFT(E0.*attvig,x0,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,'gpu');
    toc;
    attpsf = sqrt(exp(-(x1.^2 + x1.'.^2)./(2*(sigpsf)^2)));
    E1 = E1.*attpsf;
    
    % Propagate to the image plane
    tic;
    [E2,~,x2] = propFrFFT(E1,x1,x1,Inf,Inf,sm(end),sp(end),a(end),lambda,'gpu');
    toc;
    L2(:,i) = E2(:,m/2+1);
end


%% Make plots
% Get intensity
I2 = abs(L2).^2;
I70 = abs(L70).^2;

% Get the phase
% P2 = unwrap(angle(L2.*exp(-1i*angle(L2(m/2+1,:)))));
% P2 = P2 - P2(m/2+1,:);
% P70 = unwrap(angle(L70.*exp(-1i*angle(L70(m/2+1,:)))));
% P70 = P70 - P70(m/2+1,:);
P2 = angle(L2.*exp(-1i*angle(L2(m/2+1,:))));
P2 = angle(L2.*exp(-1i*(mean(P2(251:751,:)) + angle(L2(m/2+1,:)))));
P70 = angle(L70.*exp(-1i*angle(L70(m/2+1,:))));
P70 = angle(L70.*exp(-1i*(mean(P70(251:751,:)) + angle(L70(m/2+1,:)))));

% Generate legends
leg = cellfun(@num2str,num2cell(DX*1e6),'UniformOutput',false);

% Plot the intensity
figure;
plot(I2,'-');
legend(leg);
hold on;
set(gca,'ColorOrderIndex',1);
plot(I70,'--');

% Plot the phase
figure;
plot(P2,'-');
legend(leg);
hold on;
set(gca,'ColorOrderIndex',1);
plot(P70,'--');

