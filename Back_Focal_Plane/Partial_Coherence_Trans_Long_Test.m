% Initialization
clear;
close all;
clc;


%% Set up the experimental parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Get material properties
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6;
T = 1.6e-3;
Tweb = 2e-6;
N = 70;

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set the magnification
M = 10;
d2 = f.*phi.*(M + cos(N.*phi))./sin(N.*phi);
d1 = fN.*(d2 + f.*phi.*tan(N.*phi))./(d2 - fN);

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
sigma_v = arrayfun(@(x,y,z) Vignetting(R,x,mu,y,T,f,lambda,z),N,d1,sigma_p);


%% Generate the object
% Set the number of pixels
m = 256;

% Field-of-view
dx = 10e-6;

% Generate axis
x0 = ((-m/2):(m/2 - 1)).'./m.*dx;

% Make object
w = 8e-6;
E0 = recta(x0/w).*recta(x0.'/w);
%E0 = E0.*exp(10.*1i.*2.*pi.*abs(x0)./dx);


%% Calculate FrFT parameters
% Focal lengths
F = [f.*ones(N,1);Inf];

% Propagation distance
D = [d1 + T/2;T.*ones(N-1,1);T/2 + fN;d2 - fN];

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);

% Calculate the propagation parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);


%% Perform coherent propagation
% Calculate the vignetting function
v = gaussRMS(x0,sigma_v);
V = v.*v.';

% Propagate to the pupil plane
[EPc,~,xP] = propFrFT2(E0.*V,x0.',x0,Inf,Inf,s0,sp(end-2),sum(a(1:end-2)),lambda,0,'gpu');

% Calculate the pupil function
p = gaussRMS(xP,sigma_p);
P = p.*p.';

% Propagate to the BFP
[EBc,~,xB] = propFrFT2(EPc.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda,0,'gpu');

% Propagate to the image plane
[EIc,~,xI] = propFrFT2(EPc.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda,0,'gpu');

% Calculate intensities
IPc = abs(EPc).^2;
IBc = abs(EBc).^2;
IIc = abs(EIc).^2;


%% Partial coherence propagation
% Transversal coherence length: lambda*L/s, ID06: L = 18 m, s = 1 mm
L = 18;
s = 1e-3;
l = lambda*L/s;

% Partial transversal coherence parameters
sigma_f = 2.5*l;
sigma_r = sqrt(4*pi*sigma_f^4/l^2);

% Partial transversal coherence sampling spectrum
fx = ((-m/2):(m/2 - 1)).'./dx;
FF = exp(-pi^2*sigma_f^2*(fx.'.^2 + fx.^2));

% Set number of transversal repetitions
nt = 100;

% Set the number of longitudinal repetitions
nl = 35;

% Set the energy bandwidth
dE = 1e-4*E;
El = linspace(-2*dE,2*dE,nl).' + E;
lambda = 1e-10*12398.42./El;

% Get the intensity distribution
S = gaussFWHM(El - E,dE,1);

% Perform propagation
tic;
IPp = zeros(m,m,'gpuArray');
IBp = IPp;
IIp = IPp;
FF = gpuArray(FF);
E0 = gpuArray(E0);
x0 = gpuArray(x0);
xP = gpuArray(xP);
V = gpuArray(V);
P = gpuArray(P);
method = 'vec';
for i = 1:nl
    [a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda(i),R0,s0);
    for j = 1:nt
        phi = abs(ifftshift(ifft2(FF.*randn(m,m,'gpuArray'))).*(sigma_r.*m.^2./dx));
        [EPp,~,~] = propFrFT2(E0.*V.*exp(1i.*phi),x0.',x0,Inf,Inf,s0,sp(end-2),sum(a(1:end-2)),lambda(i),0,method);
        [EBp,~,~] = propFrFT2(EPp.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda(i),0,method);
        [EIp,~,~] = propFrFT2(EPp.*P,xP.',xP,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda(i),0,method);
        IPp = IPp + S(i).*abs(EPp).^2;
        IBp = IBp + S(i).*abs(EBp).^2;
        IIp = IIp + S(i).*abs(EIp).^2;
        fprintf('.');
    end
    fprintf('\n');
end
fprintf('\n');
FF = gather(FF);
EPp = gather(EPp);
EBp = gather(EBp);
EIp = gather(EIp);
IPp = gather(IPp);
IBp = gather(IBp);
IIp = gather(IIp);
E0 = gather(E0);
x0 = gather(x0);
xP = gather(xP);
V = gather(V);
P = gather(P);
phi = gather(phi);
toc;

% Average intensities
IPp = IPp/nt/sum(S);
IBp = IBp/nt/sum(S);
IIp = IIp/nt/sum(S);


%% Plots
% Plot the energy distribution
figure;
plot(El,S);
title('Energy distribution');

% Object
figure;
imagesc(x0,x0,abs(E0).^2);
axis equal tight;
set(gca,'YDir','normal');
title('Object');

% Vignetting function
figure;
imagesc(x0,x0,V);
axis equal tight;
set(gca,'YDir','normal');
title('Vignetting function');

% Pupil function
figure;
imagesc(xP,xP,P);
axis equal tight;
set(gca,'YDir','normal');
title('Pupil function');

% Coherent BFP
figure;
imagesc(xB,xB,log10(IBc),log10(max(IBc(:))) + [-6 0]);
axis equal tight;
set(gca,'YDir','normal');
title('BFP coherent');

% Coherent image
figure;
imagesc(xI,xI,IIc);
axis equal tight;
set(gca,'YDir','normal');
title('Image coherent');

% Partial coherence BFP
figure;
imagesc(xB,xB,log10(IBp),log10(max(IBp(:))) + [-6 0]);
axis equal tight;
set(gca,'YDir','normal');
title('BFP partial coherence');

% Partial coherence image
figure;
imagesc(xI,xI,IIp);
axis equal tight;
set(gca,'YDir','normal');
title('Image partial coherence');

% Make a combined image
figure;
subplot(2,2,1);
imagesc(xI,xI,IIc);
axis equal tight;
set(gca,'YDir','normal','OuterPosition',[0 0.5 0.5 0.5]);
title('Image coherent');
subplot(2,2,2);
imagesc(xB,xB,log10(IBc),log10(max(IBc(:))) + [-6 0]);
axis equal tight;
set(gca,'YDir','normal','OuterPosition',[0.5 0.5 0.5 0.5]);
title('BFP coherent');
subplot(2,2,3);
imagesc(xI,xI,IIp);
axis equal tight;
set(gca,'YDir','normal','OuterPosition',[0 0 0.5 0.5]);
title('Image partial coherence');
subplot(2,2,4);
imagesc(xB,xB,log10(IBp),log10(max(IBp(:))) + [-6 0]);
axis equal tight;
set(gca,'YDir','normal','OuterPosition',[0.5 0 0.5 0.5]);
title('BFP partial coherence');



