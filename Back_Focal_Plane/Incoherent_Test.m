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
dx = 8e-6;

% Generate axis
x0 = ((-m/2):(m/2 - 1)).'./m.*dx;

% Make object
w = 2e-6;
E0 = recta(x0/w).*recta(x0.'/w);
%E0 = E0.*exp(50.*1i.*2.*pi.*abs(x0)./dx);


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


%% Incoherent propagation
% Set number of repetitions
n = 256;

% Perform propagation
EPi = zeros(m,m,n);
EBi = EPi;
EIi = EPi;
for i = 1:n
    [EPi(:,:,i),~,~] = propFrFT2(E0.*V.*exp(1i.*2.*pi.*rand(m,m)),x0.',x0,Inf,Inf,s0,sp(end-2),sum(a(1:end-2)),lambda,0,'gpu');
    [EBi(:,:,i),~,~] = propFrFT2(EPi(:,:,i).*P,xP.',xP,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda,0,'gpu');
    [EIi(:,:,i),~,~] = propFrFT2(EPi(:,:,i).*P,xP.',xP,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda,0,'gpu');
    fprintf('.');
end
fprintf('\n');

% Calculate intensities
IPI = abs(EPi).^2;
IBI = abs(EBi).^2;
III = abs(EIi).^2;

% Average intensities
IPi = mean(IPI,3);
IBi = mean(IBI,3);
IIi = mean(III,3);


%% Plots
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

% Incoherent BFP
figure;
imagesc(xB,xB,log10(IBi),log10(max(IBi(:))) + [-6 0]);
axis equal tight;
set(gca,'YDir','normal');
title('BFP incoherent');

% Incoherent image
figure;
imagesc(xI,xI,IIi);
axis equal tight;
set(gca,'YDir','normal');
title('Image incoherent');

% Coherent slices of incoherent BFP
Slicer(log10(IBI),'name','BFP','displayRange',log10(max(IBI(:))) + [-6 0]);

% Coherent slices of incoherent image
Slicer(III,'name','Image','displayRange',[0 max(III(:))]);

% Plot the incoherent BFP on linear scale
figure;
imagesc(xB,xB,IBi);
axis equal tight;
set(gca,'YDir','normal');
title('BFP incoherent (linear)');

% Calculate the autocorrelation function of the pupil
A = xcorr2fft(P,P);

% Plot the pupil autocorrelation
figure;
imagesc(xP,xP,A);
axis equal tight;
set(gca,'YDir','normal');
title('Pupil autocorrelation');


