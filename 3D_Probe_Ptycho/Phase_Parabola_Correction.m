% Initialization
clear;
close all;
clc;


%% Load the data
% Set the data path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Lens_Aberrations_Dislocations_Phase\';

% Load the data
load([p 'Exit_Fields.mat']);


%% Propagation parameters
% Get material properties
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6;
T = 1.6e-3;
Tweb = 2e-6;
N = 70;

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set object and image positions
d1 = 0.42;
d2 = fN.*(d1 + f.*phi.*tan(N.*phi))./(d1 - fN);
M = d2.*sin(N.*phi)./(f.*phi) - cos(N.*phi);

% Set the object to detector distance
d3 = 3;

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);


%% FrFT parameters
% Set distances
D = [d1 + T./2;repmat(T,N-1,1);d2 + T./2;d3];

% Focal distances
F = [repmat(f,N,1);Inf];

% FrFT parameters
R0 = Inf;
s0 = fov./sqrt(nx);
[a,Rm,Rp,sm,sp,~,~] = FrFT_parameters(D,F,lambda,R0,s0);

% Calculate real space coordinates
x0 = ry;
x1 = sp(end-2)./s0.*x0;
x2 = sp(end-1)./s0.*x0;
x3 = sp(end)./s0.*x0;

% Define an aperture function
ap = sqrt(gaussRMS(x1,sigma_p));


%% Propagate test object
% Make test object
E0 = recta(x0./(fov.*0.05));

% Propagate to the CRL exit
E1 = propFrFT1(E0,x0,Inf,Inf,sm(1),sp(end-2),sum(a(1:end-2)),lambda,0);

% Apply the aperture
E1a = E1.*ap;

% Propagate to the image plane
E2 = propFrFT1(E1a,x1,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda,0);

% Propagate directly to the image plane
E2d = propFrFT1(E0,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,0);

% Propagate to the detector plane
E3 = propFrFT1(E2,x2,Inf,Inf,sm(end),sp(end),a(end),lambda,0);

% Propagate semi-directly to the detector plane
E3sd = propFrFT1(E1a,x1,Inf,Inf,sm(end-1),sp(end),sum(a(end-1:end)),lambda,0);

% Propagate directly to the detector plane
E3d = propFrFT1(E0,x0,Inf,Inf,sm(1),sp(end),sum(a(1:end)),lambda,0);

% Calculate reference for detector plane
E3r = flip(fftshift(fft(fftshift(E0))));


%% Plot the results
% Scale factor
s = 1e6;

% Plot the object
figure;
plot(s.*x0,E0);

% Plot the image plane
figure;
plotyy(s.*x2,abs(E2),s.*x2,angle(E2));

% Plot the direct image plane
figure;
plotyy(s.*x2,abs(E2d),s.*x2,angle(E2d));

% Detector plane direct
figure;
plotyy(s.*x3,abs(E3d),s.*x3,angle(E3d));

% Detector plane reference
figure;
plotyy(s.*x3,abs(E3r),s.*x3,angle(E3r));

% Plot the detector plane ratio
figure;
plotyy(s.*x3,abs(E3d.*conj(E3r).*y),s.*x3,angle(E3d.*conj(E3r)));


%% Fit the curvature
% Define fitting function
fun = @(a,x) 1e-5.*a.*x.^2;

% Extract the phase
ph2 = angle(E2d);
ph2 = ph2 - ph2(ny/2 + 1);

% Fit the phase
ff2 = fit(y,ph2,fun,'StartPoint',-7.5,'Exclude',abs(y) > 5);

% Plot the fit
figure;
plot(ff2,y,ph2);

% Extract phase
ph3 = angle(E3d.*conj(E3r));
ph3 = ph3 - ph3(ny/2 + 1);
w = abs(E3d.*conj(E3r).*y);

% Fit the phase
ff3 = fit(y,ph3,fun,'StartPoint',13,'Exclude',w < 0.05 | abs(y) > 120);

% Plot the fit
figure;
plot(ff3,y,ph3);


