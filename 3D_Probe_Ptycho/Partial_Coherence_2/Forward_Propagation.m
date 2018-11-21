% Initialization
clear;
close all;
clc;


%% Load the data
% Set the data path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Partial_Coherence_2\';

% Load the data
load([p 'Exit_Fields.mat']);


%% Propagation parameters
% Get material properties
[delta,mu] = Be_Prop(Es);

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
d2 = d2(iE);
M = d2.*sin(N.*phi)./(f.*phi) - cos(N.*phi);

% Set the object to detector distance
d3 = 3;

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
%sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);

% Define the slit function
slit = @(r0,ws,r) rect((r - r0)./ws);

% Define the slit width
ws = 32;

% Calculate the image plane slit positions
nps = round(ny./ws.*1.5) + 1;
rp = linspace(-ny/2,ny/2,nps);


%% FrFT parameters
% Set distances
D = [d1 + T./2;repmat(T,N-1,1);d2 + T./2;d3];

% Focal distances
F = [repmat(f,N,1);Inf(1,nE)];

% FrFT parameters
R0 = Inf;
s0 = fov./sqrt(nx);
a = zeros(length(D),nE);
Rm = a;
Rp = a;
sm = a;
sp = a;
for i = 1:nE
    [a(:,i),Rm(:,i),Rp(:,i),sm(:,i),sp(:,i),~,~] = FrFT_parameters(D,F(:,i),lambdas(i),R0,s0);
end

% Calculate real space coordinates
x0 = rx;
x1 = sp(end-2,iE)./s0.*x0;
x2 = sp(end-1,iE)./s0.*x0;
xf = sp(end,iE)./s0.*x0;

% Define an aperture function
ap = sqrt(gaussRMS(x1,sigma_p(iE)).*gaussRMS(x1.',sigma_p(iE)));

% Calculate the aperture width
wa = sigma_p(iE)./mean(diff(x1));

% Define the aberrations
ab = exp(1i.*cos(2.*pi.*sqrt(x.^2 + y.^2)./(2.*wa)));


%% Propagation
% Propagation to CRL exit
tic;
E1 = zeros(ny,nx,np2d,1,nE);
for i = 1:nE
    E1(:,:,:,:,i) = propFrFT2(E0(:,:,:,:,i),x0,x0.',Inf,Inf,sm(1,i),sp(end-2,i),sum(a(1:end-2,i)),lambdas(i),0,'gpu');
    fprintf('.');
end
fprintf('\n');
toc;
clear E0;

% Apply the pupil function
tic;
E1a = E1.*ap.*ab;
toc;
clear E1;

% Propagation to image plane
tic;
E2 = zeros(ny,nx,np2d,1,nE);
for i = 1:nE
    E2(:,:,:,:,i) = propFrFT2(E1a(:,:,:,:,i),x1,x1.',Inf,Inf,sm(end-1,i),sp(end-1,i),a(end-1,i),lambdas(i),0,'gpu');
    fprintf('.');
end
fprintf('\n');
toc;
clear E1a;

% Apply the slit
S = slit(permute(rp,[3 4 1 2]),ws,y);
I = zeros(ny,nx,nps.*np2d);
for i = 1:nE
    E2s = E2(:,:,:,:,i).*S;
    E2s = reshape(E2s,ny,nx,nps.*np2d,1,1);
    Ef = propFrFT2(E2s,x2,x2.',Inf,Inf,sm(end,i),sp(end,i),a(end,i),lambdas(i),0,'gpu');
    clear E2s;
    I = I + wE(i).*abs(Ef).^2;
    clear Ef;
    fprintf('.');
end
fprintf('\n');
clear E2;

% Weight the intensity
I = I./sum(wE);


%% Save the results
% Generate full probe positions
xp = repmat(xp,1,nps);
yp = repmat(yp,1,nps);
rp = repmat(rp,np2d,1);
xp = xp(:);
yp = yp(:);
rp = rp(:);

% Remove the empty and low intensity fields
I0 = squeeze(sum(sum(I)));
ii = find(I0 > 1e-3.*max(I0));
I = I(:,:,ii);
xp = xp(ii);
yp = yp(ii);
rp = rp(ii);
np = length(ii);

% Save the data
tic;
save([p 'Detector_Field.mat'],'xf','rp','ws','M','wa','x','y','z','xp','yp','npx','np2d','np','nx','ny','nz','ap','ab');
toc;

% Scale the maximum intensity
sc = 1e6;
sf = sc./max(I(:));
I = I.*sf;

% Apply Poisson noise
tic;
I = poissrnd(I);
toc;

% Save the data
tic;
save([p 'Simulated_Intensity.mat'],'I','sc','sf');
toc;


