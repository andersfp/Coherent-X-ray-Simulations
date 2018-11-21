% Initialization
clear;
close all;
clc;


%% Experimental parameters
% Pt, (111) reflection, no strain, 1 um object size
%E = 8e3;
% CRL objective
% ID01 geometry - object-lens distance: 2-10 cm
% Detector: Eiger 75 um pixels, 1M, 1-6 m object-detector distance
% Include detector noise
% Include objective shape errors

% Test combination of real space image + Fourier space image


%% Load data
% Load the experimental setup
load('Exp_Param.mat');

% Load the exit field
%load([p 'Exit_Field.mat']);
E0 = load_binary([p 'E0.bin'],[ny nx n_omega]);


%% Get the expected photon number
% Pt cross section at 8 keV (m^2/atom)
sig = 65099.684e-28;

% Pt atom density
rho = 4/((3.9242e-10)^3);

% Number of Pt atoms in sample
NPt = rho*((1e-6)^3);

% Total cross section
sigtot = sig*NPt;

% Scattering probability
sp = sigtot/((1e-6)^2);

% Photon flux
pf = 1e10;

% Scattered photons
I0 = pf*1/4*sp;


%% Free space propagation
% Get the FrFT parameters for free space
R0 = Inf;
s0x = mean(diff(x))*sqrt(nx);
s0y = mean(diff(y))*sqrt(ny);
[ax,Rmx,Rpx,smx,spx,~,~] = FrFT_parameters(D,[],lambda,R0,s0x);
[ay,Rmy,Rpy,smy,spy,~,~] = FrFT_parameters(D,[],lambda,R0,s0y);

% Generate full spatial axis
[X0,Y0] = meshgrid(x,y);

% Free space propagation
Ef = zeros(size(E0));
tic;
for i = 1:n_omega
    [Ef(:,:,i),Xf,Yf] = propFrFT2(E0(:,:,i),X0,Y0,[Rmx Rmy],[Rpx Rpy],[smx smy],[spx spy],[ax ay],lambda,0,'gpu');
end
toc;

% Get 1D coordinates
xf = Xf(1,:).';
yf = Yf(:,1);

% Get intensity
If0 = abs(Ef).^2;

% Scale the intensity to get photon numbers
%sf = I0/sum(sum(If(:,:,n_omega/2+1)));
sf = 1e6/max(If0(:));
If1 = If0*sf;
If2 = If1/20;

% Convert to photon number with shot noise (poisson)
If1 = poissrnd(If1);
If2 = poissrnd(If2);

% % Shift measurement space to orthogonal space
% Ifs = zeros(size(If0));
% for i = 1:n_omega
%     Ifs(:,:,i) = circshift(If0(:,:,i),round((i - n_omega/2)*shft),1);
% end
% 
% % Plot free space propagation
% figure;
% pp = patch(isosurface(Ifs,0.001*max(Ifs(:))));
% isonormals(Ifs,pp);
% pp.FaceColor = 'red';
% pp.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight(45,45);
% camlight(-135,-45);
% lighting gouraud;
% xlabel('q_x');
% ylabel('q_y');
% zlabel('q_z');


%% CRL propagation
% Set distances
d = [d1;repmat(T,N-1,1);d4];

% Focal distances
F = repmat(f,N,1);

% FrFT parameters
[ax,Rmx,Rpx,smx,spx,~,~] = FrFT_parameters(d,F,lambda,R0,s0x);
[ay,Rmy,Rpy,smy,spy,~,~] = FrFT_parameters(d,F,lambda,R0,s0y);

% Propagation
g = @(sig,x,y) sqrt(exp(-(x.^2 + y.^2)./(2*sig.^2)));
El = zeros(size(E0));
tic;
for i = 1:n_omega
    Et = g(sigma_v,X0,Y0).*E0(:,:,i); % Effective vignetting
    [Et,Xt,Yt] = propFrFT2(Et,X0,Y0,[Rmx(1) Rmy(1)],[Inf Inf],[smx(1) smy(1)],[spx(N) spy(N)],[sum(ax(1:N)) sum(ay(1:N))],lambda,0,'gpu'); % Propagation to end of lens
    Et = g(sigma_p,Xt,Yt).*Et; % Effective pupil
    Et = Et.*exp(1i*2*pi*cos(sqrt(Xt.^2 + Yt.^2)*2*pi/5e-4)); % Phase error
    [El(:,:,i),Xl,Yl] = propFrFT2(Et,Xt,Yt,[Inf Inf],[Rpx(end) Rpy(end)],[smx(end) smy(end)],[spx(end) spy(end)],[ax(end) ay(end)],lambda,0,'gpu'); % Propagation to detector
end
toc;

% Get 1D coordinates
xl = Xl(1,:).';
yl = Yl(:,1);

% Get intensity
Il0 = abs(El).^2;

% Add attenuation from Tweb thickness
Il0 = Il0*exp(-mu*N*Tweb);

% Scale the intensity to get photon numbers
Il1 = Il0*sf;
Il2 = Il1/20;

% Convert to photon number with shot noise (poisson)
Il1 = poissrnd(Il1);
Il2 = poissrnd(Il2);

% % Shift measurement space to orthogonal space
% Ils = zeros(size(Il0));
% for i = 1:n_omega
%     Ils(:,:,i) = circshift(Il0(:,:,i),round((i - n_omega/2)*shft),1);
% end
% 
% % Plot the lens based propagation
% figure;
% pp = patch(isosurface(Ils,0.001*max(Ils(:))));
% isonormals(Ils,pp);
% pp.FaceColor = 'red';
% pp.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight(45,45);
% camlight(-135,-45);
% lighting gouraud;
% xlabel('q_x');
% ylabel('q_y');
% zlabel('q_z');


%% Save the data
% Save data
%save([p 'Diffraction_Intensity.mat'],'Ef','If','El','Il','-v7.3');
save_binary([p 'If1.bin'],If1);
save_binary([p 'If2.bin'],If2);
save_binary([p 'Ef.bin'],Ef);
save_binary([p 'Il1.bin'],Il1);
save_binary([p 'Il2.bin'],Il2);
save_binary([p 'El.bin'],El);


