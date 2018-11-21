% Initialization
clear;
close all;
clc;


%% Material parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Lattice constant for Si
a = 543.102e-12;

% Momentum transfer for Si(220) reflektion
q = 2*pi/a*sqrt(8);

% Bragg angle
th_bragg = asin(q.*lambda./(4.*pi));


%% Standard BCDI setup
% Detector
det_dx = 1e-6;
nx = 512;
ny = 512;

% Sample detector distance
D = 27.4228;

% Omega scan
tth = 2*th_bragg;
omega_cen = th_bragg;
delta_omega = 5.6028e-6*pi/180;
n_omega = 512;

% Calculate q-vectors
[q1v,q2v,q3v,r1v,r2v,r3v,qxv,qyv,qzv,rxv,ryv,rzv,beta] = q_vector(D,det_dx,tth,lambda,omega_cen,delta_omega,nx,ny,n_omega);

% Calculate q-vector length
delta_q1 = norm(q1v);
delta_q2 = norm(q2v);
delta_q3 = norm(q3v);
delta_qx = norm(qxv);
delta_qy = norm(qyv);
delta_qz = norm(qzv);

% Generate q-axes
qx = (-nx/2:(nx/2-1)).'.*delta_qx;
qy = (-ny/2:(ny/2-1)).'.*delta_qy;
qz = (-n_omega/2:(n_omega/2-1)).'.*delta_qz;

% Get spatial vector length
delta_r1 = norm(r1v);
delta_r2 = norm(r2v);
delta_r3 = norm(r3v);
delta_rx = norm(rxv);
delta_ry = norm(ryv);
delta_rz = norm(rzv);

% Get field of view
FOVx = nx.*delta_rx;
FOVy = ny.*delta_ry;
FOVz = n_omega.*delta_rz;

% Generate real space axes
x = (-nx/2:(nx/2-1)).'.*delta_rx;
y = (-ny/2:(ny/2-1)).'.*delta_ry;
z = (-n_omega/2:(n_omega/2-1)).'.*delta_rz;

% Shift amount in pixels
shft = tan(beta - pi/2)*delta_rz/delta_ry;

% Polarization factor
P = cos(tth).^2;


%% CRL calculations
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
d1 = 0.1;
d2 = fN.*(d1 + f.*phi.*tan(N.*phi))./(d1 - fN);
M = d2.*sin(N.*phi)./(f.*phi) - cos(N.*phi);

% Calculate image-detector distance, lens-detector distance, and
% sample-detector distance
d3 = abs(M)*D;
d4 = d2 + d3;
d5 = d1 + N*T + d4;

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);


%% Object
% Make simple object
w = 1.5e-3;
O1 = recta(permute(x,[2 1 3])/w).*recta(y/w).*recta(permute(z,[2 3 1])/w);

% Add phase to the object
sg = 1e-2;
s = sg*y;
u = 0.5*sg*y.^2;
phi = q.*u;
O1 = O1.*exp(1i.*phi);

% Generate object in measurement coordinates
O2 = O1;
for i = 1:n_omega
    O2(:,:,i) = circshift(O1(:,:,i),round((i - n_omega/2)*shft),1);
end


%% Generate exit field
% 3D Fourier transform the object
F = ifftshift(fftn(fftshift(O2)));

% 2D Inverse Fourier transform the slices
E1 = zeros(ny,nx,n_omega);
for i = 1:n_omega
    E1(:,:,i) = ifftshift(ifft2(fftshift(F(:,:,i))));
end

% Normalize the exit field
E1 = E1/sum(O1(ny/2+1,nx/2+1,:));

% Extract the central slice
E0 = E1(:,:,n_omega/2+1);


%% Make plots
% Make colormap
cmap = hsv(256);

% Plot the strain
figure;
plot(y,s);
title('Strain');
xlabel('y [m]');
ylabel('Strain');

% Plot the displacement field
figure;
plot(y,u);
title('Displacement');
xlabel('y [m]');
ylabel('Displacement [m]');

% Plot the phase
figure;
plot(y,phi);
title('Phase');
xlabel('y [m]');
ylabel('Phase [rad]');

% Plot the object central slice
figure;
subplot(1,3,1);
image(x,y,complex2rgb(O1(:,:,n_omega/2+1),cmap));
axis equal tight;
set(gca,'YDir','normal');
subplot(1,3,2);
imagesc(x,y,abs(O1(:,:,n_omega/2+1)));
axis equal tight;
set(gca,'YDir','normal');
subplot(1,3,3);
imagesc(x,y,angle(O1(:,:,n_omega/2+1)));
axis equal tight;
set(gca,'YDir','normal');

% Plot the exit field central slice
figure;
subplot(1,3,1);
image(x,y,complex2rgb(E0,cmap));
axis equal tight;
set(gca,'YDir','normal');
subplot(1,3,2);
imagesc(x,y,abs(E0));
axis equal tight;
set(gca,'YDir','normal');
subplot(1,3,3);
imagesc(x,y,angle(E0));
axis equal tight;
set(gca,'YDir','normal');


%% Save the data
% Rename key variables
dx = FOVx;
m = nx;

% Save the central exit field
save('Si_Exit_Field.mat','E0','m','dx');


