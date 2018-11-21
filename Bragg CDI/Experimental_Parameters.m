% Initialization
clear;
close all;
clc;


%% Sample
% Saving path
p = 'C:\Users\anfils\Documents\Simulation_Results\BCDI_3D_HE\';

% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Lattice constant for Pt
a = 3.9242e-10;

% Momentum transfer for (111) reflektion
q = 2*pi/a*sqrt(3);

% Bragg angle
th_bragg = asin(q.*lambda./(4.*pi));


%% Standard BCDI setup
% Detector
det_dx = 75e-6;
nx = 1030; % 1030
ny = 1064; % 1064

% Sample detector distance
D = 4.0; % 2.73

% Omega scan
tth = 2*th_bragg;
omega_cen = th_bragg;
delta_omega = 0.003*pi/180;
n_omega = 500; % 500

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


%% Plot the measurement space
% Calculate ki and kf
ki = 2*pi/lambda*[-sin(th_bragg), 0, cos(th_bragg)];
kf = 2*pi/lambda*[sin(th_bragg), 0, cos(th_bragg)];

% Plot vectors
figure;
quiver3(-ki(2)*1e-9,-ki(1)*1e-9,-ki(3)*1e-9,ki(2)*1e-9,ki(1)*1e-9,ki(3)*1e-9,0);
hold on;
quiver3(0,0,0,kf(2)*1e-9,kf(1)*1e-9,kf(3)*1e-9,0);
quiver3(0,0,0,r1v(2)*1e9,r1v(1)*1e9,r1v(3)*1e9,0);
quiver3(0,0,0,r2v(2)*1e9,r2v(1)*1e9,r2v(3)*1e9,0);
quiver3(0,0,0,r3v(2)*1e9,r3v(1)*1e9,r3v(3)*1e9,0);
quiver3(kf(2)*1e-9,kf(1)*1e-9,kf(3)*1e-9,q1v(2)*1e-5,q1v(1)*1e-5,q1v(3)*1e-5);
quiver3(kf(2)*1e-9,kf(1)*1e-9,kf(3)*1e-9,q2v(2)*1e-5,q2v(1)*1e-5,q2v(3)*1e-5);
quiver3(kf(2)*1e-9,kf(1)*1e-9,kf(3)*1e-9,q3v(2)*1e-5,q3v(1)*1e-5,q3v(3)*1e-5);
axis equal tight;


%% CRL calculations
% Get material properties
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6;
T = 1.6e-3;
Tweb = 2e-6;
N = 52; % 21

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

% Save the data
save('Exp_Param.mat','-v7.3');


