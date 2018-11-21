% Initialization
clear;
close all;
clc;


%% CRL calculations
% Set energy
E = 17e3;
lambda = E2lambda(E);

% Get material properties
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6;
T = 2e-3;
Tweb = 2e-6;
N = 157;

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set object and image positions
d1 = 0.015;
d2 = fN.*(d1 + f.*phi.*tan(N.*phi))./(d1 - fN);
M = d2.*sin(N.*phi)./(f.*phi) - cos(N.*phi);

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);

%figure;
%semilogy(N,fN);

%figure;
%semilogy(d1,d1+d2+N*T);

%figure;
%semilogy(d1,M);


