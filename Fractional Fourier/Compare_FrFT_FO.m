% Initialization
clear;
close all;
clc;


%% Set up the object
% Set wavelength
lambda = 1e-10;
k = 2*pi/lambda;

% Set propagation distance
d = 1;

% Set the number of pizels
m = 1024;

% Set the side length
L = sqrt(m*lambda*d);

% Make the axis
x0 = ((-m/2):(m/2-1)).'.*L/m;

% Make a square
w = L/4;
E0 = recta(x0/w);

% Plot the object
figure;
plot(x0,E0);


%% FrFT propagation
% FrFT parameters
R0 = Inf;
s0 = L/sqrt(m);
[a,R0,R1,s0,s1,g0,g1] = FrFT_parameters(d,[],lambda,R0,s0);

% FrFT propagation
[E1,x1] = propFrFT1(E0,x0,R0,R1,s0,s1,a,lambda,d);

% Plot the propagation
figure;
plotyy(x1,abs(E1),x1,unwrap(angle(E1)));


%% Fourier optics propagation
% Propagation
x2 = x0;
E2 = ifftshift(fft(fftshift(E0.*exp(1i*k/(2*d)*x0.^2)))).*exp(1i*k*d).*exp(-1i*pi/4).*exp(1i*k/(2*d)*x2.^2)/sqrt(m);

% Plot the propagation
figure;
plotyy(x2,abs(E2),x2,unwrap(angle(E2)));


%% Compare the results
% Get the intensities and phases
I1 = abs(E1).^2;
I2 = abs(E2).^2;
P1 = angle(E1);
P2 = angle(E2);

% Get unwrapped phases
U1 = unwrap(angle(E1.*exp(-1i.*pi.*x1.^2./(lambda.*R1)))) + pi.*x1.^2./(lambda.*R1);
U2 = unwrap(angle(E2));
U1 = U1 - min(U1);
U2 = U2 - min(U2);

% Plot the intensities
figure;
plot(x1,I1,x2,I2);

% Plot the phases
figure;
plot(x1,P1,x2,P2);

% Plot the unwrapped phases
figure;
plot(x1,U1,x2,U2);


