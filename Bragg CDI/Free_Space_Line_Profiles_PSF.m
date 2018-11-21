% Initialization
clear;
close all;
clc;


%% Load data
% Load the line profiles
load('Free_Space_Line_Profiles.mat');

% Convert the axes to nm
x = 1e9*x;
y = 1e9*y;
z = 1e9*z;

% Normalize the amplitude
lx1 = lx1/mean(lx1(lx1 > 0));
lxa1 = lxa1/mean(lxa1(lxa1 > 0));
lx2 = lx2/mean(lx2(lx2 > 0));
lxa2 = lxa2/mean(lxa2(lxa2 > 0));
ly1 = ly1/mean(ly1(ly1 > 0));
lya1 = lya1/mean(lya1(lya1 > 0));
ly2 = ly2/mean(ly2(ly2 > 0));
lya2 = lya2/mean(lya2(lya2 > 0));
lz1 = lz1/mean(lz1(lz1 > 0));
lza1 = lza1/mean(lza1(lza1 > 0));
lz2 = lz2/mean(lz2(lz2 > 0));
lza2 = lza2/mean(lza2(lza2 > 0));

% Shift lines to center the lines
lx1 = circshift(lx1,-round(sum(lx1.*x)./sum(lx1)./mean(diff(x))));
lxa1 = circshift(lxa1,-round(sum(lxa1.*x)./sum(lxa1)./mean(diff(x))));
lx2 = circshift(lx2,-round(sum(lx2.*x)./sum(lx2)./mean(diff(x))));
lxa2 = circshift(lxa2,-round(sum(lxa2.*x)./sum(lxa2)./mean(diff(x))));
ly1 = circshift(ly1,-round(sum(ly1.*y)./sum(ly1)./mean(diff(y))));
lya1 = circshift(lya1,-round(sum(lya1.*y)./sum(lya1)./mean(diff(y))));
ly2 = circshift(ly2,-round(sum(ly2.*y)./sum(ly2)./mean(diff(y))));
lya2 = circshift(lya2,-round(sum(lya2.*y)./sum(lya2)./mean(diff(y))));
lz1 = circshift(lz1,-round(sum(lz1.*z)./sum(lz1)./mean(diff(z))));
lza1 = circshift(lza1,-round(sum(lza1.*z)./sum(lza1)./mean(diff(z))));
lz2 = circshift(lz2,-round(sum(lz2.*z)./sum(lz2)./mean(diff(z))));
lza2 = circshift(lza2,-round(sum(lza2.*z)./sum(lza2)./mean(diff(z))));


%% Plot the data
% Plot the x-lines
figure;
plot(x,[lx1 lxa1 lx2 lxa2]);
title('x-lines');
xlabel('x [nm]');
ylabel('Normalized amplitude');
legend('1000000','1000000 avg','50000','50000 avg');

% Plot the y-lines
figure;
plot(y,[ly1 lya1 ly2 lya2]);
title('y-lines');
xlabel('y [nm]');
ylabel('Normalized amplitude');
legend('1000000','1000000 avg','50000','50000 avg');

% Plot the z-lines
figure;
plot(z,[lz1 lza1 lz2 lza2]);
title('z-lines');
xlabel('z [nm]');
ylabel('Normalized amplitude');
legend('1000000','1000000 avg','50000','50000 avg');


%% Calculate PSF
% Calculate inverse space axes
nx = length(x);
qx = ((-nx/2):(nx/2 - 1)).'/mean(diff(x))/nx;
ny = length(y);
qy = ((-ny/2):(ny/2 - 1)).'/mean(diff(y))/ny;
nz = length(z);
qz = ((-nz/2):(nz/2 - 1)).'/mean(diff(z))/nz;

% Demonstrate difference between FFT and analytical Fourier transform
a = rect(x/1000);
b = abs(fftshift(fft(a)));
c = 265*abs(sinc(1000*qx));
d = c.*exp(22.31*qx.^2);
figure;
plot(qx,[b c d]);

% Calculate Fourier transforms
fx1 = abs(fftshift(fft(lx1)));
fxa1 = abs(fftshift(fft(lxa1)));
fx2 = abs(fftshift(fft(lx2)));
fxa2 = abs(fftshift(fft(lxa2)));
fy1 = abs(fftshift(fft(ly1)));
fya1 = abs(fftshift(fft(lya1)));
fy2 = abs(fftshift(fft(ly2)));
fya2 = abs(fftshift(fft(lya2)));
fz1 = abs(fftshift(fft(lz1)));
fza1 = abs(fftshift(fft(lza1)));
fz2 = abs(fftshift(fft(lz2)));
fza2 = abs(fftshift(fft(lza2)));

% Normalize the Fourier transforms
fx1 = fx1./fx1(nx/2 + 1);
fxa1 = fxa1./fxa1(nx/2 + 1);
fx2 = fx2./fx2(nx/2 + 1);
fxa2 = fxa2./fxa2(nx/2 + 1);
fy1 = fy1./fy1(ny/2 + 1);
fya1 = fya1./fya1(ny/2 + 1);
fy2 = fy2./fy2(ny/2 + 1);
fya2 = fya2./fya2(ny/2 + 1);
fz1 = fz1./fz1(nz/2 + 1);
fza1 = fza1./fza1(nz/2 + 1);
fz2 = fz2./fz2(nz/2 + 1);
fza2 = fza2./fza2(nz/2 + 1);

% Generate an ideal object fft
ox = rect(x/1000);
oy = rect(y/1000);
oz = rect(z/1000);
fox = abs(fftshift(fft(ox)));
foy = abs(fftshift(fft(oy)));
foz = abs(fftshift(fft(oz)));
fox = fox./fox(nx/2 + 1);
foy = foy./foy(ny/2 + 1);
foz = foz./foz(nz/2 + 1);

% Calculate evelopes
efox = envelope(fox,1,'peak');
efoy = envelope(foy,1,'peak');
efoz = envelope(foz,1,'peak');

% Fit the PSF in Fourier space
sig0 = 6;
w0 = 1000;

vx1 = PSF_Fitting(x,fx1,1./efox.*(abs(qx) < 0.036),2*550,sig0);
figure;
semilogy(qx,[fx1 MTF_Fitting(x,[w0;2*550;vx1])]);

vxa1 = PSF_Fitting(x,fxa1,1./efox.*(abs(qx) < 0.035),2*550,sig0);
figure;
semilogy(qx,[fxa1 MTF_Fitting(x,[w0;2*550;vxa1])]);

vx2 = PSF_Fitting(x,fx2,1./efox.*(abs(qx) < 0.055),2*550,sig0);
figure;
semilogy(qx,[fx2 MTF_Fitting(x,[w0;2*550;vx2])]);

vxa2 = PSF_Fitting(x,fxa2,1./efox.*(abs(qx) < 0.048),2*550,sig0);
figure;
semilogy(qx,[fxa2 MTF_Fitting(x,[w0;2*550;vxa2])]);

vx = [vx1 vxa1 vx2 vxa2];

vy1 = PSF_Fitting(y,fy1,1./efoy.*(abs(qy) < 0.026),2*550,sig0);
figure;
semilogy(qy,[fy1 MTF_Fitting(y,[w0;2*550;vy1])]);

vya1 = PSF_Fitting(y,fya1,1./efoy.*(abs(qy) < 0.043),2*550,sig0);
figure;
semilogy(qy,[fya1 MTF_Fitting(y,[w0;2*550;vya1])]);

vy2 = PSF_Fitting(y,fy2,1./efoy.*(abs(qy) < 0.026),2*550,sig0);
figure;
semilogy(qy,[fy2 MTF_Fitting(y,[w0;2*550;vy2])]);

vya2 = PSF_Fitting(y,fya2,1./efoy.*(abs(qy) < 0.029),2*550,sig0);
figure;
semilogy(qy,[fya2 MTF_Fitting(y,[w0;2*550;vya2])]);

vy = [vy1 vya1 vy2 vya2];

vz1 = PSF_Fitting(z,fz1,1./efoz.*(abs(qz) < 0.020),2*550,14);
figure;
semilogy(qz,[fz1 MTF_Fitting(z,[w0;2*550;vz1])]);

vza1 = PSF_Fitting(z,fza1,1./efoz.*(abs(qz) < 0.031),2*550,10);
figure;
semilogy(qz,[fza1 MTF_Fitting(z,[w0;2*550;vza1])]);

vz2 = PSF_Fitting(z,fz2,1./efoz.*(abs(qz) < 0.011),2*550,29);
figure;
semilogy(qz,[fz2 MTF_Fitting(z,[w0;2*550;vz2])]);

vza2 = PSF_Fitting(z,fza2,1./efoz.*(abs(qz) < 0.011),2*550,24);
figure;
semilogy(qz,[fza2 MTF_Fitting(z,[w0;2*550;vza2])]);

vz = [vz1 vza1 vz2 vza2];

v = [vx;vy;vz];

dr = [mean(diff(x));mean(diff(y));mean(diff(z))];


