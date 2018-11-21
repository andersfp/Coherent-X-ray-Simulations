% Initialization
clear;
close all;
clc;


%% Load data
% Load the line profiles
load('Lens_Line_Profiles.mat');

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


%% Calculate MTF
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
c = 167*abs(sinc(1005*qx));
d = c.*exp(63.18*qx.^2);
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
efx1 = envelope(fx1,1,'peak');
efxa1 = envelope(fxa1,1,'peak');
efx2 = envelope(fx2,1,'peak');
efxa2 = envelope(fxa2,1,'peak');
efy1 = envelope(fy1,1,'peak');
efya1 = envelope(fya1,1,'peak');
efy2 = envelope(fy2,1,'peak');
efya2 = envelope(fya2,1,'peak');
efz1 = envelope(fz1,1,'peak');
efza1 = envelope(fza1,1,'peak');
efz2 = envelope(fz2,1,'peak');
efza2 = envelope(fza2,1,'peak');

% Calculate MTFs
mx1 = efx1./efox;
mxa1 = efxa1./efox;
mx2 = efx2./efox;
mxa2 = efxa2./efox;
my1 = efy1./efoy;
mya1 = efya1./efoy;
my2 = efy2./efoy;
mya2 = efya2./efoy;
mz1 = efz1./efoz;
mza1 = efza1./efoz;
mz2 = efz2./efoz;
mza2 = efza2./efoz;

% Plot the x-direction MTFs
figure;
plot(qx,[mx1 mxa1 mx2 mxa2]);
title('MTF x-direction');
xlabel('qx [nm^{-1}]');
ylabel('Amplitude');
legend('1000000','1000000 avg','50000','50000 avg');
set(gca,'YLim',[0 1.5]);

% Plot the y-direction MTFs
figure;
plot(qy,[my1 mya1 my2 mya2]);
title('MTF y-direction');
xlabel('qy [nm^{-1}]');
ylabel('Amplitude');
legend('1000000','1000000 avg','50000','50000 avg');
set(gca,'YLim',[0 1.5]);

% Plot the z-direction MTFs
figure;
plot(qz,[mz1 mza1 mz2 mza2]);
title('MTF z-direction');
xlabel('qz [nm^{-1}]');
ylabel('Amplitude');
legend('1000000','1000000 avg','50000','50000 avg');
set(gca,'YLim',[0 1.5]);

% Smooth MTFs
s = 0.03;
m = 'lowess';
mx1s = smooth(mx1,s,m);
mxa1s = smooth(mxa1,s,m);
mx2s = smooth(mx2,s,m);
mxa2s = smooth(mxa2,s,m);
my1s = smooth(my1,s,m);
mya1s = smooth(mya1,s,m);
my2s = smooth(my2,s,m);
mya2s = smooth(mya2,s,m);
mz1s = smooth(mz1,s,m);
mza1s = smooth(mza1,s,m);
mz2s = smooth(mz2,s,m);
mza2s = smooth(mza2,s,m);

% Plot the x-direction MTFs
figure;
plot(qx,[mx1s mxa1s mx2s mxa2s]);
title('Smoothed MTF x-direction');
xlabel('qx [nm^{-1}]');
ylabel('Amplitude');
legend('1000000','1000000 avg','50000','50000 avg');
set(gca,'YLim',[0 1.1]);

% Plot the y-direction MTFs
figure;
plot(qy,[my1s mya1s my2s mya2s]);
title('Smoothed MTF y-direction');
xlabel('qy [nm^{-1}]');
ylabel('Amplitude');
legend('1000000','1000000 avg','50000','50000 avg');
set(gca,'YLim',[0 1.1]);

% Plot the z-direction MTFs
figure;
plot(qz,[mz1s mza1s mz2s mza2s]);
title('Smoothed MTF z-direction');
xlabel('qz [nm^{-1}]');
ylabel('Amplitude');
legend('1000000','1000000 avg','50000','50000 avg');
set(gca,'YLim',[0 1.1]);






