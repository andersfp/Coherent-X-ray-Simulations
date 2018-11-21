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


%% Calculate MTF
% Calculate inverse space axes
nx = length(x);
qx = ((-nx/2):(nx/2 - 1)).'/mean(diff(x))/nx*2;
ny = length(y);
qy = ((-ny/2):(ny/2 - 1)).'/mean(diff(y))/ny*2;
nz = length(z);
qz = ((-nz/2):(nz/2 - 1)).'/mean(diff(z))/nz*2;

% Demonstrate difference between FFT and analytical Fourier transform
a = rect(x/1000);
b = abs(fftshift(fft(a)));
c = 167*abs(sinc(502.5*qx));
figure;
plot(qx,[b c]);

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

% Fit the frequency of the sinc
% fun = @(a,b,x) a.*abs(sinc(b.*x));
% f = fit(qx,fx1,fun,'StartPoint',[500/3 500]);

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

% Calculate MTFs
efox = envelope(fox,1,'peak');
efoy = envelope(foy,1,'peak');
efoz = envelope(foz,1,'peak');
mx1 = envelope(fx1./fx1(nx/2 + 1),1,'peak')./efox;
mxa1 = envelope(fxa1./fxa1(nx/2 + 1),1,'peak')./efox;
mx2 = envelope(fx2./fx2(nx/2 + 1),1,'peak')./efox;
mxa2 = envelope(fxa2./fxa2(nx/2 + 1),1,'peak')./efox;
my1 = envelope(fy1./fy1(ny/2 + 1),1,'peak')./efoy;
mya1 = envelope(fya1./fya1(ny/2 + 1),1,'peak')./efoy;
my2 = envelope(fy2./fy2(ny/2 + 1),1,'peak')./efoy;
mya2 = envelope(fya2./fya2(ny/2 + 1),1,'peak')./efoy;
mz1 = envelope(fz1./fz1(nz/2 + 1),1,'peak')./efoz;
mza1 = envelope(fza1./fza1(nz/2 + 1),1,'peak')./efoz;
mz2 = envelope(fz2./fz2(nz/2 + 1),1,'peak')./efoz;
mza2 = envelope(fza2./fza2(nz/2 + 1),1,'peak')./efoz;
% mx1 = fx1./fx1(nx/2 + 1)./fox;
% mxa1 = fxa1./fxa1(nx/2 + 1)./fox;
% mx2 = fx2./fx2(nx/2 + 1)./fox;
% mxa2 = fxa2./fxa2(nx/2 + 1)./fox;
% my1 = fy1./fy1(ny/2 + 1)./foy;
% mya1 = fya1./fya1(ny/2 + 1)./foy;
% my2 = fy2./fy2(ny/2 + 1)./foy;
% mya2 = fya2./fya2(ny/2 + 1)./foy;
% mz1 = fz1./fz1(nz/2 + 1)./foz;
% mza1 = fza1./fza1(nz/2 + 1)./foz;
% mz2 = fz2./fz2(nz/2 + 1)./foz;
% mza2 = fza2./fza2(nz/2 + 1)./foz;
% mx1 = fx1./fox;
% mxa1 = fxa1./fox;
% mx2 = fx2./fox;
% mxa2 = fxa2./fox;
% my1 = fy1./foy;
% mya1 = fya1./foy;
% my2 = fy2./foy;
% mya2 = fya2./foy;
% mz1 = fz1./foz;
% mza1 = fza1./foz;
% mz2 = fz2./foz;
% mza2 = fza2./foz;
% mx1 = fx1./fx1(nx/2 + 1)./abs(sinc(f.b*qx));
% mxa1 = fxa1./fxa1(nx/2 + 1)./abs(sinc(f.b*qx));
% mx2 = fx2./fx2(nx/2 + 1)./abs(sinc(f.b*qx));
% mxa2 = fxa2./fxa2(nx/2 + 1)./abs(sinc(f.b*qx));
% my1 = fy1./fy1(ny/2 + 1)./abs(sinc(f.b*qy));
% mya1 = fya1./fya1(ny/2 + 1)./abs(sinc(f.b*qy));
% my2 = fy2./fy2(ny/2 + 1)./abs(sinc(f.b*qy));
% mya2 = fya2./fya2(ny/2 + 1)./abs(sinc(f.b*qy));
% mz1 = fz1./fz1(nz/2 + 1)./abs(sinc(f.b*qz));
% mza1 = fza1./fza1(nz/2 + 1)./abs(sinc(f.b*qz));
% mz2 = fz2./fz2(nz/2 + 1)./abs(sinc(f.b*qz));
% mza2 = fza2./fza2(nz/2 + 1)./abs(sinc(f.b*qz));

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
s = 0.2;
m = 'rloess';
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


% funx = @(x) sum((fox.*(x(1)*cos(qx/x(2)*2*pi).^2 + (1 - x(1))) - fxa1./fxa1(nx/2 + 1)).^2);
% tic;
% k = fminsearch(funx,[0.8 0.32]);
% toc;

