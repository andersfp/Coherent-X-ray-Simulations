% Initialization
clear;
close all;
clc;


%% Load data
% Load the .mat file
load('Pink_Beam_PSF.mat');

% Find the 1.3% FWHM spectrum
w0 = 1.3e-2/(2*sqrt(2*log(2)));
[~,ii] = min(abs(w - w0));

% Retrieve the 1.3% bandwidth data
p = psf(:,ii);


%% Load the measurement
% File name
f1 = 'Falch_Pink_Resolution.txt';

% Load file
dat1 = dlmread(f1,'\t');

% Extract the data
[x1,ii,~] = unique(dat1(:,1));
y1 = dat1(ii,2);

% File name
f2 = 'Falch_Pink_Source.txt';

% Load file
dat2 = dlmread(f2,'\t');

% Extract the data
[x2,ii,~] = unique(dat2(:,1));
y2 = dat2(ii,2);


%% Generate Gaussian profile
% Find FWHM
[~,ii] = min(abs(p - 0.5));
fwhm = 2*abs(x(ii));

% Generate Gaussian with same FWHM
sig1 = fwhm./(2*sqrt(2*log(2)));
g1 = exp(-x.^2./(2*sig1.^2));

% Find Gaussian width of pink spectrum
sig2 = sqrt(sum(p.*x.^2)./sum(p));
g2 = exp(-x.^2./(2*sig2.^2));

% Plot the two profiles
figure;
plot(x,[p g1 g2]);
set(gca,'XLim',[-3000 3000]);
xlabel('Position [nm]');
ylabel('PSF');
legend('Pink','Gauss (FWHM)','Gauss (sigma)');


%% Generate sample images
% Make sample
s = rect(x/8e3);

% Convolute with PSFs
sp = ifftshift(ifft(fft(fftshift(s)).*fft(fftshift(p))));
sg1 = ifftshift(ifft(fft(fftshift(s)).*fft(fftshift(g1))));
sg2 = ifftshift(ifft(fft(fftshift(s)).*fft(fftshift(g2))));

% Normalize
sp = sp./max(sp);
sg1 = sg1./max(sg1);
sg2 = sg2./max(sg2);

% Plot the sample and images
figure;
plot(x,[s sp sg1 sg2]);
set(gca,'YLim',[0 1.1]);
xlabel('Position [nm]');
ylabel('Intensity [a.u.]');
legend('Sample','Pink','Gauss (FWHM)','Gauss (sigma)');

% Plot the sample and measurement
figure;
plot(x1*1e3,y1,x+12300,sp*0.16+0.84);

% Fit the image with erf
fun = @(b,sig,x) (0.5*erfc((x - b)./(sqrt(2).*sig))).*(1 - 0.5*erfc((x + b)./(sqrt(2).*sig)));
ff = fit(x,sp,fun,'StartPoint',[4e3 800]);
figure;
plot(ff,x,sp);

% Plot the energy distributions
E0 = 17.15;
figure;
plot(x2,y2,x2,max(y2)*exp(-(x2 - E0).^2./(2.*(w0.*E0).^2)));


