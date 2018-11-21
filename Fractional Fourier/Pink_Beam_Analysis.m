% Initialization
clear;
close all;
clc;


%% Load the simulation results
% Load the .mat file
load('Pink_Beam_Profiles.mat');

% Save the plots?
s = 0;


%% Process the data
% Generate intensity
I = abs(L).^2;

% Plot the raw spectra intensity
figure;
pcolor(1e-3*Es,1e6*x,I/max(max(I))*1e3);
shading flat;
xlabel('Energy (keV)');
ylabel('Position (\mum)');

% Calculate the RMS width
sig = sqrt(sum(x.^2.*I)./sum(I));

% Plot the width of each PSF
figure;
plot(1e-3*Es,1e6*sig);
xlabel('Energy (keV)');
ylabel('PSF width (\mum)');

% Generate new axis
m2 = 20000;
dx2 = 5e-4;
x2 = ((-m2/2):(m2/2 - 1)).'/m2*dx2;

% Generate new dataset
fun = @(sig,x) exp(-x.^2./(2*(sig).^2));
I2 = I(m/2+1,:).*fun(sig,x2);

% Plot the new dataset
figure;
imagesc(1e-3*Es,1e6*x2,I2/max(max(I2))*1e3);
xlabel('Energy (keV)');
ylabel('Position (\mum)');

% Generate the bandwidth array
%w = [1e-10 1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3 8e-3 9e-3 1e-2];
w = [1e-10 (1:100)*1e-4];
w = permute(w,[3 1 2]);
dE = w.*E;
n = length(w);

% Generate weigthing factor
S = fun(dE,Es - E);

% Plot the weighting factor
figure;
imagesc(squeeze(w),1e-3*Es,squeeze(S));
xlabel('Relative bandwidth');
ylabel('Energy (keV)');

% Generate total intensities
It = squeeze(sum(I2.*S,2));

% Plot the total intensities
figure;
imagesc(squeeze(w),1e6*x2,It./max(max(It))*1e3);
xlabel('Relative bandwidth');
ylabel('Position (\mum)');
set(gca,'YLim',[-5 5]);

% Plot the total normalized intensities
figure;
imagesc(squeeze(w),1e6*x2,It./max(It));
xlabel('Relative bandwidth');
ylabel('Position (\mum)');
set(gca,'YLim',[-5 5]);

% Calculate the RMS width
sigt = sqrt(sum(x2.^2.*It)./sum(It));

% Plot the RMS width
figure;
plot(1e2*squeeze(w),1e6*sigt);
xlabel('Relative bandwidth (%)');
ylabel('RMS width (\mum)');


%% Geometrical optics
% Calculate parameters
i0 = (nE + 1)/2;
dch = N*phi(i0)*(d0*dd/(f(i0)*phi(i0)) - f(i0)*phi(i0))*cos(N*phi(i0)) + (d0*dd/(f(i0)*phi(i0)) + f(i0)*phi(i0) + N*phi(i0)*(d0 + dd))*sin(N*phi(i0));

% Calculate the chromatic intensity distribution
Ig = squeeze(exp(-abs(x2)./(sigA(i0).*w.*dch)));

% Convolute the intensity distribution with aperture PSF
PSFa = exp(-x2.^2./(2*(sigt(1)).^2));
PSFa = PSFa/sum(PSFa);
for i = 1:n
    Ig(:,i) = fconv1(Ig(:,i),PSFa);
end

% Plot the total intensities
figure;
imagesc(squeeze(w),1e6*x2,Ig./max(max(Ig))*1e3);
xlabel('Relative bandwidth');
ylabel('Position (\mum)');
set(gca,'YLim',[-5 5]);

% Plot the normalized intensities
figure;
imagesc(squeeze(w),1e6*x2,Ig./max(Ig));
xlabel('Relative bandwidth');
ylabel('Position (\mum)');
set(gca,'YLim',[-5 5]);

% Plots
lm = 10;
%c1 = winter(n);
%c2 = summer(n);
c1 = [zeros(n,2) linspace(0,1,n).'];
c2 = [linspace(1,0,n).' zeros(n,2)];

figure;
plotyy(squeeze(w),sum(It),squeeze(w),sum(Ig));
xlabel('Relative bandwidth');

figure;
hold on;
for i = 1:n
    plot(1e6*x2,It(:,i)./It(m2/2+1,i) + (i - 1)*0.2,'Color',c1(i,:));
end
hold off;
xlabel('Position (\mum)');
ylabel('Intensity');
set(gca,'XLim',[-lm lm]);

figure;
hold on;
for i = 1:n
    plot(1e6*x2,Ig(:,i)./Ig(m2/2+1,i) + (i - 1)*0.1,'Color',c2(i,:));
end
hold off;
xlabel('Position (\mum)');
ylabel('Intensity');
set(gca,'XLim',[-lm lm]);

figure;
hold on;
for i = 1:n
    plot(1e6*x2,It(:,i)./It(m2/2+1,i) + (i - 1)*0.5,'Color',c1(i,:));
    plot(1e6*x2,Ig(:,i)./Ig(m2/2+1,i) + (i - 1)*0.5,'Color',c2(i,:));
end
hold off;
xlabel('Position (\mum)');
ylabel('Intensity');
set(gca,'XLim',[-lm lm]);

% Plot the intensity profiles vs bandwidth
fs = 14;
figure;
set(gcf,'Position',[150 370 1200 500]);
subplot(1,2,1);
imagesc(100*squeeze(w),1e9*x2/M,It./max(max(It))*1e3);
hold on;
[C,h] = contour(100*squeeze(w),1e9*x2/M,It./max(It),[0.25 0.5 0.75],'Color','black');
clabel(C,h,'FontSize',12,'LabelSpacing',170);
hold off;
xlabel('Relative bandwidth [%]');
ylabel('Position [nm]');
set(gca,'YLim',[-500 500],'YTick',-400:200:400,'YDir','normal','FontSize',fs,'OuterPosition',[0 0 0.5 1]);
text(5e-4,450,'(a)','FontSize',22,'Color','white');
subplot(1,2,2);
imagesc(100*squeeze(w),1e9*x2/M,Ig./max(max(Ig))*1e3);
hold on;
[C,h] = contour(100*squeeze(w),1e9*x2/M,Ig./max(Ig),[0.25 0.5 0.75],'Color','black');
clabel(C,h,'FontSize',12,'LabelSpacing',110);
hold off;
xlabel('Relative bandwidth [%]');
ylabel('Position [nm]');
set(gca,'YLim',[-500 500],'YTick',-400:200:400,'YDir','normal','FontSize',fs,'OuterPosition',[0.5 0 0.5 1]);
text(5e-4,450,'(b)','FontSize',22,'Color','white');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Pink_Beam.pdf','-dpdf','-r0');
    print('Pink_Beam.png','-dpng','-r600');
end

% Plot the normalized intensity profiles vs bandwidth
figure;
set(gcf,'Position',[150 370 1200 500]);
subplot(1,2,1);
imagesc(100*squeeze(w),1e9*x2/M,It./max(It));
hold on;
[C,h] = contour(100*squeeze(w),1e9*x2/M,It./max(It),[0.25 0.5 0.75],'Color','black');
clabel(C,h,'FontSize',12,'LabelSpacing',170);
hold off;
xlabel('Relative bandwidth [%]');
ylabel('Position [nm]');
set(gca,'YLim',[-500 500],'YTick',-400:200:400,'YDir','normal','FontSize',fs,'OuterPosition',[0 0 0.5 1]);
text(5e-4,450,'(a)','FontSize',22,'Color','white');
subplot(1,2,2);
imagesc(100*squeeze(w),1e9*x2/M,Ig./max(Ig));
hold on;
[C,h] = contour(100*squeeze(w),1e9*x2/M,Ig./max(Ig),[0.25 0.5 0.75],'Color','black');
clabel(C,h,'FontSize',12,'LabelSpacing',110);
hold off;
xlabel('Relative bandwidth [%]');
ylabel('Position [nm]');
set(gca,'YLim',[-500 500],'YTick',-400:200:400,'YDir','normal','FontSize',fs,'OuterPosition',[0.5 0 0.5 1]);
text(5e-4,450,'(b)','FontSize',22,'Color','white');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Pink_Beam_Normalized.pdf','-dpdf','-r0');
    print('Pink_Beam_Normalized.png','-dpng','-r600');
end

% Plot selected PSFs
ii = [1 2 11 101];
figure;
plot(1e9*x2/M,It(:,ii)./It(m2/2+1,ii),'LineWidth',2);
xlabel('Position [nm]');
ylabel('Normalized intensity');
legend('\DeltaE/E = 0','\DeltaE/E = 10^{-4}','\DeltaE/E = 10^{-3}','\DeltaE/E = 10^{-2}');
set(gca,'XLim',[-300 300],'FontSize',fs);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Pink_Beam_Profiles.pdf','-dpdf','-r0');
    print('Pink_Beam_Profiles.png','-dpng','-r600');
end

% sse = zeros(nE,1);
% rsquare = sse;
% dfe = sse;
% adjrsquare = sse;
% rmse = sse;
% sig2 = sse;
% for i = 1:nE
%     [ff,gg] = fit(x(:,i),I(:,i)./max(I(:,i)),fun,'StartPoint',sig(i));
%     figure;plot(ff,x(:,i),I(:,i)./max(I(:,i)));
%     sse(i) = gg.sse;
%     rsquare(i) = gg.rsquare;
%     dfe(i) = gg.dfe;
%     adjrsquare(i) = gg.adjrsquare;
%     rmse(i) = gg.rmse;
%     sig2(i) = ff.sig;
% end
% figure;
% plot(Es,sse);
% figure;
% plot(Es,rsquare);
% figure;
% plot(Es,dfe);
% figure;
% plot(Es,adjrsquare);
% figure;
% plot(Es,rmse);
% figure;
% plot(Es,sig,Es,sig2);


