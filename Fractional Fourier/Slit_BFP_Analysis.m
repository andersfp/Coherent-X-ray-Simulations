% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters and results
load('Slit_BFP_Propagation_Profiles.mat');

% Save figure?
s = 0;
pl = 0;


%% Load results
% Calibrate the slit sizes
sl = sum(rect(x4./sl.')).*mean(diff(x4)).';

% Calculate intensities
I = abs(L).^2;

% Normalize the intensity
In = I./max(I);

% Plot the profiles
figure;
plot(xd,I);

% Plot selected profiles
figure;
semilogy(xd,I(:,1:9:end));
set(gca,'YLim',[1e-15 1e-5]);

% Plot selected profiles
figure;
plot(xd,In(:,1:9:end));
set(gca,'XLim',[-2e-5 2e-5]);


%% Fit the line profiles
% Make fitting model
fun = @(sig,x) exp(-x.^2./(2*(1e-6*sig).^2));

% Make the fits
res = zeros(n,1);
rese = res;
for i = 1:n
    fr = fit(xd,In(:,i),fun,'StartPoint',1);
    ci = confint(fr,0.682689492);
    ci = diff(ci);
    res(i) = fr.sig*1e-6/M;
    rese(i) = ci(1)*1e-6/M;
    if pl == 1
        figure;
        plot(fr,xd,In(:,i));
        set(gca,'XLim',[-5*M*res(i) 5*M*res(i)]);
    end
end

% Plot the resolution
figure;
errorbar(1e3*sl,1e9*res,1e9*rese);
set(gca,'XScale','log','YScale','log');
title('On-axis resolution vs. slit size');
xlabel('Slit size (mm)');
ylabel('Resolution (nm)');


%% Compare to geometrical optics results
% Calculate CRL resolution based on critical angle
res_CRL = lambda/(4*pi*sigA);

% Calculate slit resolution
res_sl = 0.3645*lambda*fN./(sl*cos(N*phi));

% Calculate total resolution
res_tot = sqrt(res_CRL.^2 + res_sl.^2);

% Plot fitted resolution with predicted resolutions
fs = 14;
lw = 2;
figure;
plot(1e3*[sl(1) sl(end)],1e9*res_CRL*[1 1],'LineWidth',lw);
hold on;
plot(1e3*sl,1e9*res_sl,'LineWidth',lw);
plot(1e3*sl,1e9*res_tot,'LineWidth',lw);
%errorbar(1e3*sl,1e9*res,1e9*rese,'LineWidth',lw);
plot(1e3*sl,1e9*res,'LineWidth',lw);
hold off;
set(gca,'XScale','log','YScale','log','XLim',[0.028 2.5],'YLim',[10 300]);
xlabel('Aperture size [mm]');
ylabel('Resolution [nm]');
%legend('RTM CRL','RTM slit','RTM combined','FrFT simulation');
legend('\sigma_{CRL}','\sigma_{slit}','\sigma_{tot}','\sigma_{FrFT}');
set(gca,'FontSize',fs,'XTick',[0.03 0.1 0.3 1],'YTick',[10 30 100 300]);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Slit_BFP.pdf','-dpdf','-r0');
    print('Slit_BFP.png','-dpng','-r600');
end


