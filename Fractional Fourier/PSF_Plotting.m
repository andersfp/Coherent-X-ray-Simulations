% Initialization
clear;
close all;
clc;


%% Load data
% Load the large object results
load('PSF_Propagation_Profiles.mat');

% Save the figures?
s = 0;


%% Process data
% Coordinate
x = x2;

% Intensity
I2 = abs(L2).^2;
IM = abs(LM).^2;

% Normalize the intensity
I2 = I2/max(IM);
IM = IM/max(IM);

% Phase
% P2 = unwrap(angle(L2)) + pi*x.^2/(lambda*Rp(end));
% PM = unwrap(angle(LM)) + pi*x.^2/(lambda*Rp(end));
P2 = unwrap(angle(L2));
PM = unwrap(angle(LM));
%P2 = unwrap(angle(L2));
%PM = unwrap(angle(LM));

% Set the phase to 0 at x=0
[~,i0] = min(abs(x));
P2 = P2 - P2(i0);
PM = PM - PM(i0);

% Intensity differences
Ida = I2 - IM;
Idr = Ida./IM;

% Phase differences
Pda = P2 - PM;
Pdr = Pda./PM;


%% Make plots
% Plot the image plane for large object
figure;
plot(1e6*x,I2,'-',1e6*x,IM,'--');
title('Image intensity');
xlabel('Position [\mum]');
ylabel('Intensity [a.u.]');
legend('Effective aperture','Lens by lens');

% Plot the image phase for large object
figure;
plot(1e6*x,P2,'-',1e6*x,PM,'--');
title('Image phase');
xlabel('Position [\mum]');
ylabel('Phase [rad]');
legend('Effective aperture','Lens by lens');

% Plot intensity and intensity difference
figure;
subplot(2,2,1);
plot(1e6*x,I2,'-',1e6*x,IM,'--');
title('Image intensity');
xlabel('Position [\mum]');
ylabel('Intensity [a.u.]');
legend('Effective aperture','Lens by lens');
subplot(2,2,2);
yyaxis left;
plot(1e6*x,Ida);
title('Intensity error');
xlabel('Position [\mum]');
ylabel('Absolute error [a.u.]');
yyaxis right;
plot(1e6*x,100*Idr);
ylabel('Relative error [%]');
subplot(2,2,3);
plot(1e6*x,P2,'-',1e6*x,PM,'--');
title('Image phase');
xlabel('Position [\mum]');
ylabel('Phase [rad]');
legend('Effective aperture','Lens by lens');
subplot(2,2,4);
yyaxis left;
plot(1e6*x,Pda);
title('Phase error');
xlabel('Position [\mum]');
ylabel('Absolute error [rad]');
yyaxis right;
plot(1e6*x,100*Pdr);
ylabel('Relative error [%]');


%% Plots for paper
% Plotting parameters
F = 14;
lw = 2;

% Plot the objects
figure;
set(gcf,'Position',[200 100 1400 1000]);
subplot(2,2,1);
plot(1e9*x/M,I2,'-',1e9*x/M,IM,'--','LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0 0.5 0.5 0.5],'XLim',[-150 150],'YLim',[0 1],'XTick',-150:50:150,'YTick',0:0.2:1);
xlabel('Position [nm]');
ylabel('Intensity [a.u.]');
text(-145,0.95,'(a)','FontSize',22,'Color','black');
legend('Effective aperture','Lens by lens');
subplot(2,2,2);
yyaxis left;
plot(1e9*x/M,Ida,'LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0.5 0.5 0.5 0.5],'XLim',[-150 150],'YLim',[-6e-8 6e-8],'XTick',-150:50:150,'YTick',-6e-8:2e-8:6e-8);
xlabel('Position [nm]');
ylabel('Absolute error [a.u.]');
text(-145,5.4e-8,'(b)','FontSize',22,'Color','black');
yyaxis right;
plot(1e9*x/M,100*Idr,'LineWidth',lw);
set(gca,'YLim',[-5e-4 5e-4],'YTick',-4e-4:2e-4:4e-4);
ylabel('Relative error [%]');
subplot(2,2,3);
plot(1e9*x/M,P2,'-',1e9*x/M,PM,'--','LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0 0 0.5 0.5],'XLim',[-150 150],'YLim',[-1e-3 1e-3],'XTick',-150:50:150,'YTick',-1e-3:0.5e-3:1e-3);
xlabel('Position [nm]');
ylabel('Phase [rad]');
text(-145,9e-4,'(c)','FontSize',22,'Color','black');
legend('Effective aperture','Lens by lens');
subplot(2,2,4);
% yyaxis left;
plot(1e9*x/M,Pda,'LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0.5 0 0.5 0.5],'XLim',[-150 150],'YLim',[-5e-4 5e-4],'XTick',-150:50:150,'YTick',-4e-4:2e-4:4e-4);
xlabel('Position [nm]');
ylabel('Absolute error [rad]');
text(-145,4.5e-4,'(d)','FontSize',22,'Color','black');
% yyaxis right;
% plot(1e9*x/M,100*Pdr,'LineWidth',lw);
% set(gca,'YLim',[-200 200]);
% ylabel('Relative error [%]');
if s == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('PSF.pdf','-dpdf','-r0');
    print('PSF.png','-dpng','-r600');
end



