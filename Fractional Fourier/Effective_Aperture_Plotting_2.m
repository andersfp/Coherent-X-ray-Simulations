% Initialization
clear;
close all;
clc;


%% Load data
% Run the simulation scripts?
%run('Effective_Aperture_Small');
%run('Effective_Aperture_Large');

% Load the simulation parameters
load('Sim_Parameters.mat');

% Load the large object results
load('Effective_Aperture_Large_Profiles_2.mat');

% Save the figures?
s = 0;


%% Process data
% Coordinate
x = x2;

% Intensity
% I2 = 100*abs(L2).^2;
% IM = 100*abs(LM).^2;
I2 = abs(L2).^2;
IM = abs(LM).^2;
mm = max([max(I2) max(IM)]);
I2 = I2./mm;
IM = IM./mm;

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
plot(1e3*x,I2,'-',1e3*x,IM,'--');
title('Image intensity');
xlabel('Position [mm]');
ylabel('Intensity [a.u.]');
legend('Effective aperture','Lens by lens');

% Plot the image phase for large object
figure;
plot(1e3*x,P2,'-',1e3*x,PM,'--');
title('Image phase');
xlabel('Position [mm]');
ylabel('Phase [rad]');
legend('Effective aperture','Lens by lens');

% Plot intensity and intensity difference
figure;
subplot(2,2,1);
plot(1e3*x,I2,'-',1e3*x,IM,'--');
title('Image intensity');
xlabel('Position [mm]');
ylabel('Intensity [a.u.]');
legend('Effective aperture','Lens by lens');
subplot(2,2,2);
yyaxis left;
plot(1e3*x,Ida);
title('Intensity error');
xlabel('Position [mm]');
ylabel('Absolute error [a.u.]');
yyaxis right;
plot(1e3*x,100*Idr);
ylabel('Relative error [%]');
subplot(2,2,3);
plot(1e3*x,P2,'-',1e3*x,PM,'--');
title('Image phase');
xlabel('Position [mm]');
ylabel('Phase [rad]');
legend('Effective aperture','Lens by lens');
subplot(2,2,4);
yyaxis left;
plot(1e3*x,Pda);
title('Phase error');
xlabel('Position [mm]');
ylabel('Absolute error [rad]');
yyaxis right;
plot(1e3*x,100*Pdr);
ylabel('Relative error [%]');


%% Plots for paper
% Plotting parameters
F = 14;
lw = 2;

% Plot the objects
figure;
set(gcf,'Position',[200 100 1400 1000]);
subplot(2,2,1);
plot(1e3*x,I2,'-',1e3*x,IM,'--','LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0 0.5 0.5 0.5],'XLim',[-5 5],'YLim',[0 1],'XTick',-5:5);
xlabel('Position [mm]');
ylabel('Intensity [a.u.]');
text(-4.5,0.9,'(a)','FontSize',22,'Color','black');
legend('Effective aperture','Lens by lens');
subplot(2,2,2);
yyaxis left;
plot(1e3*x,Ida,'LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0.5 0.5 0.5 0.5],'XLim',[-5 5],'YLim',[-1.5e-4 1.5e-4],'XTick',-5:5);
xlabel('Position [mm]');
ylabel('Absolute error [a.u.]');
text(-4.5,1.2e-4,'(b)','FontSize',22,'Color','black');
yyaxis right;
plot(1e3*x,100*Idr,'LineWidth',lw);
set(gca,'YLim',[-0.15 0.15]);
ylabel('Relative error [%]');
subplot(2,2,3);
plot(1e3*x,P2,'-',1e3*x,PM,'--','LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0 0 0.5 0.5],'XLim',[-5 5],'YLim',[-pi 6*pi],'XTick',-5:5);
xlabel('Position [mm]');
ylabel('Phase [rad]');
text(-4.5,2.97e5,'(c)','FontSize',22,'Color','black');
legend('Effective aperture','Lens by lens');
subplot(2,2,4);
% yyaxis left;
plot(1e3*x,Pda,'LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0.5 0 0.5 0.5],'XLim',[-5 5],'YLim',[-2e-3 2e-3],'XTick',-5:5);
xlabel('Position [mm]');
ylabel('Absolute error [rad]');
text(-4.5,0.175,'(d)','FontSize',22,'Color','black');
% yyaxis right;
% plot(1e3*x,100*Pdr,'LineWidth',lw);
% set(gca,'YLim',[-1e-3 4e-3]);
% ylabel('Relative error [%]');
if s == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('Figure_3.pdf','-dpdf','-r0');
    print('Figure_3.png','-dpng','-r600');
end



