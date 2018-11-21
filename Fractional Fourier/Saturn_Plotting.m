% Initialization
clear;
close all;
clc;


%% Load data
% Load nano
N = load('Nano_Saturn_Propagation.mat');

% Load micro
M = load('Micro_Saturn_Propagation.mat');

% Save the figures?
s = 0;


%% Plots for paper
% Figure parameters
m = max(max(N.I0));

% Plotting parameters
F = 14;

% Plot the objects
figure;
set(gcf,'Position',[200 100 1300 900]);
subplot(2,3,1);
imagesc(1e6*N.x0,1e6*N.y0,N.I0,[0 m]);
axis equal tight;
colormap gray;
set(gca,'YDir','normal','FontSize',F,'OuterPosition',[0 1/2 1/3 1/2],'XTick',-2:2,'YTick',-3:3);
xlabel('x [\mum]');
ylabel('y [\mum]');
text(-2.8,3.2,'(a)','FontSize',22,'Color','white');
subplot(2,3,2);
imagesc(1e6*N.xM,1e6*N.yM,N.IM,[0 m/100]);
axis equal tight;
colormap gray;
set(gca,'YDir','normal','FontSize',F,'OuterPosition',[1/3 1/2 1/3 1/2],'XTick',-20:10:20,'YTick',-30:10:30);
xlabel('x [\mum]');
ylabel('y [\mum]');
text(-28,32,'(b)','FontSize',22,'Color','white');
subplot(2,3,3);
imagesc(1e6*N.x2,1e6*N.y2,N.I2,[0 m/100]);
axis equal tight;
colormap gray;
set(gca,'YDir','normal','FontSize',F,'OuterPosition',[2/3 1/2 1/3 1/2],'XTick',-20:10:20,'YTick',-30:10:30);
xlabel('x [\mum]');
ylabel('y [\mum]');
text(-28,32,'(c)','FontSize',22,'Color','white');
subplot(2,3,4);
imagesc(1e3*M.x0,1e3*M.y0,M.I0,[0 m]);
axis equal tight;
colormap gray;
set(gca,'YDir','normal','FontSize',F,'OuterPosition',[0 0 1/3 1/2],'XTick',-0.2:0.1:0.2,'YTick',-0.2:0.1:0.2);
xlabel('x [mm]');
ylabel('y [mm]');
text(-0.224,0.256,'(d)','FontSize',22,'Color','white');
subplot(2,3,5);
imagesc(1e3*M.xM,1e3*M.yM,M.IM,[0 m/100]);
axis equal tight;
colormap gray;
set(gca,'YDir','normal','FontSize',F,'OuterPosition',[1/3 0 1/3 1/2],'XTick',-2:2,'YTick',-2:2);
xlabel('x [mm]');
ylabel('y [mm]');
text(-2.24,2.56,'(e)','FontSize',22,'Color','white');
subplot(2,3,6);
imagesc(1e3*M.x2,1e3*M.y2,M.I2,[0 m/100]);
axis equal tight;
colormap gray;
set(gca,'YDir','normal','FontSize',F,'OuterPosition',[2/3 0 1/3 1/2],'XTick',-2:2,'YTick',-2:2);
xlabel('x [mm]');
ylabel('y [mm]');
text(-2.24,2.56,'(f)','FontSize',22,'Color','white');

% Save the figure
if s == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('Saturn.pdf','-dpdf','-r0');
    print('Saturn.png','-dpng','-r600');
end



