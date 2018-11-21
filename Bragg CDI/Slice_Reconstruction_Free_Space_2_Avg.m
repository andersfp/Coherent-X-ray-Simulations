% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the original object
F2 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_avg_xyz.bin'],[ny nx n_omega]);

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures\';
sav = 0;


%% Plot the data
% Scale the data
F2 = 200*F2;

% Offset the phase
mp = sum(angle(F2(F2 ~= 0)))/sum(sum(sum(F2 ~= 0)));
F2 = F2.*exp(-1i*mp);

% Extract the 3 slices
xy = F2(:,:,235);
zy = squeeze(F2(:,528,:));
xz = squeeze(F2(495,:,:)).';

% Make the color slices
cmap = hsv(256);
xy = complex2rgb(xy,cmap);
zy = complex2rgb(zy,cmap);
xz = complex2rgb(xz,cmap);

% Set plotting font size
F = 16;
s = 1e6;
u = '[\mum]';
xt = -3:3;
yt = -3:3;
zt = -2:2;

% Plot the xy slice
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,xy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Reconstruction_Free_Space_2_Avg_XY.pdf'],'-dpdf','-r0');
    print([ps 'Reconstruction_Free_Space_2_Avg_XY.png'],'-dpng','-r600');
end

% Plot the zy slice
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,zy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Reconstruction_Free_Space_2_Avg_ZY.pdf'],'-dpdf','-r0');
    print([ps 'Reconstruction_Free_Space_2_Avg_ZY.png'],'-dpng','-r600');
end

% Plot the xz slice
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*z,xz);
axis equal tight;
xlabel(['x ' u]);
ylabel(['z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Reconstruction_Free_Space_2_Avg_XZ.pdf'],'-dpdf','-r0');
    print([ps 'Reconstruction_Free_Space_2_Avg_XZ.png'],'-dpng','-r600');
end

% Generate color-wheel
[X,Y] = meshgrid(linspace(-1,1,1000),linspace(-1,1,1000));
A = double(hypot(X,Y) <= 0.95).*exp(1i*atan2(Y,X));
B = complex2rgb(A,cmap);
figure;
image(B);
ax1 = gca;
axis equal tight;
ax2 = axes('Position',ax1.Position,'Color', 'none');
axis equal tight;
set(ax2,'XAxisLocation','top','YAxisLocation','Right');
set(ax2,'XLim',ax1.XLim,'YLim',ax1.YLim);
set(ax1,'YDir','normal','XTick',500.5,'XTickLabel',{'-\pi/2'},'YTick',500.5,'YTickLabel',{'\pi'},'Color','none','FontSize',F);
set(ax2,'YDir','normal','XTick',500.5,'XTickLabel',{'\pi/2'},'YTick',500.5,'YTickLabel',{'0'},'Color','none','FontSize',F);

% Make a combined figure
l = 500;
ll = 73;
b = 15;
o = 0.06;
sy = (max(y) - min(y))/(max(x) - min(x));
sz = (max(z) - min(z))/(max(x) - min(x));
figure;
set(gcf,'Position',[460 50 l+sz*l+2*ll+b sy*l+sz*l+2*ll+b]);
axes('Units','pixels','Position',[ll ll+b+sz*l l sy*l]);
image(s*x,s*y,xy);
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-3,2.8,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+b+sz*l sz*l sy*l]);
image(s*z,s*y,zy);
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-2.2,2.8,'(b)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll l sz*l]);
image(s*x,s*z,xz);
xlabel(['x ' u]);
ylabel(['z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-3,1.98,'(c)','FontSize',20,'Color','white');
ax1 = axes('Units','pixels','Position',[ll+l+b+o*sz*l ll+o*sz*l (1-2*o)*sz*l (1-2*o)*sz*l]);
image(B);
ax2 = axes('Units','pixels','Position',[ll+l+b+o*sz*l ll+o*sz*l (1-2*o)*sz*l (1-2*o)*sz*l],'Color','none');
set(ax2,'XAxisLocation','top','YAxisLocation','Right','XLim',ax1.XLim,'YLim',ax1.YLim);
set(ax1,'YDir','normal','XTick',500.5,'XTickLabel',{'-\pi/2'},'YTick',500.5,'YTickLabel',{'\pi'},'FontSize',F);
set(ax2,'YDir','normal','XTick',500.5,'XTickLabel',{'\pi/2'},'YTick',500.5,'YTickLabel',{'0'},'FontSize',F);
text(ax1,30,930,'(d)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Reconstruction_Free_Space_2_Avg_XYZ.pdf'],'-dpdf','-r0');
    print([ps 'Reconstruction_Free_Space_2_Avg_XYZ.png'],'-dpng','-r600');
end

