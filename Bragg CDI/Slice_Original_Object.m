% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the original object
O = load_binary([p 'Oxyz.bin'],[ny nx n_omega]);

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Plot the data
% Extract the 3 slices
xy = O(:,:,n_omega/2 + 1);
zy = squeeze(O(:,nx/2 + 1,:));
xz = squeeze(O(ny/2 + 1,:,:)).';

% Make the color slices
cmap = hsv(256);
xy = complex2rgb(xy,cmap);
zy = complex2rgb(zy,cmap);
xz = complex2rgb(xz,cmap);

% Set plotting font size
F = 50;
s = 1e6;
u = '[\mum]';
xt = -1:1:1;
yt = -1:1:1;
zt = -2:1:2;
dcm = '0';
ofsx = 0.20;
ofsy = 0.20;
xl = 1.5;
yl = 1.5;
zl = 1.5;

% Plot the xy slice
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,xy);
xlabel(['x ' u]);
ylabel(['y ' u]);
hxy = gca;
hxy.FontSize = F;
hxy.YDir = 'normal';
hxy.XTick = xt;
hxy.XTickLabel = num2str(xt.',['%.' dcm 'f']);
hxy.YTick = yt;
hxy.YTickLabel = num2str(yt.',['%.' dcm 'f']);
hxy.Position = [ofsx ofsy 0.95-ofsx 0.95-ofsy];
hxy.XLim = [-xl xl];
hxy.YLim = [-yl yl];
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Original_Object_XY.pdf'],'-dpdf','-r0');
    print([ps 'Original_Object_XY.png'],'-dpng','-r600');
end

% Plot the zy slice
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,zy);
xlabel(['z ' u]);
ylabel(['y ' u]);
hzy = gca;
hzy.FontSize = F;
hzy.YDir = 'normal';
hzy.XTick = zt;
hzy.XTickLabel = num2str(zt.',['%.' dcm 'f']);
hzy.YTick = yt;
hzy.YTickLabel = num2str(yt.',['%.' dcm 'f']);
hzy.Position = [ofsx ofsy 0.95-ofsx 0.95-ofsy];
hzy.XLim = [-zl zl];
hzy.YLim = [-yl yl];
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Original_Object_ZY.pdf'],'-dpdf','-r0');
    print([ps 'Original_Object_ZY.png'],'-dpng','-r600');
end

% Plot the xz slice
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*z,xz);
xlabel(['x ' u]);
ylabel(['z ' u]);
hxz = gca;
hxz.FontSize = F;
hxz.YDir = 'normal';
hxz.XTick = xt;
hxz.XTickLabel = num2str(xt.',['%.' dcm 'f']);
hxz.YTick = zt;
hxz.YTickLabel = num2str(zt.',['%.' dcm 'f']);
hxz.Position = [ofsx ofsy 0.95-ofsx 0.95-ofsy];
hxz.XLim = [-xl xl];
hxz.YLim = [-zl zl];
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Original_Object_XZ.pdf'],'-dpdf','-r0');
    print([ps 'Original_Object_XZ.png'],'-dpng','-r600');
end

% Generate color-wheel
[X,Y] = meshgrid(linspace(-1,1,1000),linspace(-1,1,1000));
A = sqrt(X.^2 + Y.^2).*double(hypot(X,Y) <= 0.95).*exp(1i*atan2(Y,X));
B = complex2rgb(A,cmap);
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(B);
ax1 = gca;
axis equal tight;
ax2 = axes('Position',ax1.Position,'Color', 'none');
axis equal tight;
set(ax2,'XAxisLocation','top','YAxisLocation','Right');
set(ax2,'XLim',ax1.XLim,'YLim',ax1.YLim);
set(ax1,'YDir','normal','XTick',500.5,'XTickLabel',{'-\pi/2'},'YTick',500.5,'YTickLabel',{'\pi'},'Color','none','FontSize',F);
set(ax2,'YDir','normal','XTick',500.5,'XTickLabel',{'\pi/2'},'YTick',500.5,'YTickLabel',{'0'},'Color','none','FontSize',F);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Color_Wheel.pdf'],'-dpdf','-r0');
    print([ps 'Color_Wheel.png'],'-dpng','-r600');
end

% Make a combined figure
l = 400;
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
text(-1.85,1.75,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+b+sz*l sz*l sy*l]);
image(s*z,s*y,zy);
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-2.1,1.75,'(b)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll l sz*l]);
image(s*x,s*z,xz);
xlabel(['x ' u]);
ylabel(['z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-1.85,1.98,'(c)','FontSize',20,'Color','white');
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
    %print([ps 'Original_Object_XYZ.pdf'],'-dpdf','-r0');
    print([ps 'Original_Object_XYZ.png'],'-dpng','-r600');
end

