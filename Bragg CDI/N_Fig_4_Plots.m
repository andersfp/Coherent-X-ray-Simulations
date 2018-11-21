% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Get the averaged object reconstructed from 1000000 counts
tic;
F1 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 50000 counts
tic;
F2 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 1000000 counts
tic;
L1 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Get the averaged object reconstructed from 50000 counts
tic;
L2 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_2_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Process data
% Scale the objects to have comparable amplitudes
s0 = 100;
tic;
F1 = s0*F1;
F2 = (s0*sqrt(20))*F2;
L1 = s0*L1;
L2 = (s0*sqrt(20))*L2;
toc;

% Center the mean phase at 0
tic;
F1 = F1.*exp(1i*1.74);
F2 = F2.*exp(1i*2.13);
L1 = L1.*exp(1i*0.19);
L2 = L2.*exp(1i*2.04);
mpf1 = mean(angle(F1(F1 ~= 0)));
mpf2 = mean(angle(F2(F2 ~= 0)));
mpl1 = mean(angle(L1(L1 ~= 0)));
mpl2 = mean(angle(L2(L2 ~= 0)));
F1 = F1.*exp(-1i*mpf1);
F2 = F2.*exp(-1i*mpf2);
L1 = L1.*exp(-1i*mpl1);
L2 = L2.*exp(-1i*mpl2);
toc;

% Find the center of mass
tic;
xf1 = sum(sum(sum(abs(F1).*(x.'))))./sum(abs(F1(:)));
yf1 = sum(sum(sum(abs(F1).*(y))))./sum(abs(F1(:)));
zf1 = sum(sum(sum(abs(F1).*(permute(z,[3 2 1])))))./sum(abs(F1(:)));
xf2 = sum(sum(sum(abs(F2).*(x.'))))./sum(abs(F2(:)));
yf2 = sum(sum(sum(abs(F2).*(y))))./sum(abs(F2(:)));
zf2 = sum(sum(sum(abs(F2).*(permute(z,[3 2 1])))))./sum(abs(F2(:)));
xl1 = sum(sum(sum(abs(L1).*(x.'))))./sum(abs(L1(:)));
yl1 = sum(sum(sum(abs(L1).*(y))))./sum(abs(L1(:)));
zl1 = sum(sum(sum(abs(L1).*(permute(z,[3 2 1])))))./sum(abs(L1(:)));
xl2 = sum(sum(sum(abs(L2).*(x.'))))./sum(abs(L2(:)));
yl2 = sum(sum(sum(abs(L2).*(y))))./sum(abs(L2(:)));
zl2 = sum(sum(sum(abs(L2).*(permute(z,[3 2 1])))))./sum(abs(L2(:)));
toc;

% Find the amount to shift
ixf1 = -round(xf1/mean(diff(x)));
iyf1 = -round(yf1/mean(diff(y)));
izf1 = -round(zf1/mean(diff(z)));
ixf2 = -round(xf2/mean(diff(x)));
iyf2 = -round(yf2/mean(diff(y)));
izf2 = -round(zf2/mean(diff(z)));
ixl1 = -round(xl1/mean(diff(x)));
iyl1 = -round(yl1/mean(diff(y)));
izl1 = -round(zl1/mean(diff(z)));
ixl2 = -round(xl2/mean(diff(x)));
iyl2 = -round(yl2/mean(diff(y)));
izl2 = -round(zl2/mean(diff(z)));

% Shift the images
tic;
F1 = circshift(F1,[iyf1 ixf1 izf1]);
F2 = circshift(F2,[iyf2 ixf2 izf2]);
L1 = circshift(L1,[iyl1 ixl1 izl1]);
L2 = circshift(L2,[iyl2 ixl2 izl2]);
toc;

% Plot the images
Slicer(abs(F1),'displayRange',[0 1]);
Slicer(abs(F2),'displayRange',[0 1]);
Slicer(abs(L1),'displayRange',[0 1]);
Slicer(abs(L2),'displayRange',[0 1]);


%% Make plots
% Get the slices to plot
F1xy = F1(:,:,n_omega/2 + 1);
F1zy = squeeze(F1(:,nx/2 + 1,:));
F2xy = F2(:,:,n_omega/2 + 1);
F2zy = squeeze(F2(:,nx/2 + 1,:));
L1xy = L1(:,:,n_omega/2 + 1);
L1zy = squeeze(L1(:,nx/2 + 1,:));
L2xy = L2(:,:,n_omega/2 + 1);
L2zy = squeeze(L2(:,nx/2 + 1,:));

% Scale the plots
mm = max([max(abs(F1xy(:))) max(abs(F1zy(:))) max(abs(F2xy(:))) max(abs(F2zy(:))) max(abs(L1xy(:))) max(abs(L1zy(:))) max(abs(L2xy(:))) max(abs(L2zy(:)))]);
F1xy = F1xy./mm;
F1zy = F1zy./mm;
F2xy = F2xy./mm;
F2zy = F2zy./mm;
L1xy = L1xy./mm;
L1zy = L1zy./mm;
L2xy = L2xy./mm;
L2zy = L2zy./mm;

% Center the slice phase at 0
% mpf1xy = mean(angle(F1xy(F1xy ~= 0)));
% mpf1zy = mean(angle(F1zy(F1zy ~= 0)));
% mpf2xy = mean(angle(F2xy(F2xy ~= 0)));
% mpf2zy = mean(angle(F2zy(F2zy ~= 0)));
% mpl1xy = mean(angle(L1xy(L1xy ~= 0)));
% mpl1zy = mean(angle(L1zy(L1zy ~= 0)));
% mpl2xy = mean(angle(L2xy(L2xy ~= 0)));
% mpl2zy = mean(angle(L2zy(L2zy ~= 0)));
% F1xy = F1xy.*exp(-1i*mpf1xy);
% F1zy = F1zy.*exp(-1i*mpf1zy);
% F2xy = F2xy.*exp(-1i*mpf1xy);
% F2zy = F2zy.*exp(-1i*mpf1zy);
% L1xy = L1xy.*exp(-1i*mpf1xy);
% L1zy = L1zy.*exp(-1i*mpf1zy);
% L2xy = L2xy.*exp(-1i*mpf1xy);
% L2zy = L2zy.*exp(-1i*mpf1zy);

% Make the color slices
cmap = hsv(256);
F1xy = complex2rgb(F1xy,cmap);
F1zy = complex2rgb(F1zy,cmap);
F2xy = complex2rgb(F2xy,cmap);
F2zy = complex2rgb(F2zy,cmap);
L1xy = complex2rgb(L1xy,cmap);
L1zy = complex2rgb(L1zy,cmap);
L2xy = complex2rgb(L2xy,cmap);
L2zy = complex2rgb(L2zy,cmap);

% Set plotting font size
F = 51;
s = 1e6;
u = '[\mum]';
xl = [-1 1];
yl = [-1 1];
zl = [-1 1];
% xt = -1:0.5:1;
% yt = -1:0.5:1;
% zt = -1:0.5:1;
xt = -1:1:1;
yt = -1:1:1;
zt = -1:1:1;
dcm = '0';
pos1 = [0.1743 0.21 0.7307 0.7504];

% Plot F1 planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,F1xy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'F1xy.pdf'],'-dpdf','-r0');
    print([ps 'F1xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,F1zy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'F1zy.pdf'],'-dpdf','-r0');
    print([ps 'F1zy.png'],'-dpng','-r600');
end

% Plot F2 planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,F2xy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'F2xy.pdf'],'-dpdf','-r0');
    print([ps 'F2xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,F2zy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'F2zy.pdf'],'-dpdf','-r0');
    print([ps 'F2zy.png'],'-dpng','-r600');
end

% Plot L1 planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,L1xy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'L1xy.pdf'],'-dpdf','-r0');
    print([ps 'L1xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,L1zy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'L1zy.pdf'],'-dpdf','-r0');
    print([ps 'L1zy.png'],'-dpng','-r600');
end

% Plot L2 planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,L2xy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'L2xy.pdf'],'-dpdf','-r0');
    print([ps 'L2xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,L2zy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'L2zy.pdf'],'-dpdf','-r0');
    print([ps 'L2zy.png'],'-dpng','-r600');
end


