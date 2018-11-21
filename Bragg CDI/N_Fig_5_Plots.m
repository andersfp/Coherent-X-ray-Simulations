% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Get the true object
tic;
O = load_binary([p 'Oxyz.bin'],[ny nx n_omega]);
toc;

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

% Find the center of mass
tic;
xo = sum(sum(sum(abs(O).*(x.'))))./sum(abs(O(:)));
yo = sum(sum(sum(abs(O).*(y))))./sum(abs(O(:)));
zo = sum(sum(sum(abs(O).*(permute(z,[3 2 1])))))./sum(abs(O(:)));
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
ixo = -round(xo/mean(diff(x)));
iyo = -round(yo/mean(diff(y)));
izo = -round(zo/mean(diff(z)));
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
O = circshift(O,[iyo ixo izo]);
F1 = circshift(F1,[iyf1 ixf1 izf1]);
F2 = circshift(F2,[iyf2 ixf2 izf2]);
L1 = circshift(L1,[iyl1 ixl1 izl1]);
L2 = circshift(L2,[iyl2 ixl2 izl2]);
toc;

% Find the center volume
[X,Y,Z] = meshgrid(x,y,z);
ii = ((abs(X) <= 0.4e-6) & (abs(Y) <= 0.4e-6) & (abs(Z) <= 0.4e-6) & (Y <= 0) | (X.^2 + Y.^2 <= (0.4e-6).^2) & (abs(Z) <= 0.4e-6) & (Y >= 0));

% Center the mean phase at 0
tic;
F1 = F1.*exp(1i*1.74);
F2 = F2.*exp(1i*2.13);
L1 = L1.*exp(1i*0.19);
L2 = L2.*exp(1i*2.04);
% mpf1 = mean(angle(F1(F1 ~= 0)));
% mpf2 = mean(angle(F2(F2 ~= 0)));
% mpl1 = mean(angle(L1(L1 ~= 0)));
% mpl2 = mean(angle(L2(L2 ~= 0)));
mpf1 = mean(angle(F1(ii)));
mpf2 = mean(angle(F2(ii)));
mpl1 = mean(angle(L1(ii)));
mpl2 = mean(angle(L2(ii)));
F1 = F1.*exp(-1i*mpf1);
F2 = F2.*exp(-1i*mpf2);
L1 = L1.*exp(-1i*mpl1);
L2 = L2.*exp(-1i*mpl2);
toc;

% Set the average amplitude
tic;
mF1 = mean(abs(F1(ii)));
mF2 = mean(abs(F2(ii)));
mL1 = mean(abs(L1(ii)));
mL2 = mean(abs(L2(ii)));
F1 = F1./mF1;
F2 = F2./mF2;
L1 = L1./mL1;
L2 = L2./mL2;
toc;

% Plot the images
Slicer(O,'displayRange',[0 1]);
Slicer(abs(F1),'displayRange',[0 1]);
Slicer(abs(F2),'displayRange',[0 1]);
Slicer(abs(L1),'displayRange',[0 1]);
Slicer(abs(L2),'displayRange',[0 1]);


%% Make plots
% Get the slices to plot
Ozy = squeeze(O(:,nx/2 + 1,:));
F1zy = squeeze(F1(:,nx/2 + 1,:));
F2zy = squeeze(F2(:,nx/2 + 1,:));
L1zy = squeeze(L1(:,nx/2 + 1,:));
L2zy = squeeze(L2(:,nx/2 + 1,:));

% Scale the plots
%mm = mean([mean(abs(F1zy(abs(F1zy) > 0))) mean(abs(F2zy(abs(F2zy) > 0))) mean(abs(L1zy(abs(L1zy) > 0))) mean(abs(L2zy(abs(L2zy) > 0)))]);
% mm = mean([sum(abs(F1zy(:))) sum(abs(F1zy(:))) sum(abs(F1zy(:))) sum(abs(F1zy(:)))])./sum(Ozy(:));
% F1zy = F1zy./mm;
% F2zy = F2zy./mm;
% L1zy = L1zy./mm;
% L2zy = L2zy./mm;

% Center the slice phase at 0
% mpf1zy = mean(angle(F1zy(F1zy ~= 0)));
% mpf2zy = mean(angle(F2zy(F2zy ~= 0)));
% mpl1zy = mean(angle(L1zy(L1zy ~= 0)));
% mpl2zy = mean(angle(L2zy(L2zy ~= 0)));
% F1zy = F1zy.*exp(-1i*mpf1zy);
% F2zy = F2zy.*exp(-1i*mpf1zy);
% L1zy = L1zy.*exp(-1i*mpf1zy);
% L2zy = L2zy.*exp(-1i*mpf1zy);

% Make the difference plots
af1 = abs(F1zy) - Ozy;
pf1 = angle(F1zy);
af2 = abs(F2zy) - Ozy;
pf2 = angle(F2zy);
al1 = abs(L1zy) - Ozy;
pl1 = angle(L1zy);
al2 = abs(L2zy) - Ozy;
pl2 = angle(L2zy);

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
al = [-0.1 0.1];
pl = [-pi/20 pi/20];
w = 1100;
h = 1000;
pos1 = [0.1412 0.21 0.7638 0.7304];
pos2 = [0.1278 0.21 0.7054 0.7304];

% Plot F1 differences
figure;
set(gcf,'Position',[460 100 w h]);
imagesc(s*z,s*y,af1,al);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'AF1.pdf'],'-dpdf','-r0');
    print([ps 'AF1.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 w h]);
imagesc(s*z,s*y,pf1,pl);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'PF1.pdf'],'-dpdf','-r0');
    print([ps 'PF1.png'],'-dpng','-r600');
end

% Plot F2 differences
figure;
set(gcf,'Position',[460 100 w h]);
imagesc(s*z,s*y,af2,al);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'AF2.pdf'],'-dpdf','-r0');
    print([ps 'AF2.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 w h]);
imagesc(s*z,s*y,pf2,pl);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'PF2.pdf'],'-dpdf','-r0');
    print([ps 'PF2.png'],'-dpng','-r600');
end

% Plot L1 differences
figure;
set(gcf,'Position',[460 100 w h]);
imagesc(s*z,s*y,al1,al);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'AL1.pdf'],'-dpdf','-r0');
    print([ps 'AL1.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 w h]);
imagesc(s*z,s*y,pl1,pl);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'PL1.pdf'],'-dpdf','-r0');
    print([ps 'PL1.png'],'-dpng','-r600');
end

% Plot L2 differences
figure;
set(gcf,'Position',[460 100 w+100 h]);
imagesc(s*z,s*y,al2,al);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos2);
colorbar(gca,'Ticks',linspace(al(1),al(2),5));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'AL2.pdf'],'-dpdf','-r0');
    print([ps 'AL2.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 w+100 h]);
imagesc(s*z,s*y,pl2,pl);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',F,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos2);
%colorbar(gca,'Ticks',linspace(pl(1),pl(2),9),'TickLabels',{'-\pi/20','-3\pi/80','-\pi/40','-\pi/80','0','\pi/80','\pi/40','3\pi/80','\pi/20'});
colorbar(gca,'Ticks',linspace(pl(1),pl(2),5),'TickLabels',{'-\pi/20','-\pi/40','0','\pi/40','\pi/20'});
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'PL2.pdf'],'-dpdf','-r0');
    print([ps 'PL2.png'],'-dpng','-r600');
end


