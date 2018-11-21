% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Get the object from free space propagation
tic;
F = load_binary([p 'Inversion_Object_F.bin'],[ny nx n_omega]);
toc;

% Get the object from the lens-based propagation
tic;
L = load_binary([p 'Inversion_Object_L.bin'],[ny nx n_omega]);
toc;

% Get the object from the hybrid simulation
tic;
H = load_binary([p 'Inversion_Object_H.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Process data
% Scale the objects to have comparable amplitudes
s0 = 2e7;
tic;
F = s0*F;
L = s0*L;
H = s0*H;
toc;

% Center the mean phase at 0
tic;
F = F.*exp(1i*pi/2);
L = L.*exp(1i*pi/2);
H = H.*exp(1i*pi/2);
% mpf = mean(angle(F(abs(F) > 0.5)));
% mpl = mean(angle(L(abs(L) > 0.5)));
% mph = mean(angle(H(abs(H) > 0.5)));
% F = F.*exp(-1i*mpf);
% L = L.*exp(-1i*mpl);
% H = H.*exp(-1i*mph);
toc;

% Plot the images
Slicer(abs(F),'displayRange',[0 1]);
Slicer(abs(L),'displayRange',[0 1]);
Slicer(abs(H),'displayRange',[0 1]);


%% Make plots
% Get the slices to plot
Fxy = F(:,:,n_omega/2 + 1);
Fzy = squeeze(F(:,nx/2 + 1,:));
Lxy = L(:,:,n_omega/2 + 1);
Lzy = squeeze(L(:,nx/2 + 1,:));
Hxy = H(:,:,n_omega/2 + 1);
Hzy = squeeze(H(:,nx/2 + 1,:));

% Scale the plots
mm = max([max(abs(Fxy(:))) max(abs(Fzy(:))) max(abs(Lxy(:))) max(abs(Lzy(:))) max(abs(Hxy(:))) max(abs(Hzy(:)))]);
Fxy = Fxy./mm;
Fzy = Fzy./mm;
Lxy = Lxy./mm;
Lzy = Lzy./mm;
Hxy = Hxy./mm;
Hzy = Hzy./mm;

% Make the color slices
cmap = hsv(256);
Fxy = complex2rgb(Fxy,cmap);
Fzy = complex2rgb(Fzy,cmap);
Lxy = complex2rgb(Lxy,cmap);
Lzy = complex2rgb(Lzy,cmap);
Hxy = complex2rgb(Hxy,cmap);
Hzy = complex2rgb(Hzy,cmap);

% Set plotting font size
fnt = 51;
s = 1e6;
u = '[\mum]';
xl = [-2 2];
yl = [-2 2];
zl = [-2 2];
% xt = -2:1.0:2;
% yt = -2:1.0:2;
% zt = -2:1.0:2;
xt = -2:2:2;
yt = -2:2:2;
zt = -2:2:2;
dcm = '0';
pos1 = [0.1743 0.21 0.7307 0.7504];

% Plot F planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,Fxy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Inversion_Fxy.pdf'],'-dpdf','-r0');
    print([ps 'Inversion_Fxy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,Fzy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Inversion_Fzy.pdf'],'-dpdf','-r0');
    print([ps 'Inversion_Fzy.png'],'-dpng','-r600');
end

% Plot L planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,Lxy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Inversion_Lxy.pdf'],'-dpdf','-r0');
    print([ps 'Inversion_Lxy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,Lzy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Inversion_Lzy.pdf'],'-dpdf','-r0');
    print([ps 'Inversion_Lzy.png'],'-dpng','-r600');
end

% Plot H planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*x,s*y,Hxy);
axis equal tight;
xlabel(['x ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Inversion_Hxy.pdf'],'-dpdf','-r0');
    print([ps 'Inversion_Hxy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*z,s*y,Hzy);
axis equal tight;
xlabel(['z ' u]);
ylabel(['y ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',zl,'XTick',zt,'XTickLabel',num2str(zt.',['%.' dcm 'f']),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.',['%.' dcm 'f']),'Position',pos1);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Inversion_Hzy.pdf'],'-dpdf','-r0');
    print([ps 'Inversion_Hzy.png'],'-dpng','-r600');
end



