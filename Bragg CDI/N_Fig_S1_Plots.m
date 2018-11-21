% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the true free space diffraction
tic;
Ef = load_binary([p 'Ef.bin'],[ny nx n_omega]);
toc;

% Load the true lens-based diffraction
tic;
El = load_binary([p 'El.bin'],[ny nx n_omega]);
toc;

% Load the high SNR free space diffraction
tic;
If1 = load_binary([p 'If1.bin'],[ny nx n_omega]);
toc;

% Load the low SNR free space diffraction
tic;
If2 = load_binary([p 'If2.bin'],[ny nx n_omega]);
toc;

% Load the high SNR lens-based diffraction
tic;
Il1 = load_binary([p 'Il1.bin'],[ny nx n_omega]);
toc;

% Load the low SNR lens-based diffraction
tic;
Il2 = load_binary([p 'Il2.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Process the data
% Calculate the true intensities
tic;
If0 = abs(Ef).^2;
Il0 = abs(El).^2;
toc;

% Scale the true intensities
tic;
Il0 = Il0./max(If0(:)).*1e6;
If0 = If0./max(If0(:)).*1e6;
toc;

% Extract the slices to plot
if0 = log10(If0(:,:,n_omega/2+1));
if1 = log10(If1(:,:,n_omega/2+1));
if2 = log10(If2(:,:,n_omega/2+1));
il0 = log10(Il0(:,:,n_omega/2+1));
il1 = log10(Il1(:,:,n_omega/2+1));
il2 = log10(Il2(:,:,n_omega/2+1));


%% Make the plots
% Set plotting font size
fnt = 51;
s = 1e-9;
u = '[nm^{-1}]';
xl = [-0.6 0.6];
yl = [-0.6 0.6];
xt = -0.6:0.6:0.6;
yt = -0.6:0.6:0.6;
cl = [-1 6];

% Plot infinite SNR free space intensity
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,if0,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_if0.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_if0.png'],'-dpng','-r600');
end

% Plot high SNR free space intensity
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,if1,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_if1.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_if1.png'],'-dpng','-r600');
end

% Plot low SNR free space intensity
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,if2,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_if2.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_if2.png'],'-dpng','-r600');
end

% Plot infinite SNR lens-based intensity
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,il0,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_il0.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_il0.png'],'-dpng','-r600');
end

% Plot high SNR lens-based intensity
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,il1,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_il1.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_il1.png'],'-dpng','-r600');
end

% Plot low SNR lens-based intensity
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,il2,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_il2.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_il2.png'],'-dpng','-r600');
end


