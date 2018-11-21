% Initialization
clear;
close all;
clc;


%% Load the data
% Load .mat file
load('Si_Single_Crystal_Partial.mat');

% Save location
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\BFP\Figures\';

% Save?
sav = 0;


%% Make the plots
% Plot the image plane
s = 1e3;
u = ' [mm]';
F = 26;

figure;
set(gcf,'Position',[460 100 1060 1000]);
imagesc(s*xI,s*xI,IIc,[0 max(IIc(:))]);
axis equal tight;
xlabel(['x' u]);
ylabel(['y' u]);
set(gca,'FontSize',F,'XLim',[-3 3],'XTick',-2:1:2,'YLim',[-3 3],'YTick',-2:1:2,'YDir','normal');
if sav == 1
    print([ps 'IP_Coherent.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1060 1000]);
imagesc(s*xI,s*xI,IIp,[0 max(IIc(:))]);
axis equal tight;
xlabel(['x' u]);
ylabel(['y' u]);
set(gca,'FontSize',F,'XLim',[-3 3],'XTick',-2:1:2,'YLim',[-3 3],'YTick',-2:1:2,'YDir','normal');
if sav == 1
    print([ps 'IP_Partial.png'],'-dpng','-r600');
end

% Plot the BFP
s = 1e6;
u = ' [\mum]';

figure;
set(gcf,'Position',[460 100 1130 1000]);
imagesc(s*xB,s*xB,log10(IBc),log10(max(IBc(:))) + [-6 0]);
axis equal tight;
xlabel(['x' u]);
ylabel(['y' u]);
set(gca,'FontSize',F,'XLim',[-3 3],'XTick',-2:1:2,'YLim',[-3 3],'YTick',-2:1:2,'YDir','normal');
cb = colorbar;
%set(cb,'YTick',1:6,'YTickLabel',{'10^1','10^2','10^3','10^4','10^5','10^6'});
set(cb,'YTick',2:2:6,'YTickLabel',{'10^2','10^4','10^6'});
if sav == 1
    print([ps 'BFP_Coherent.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1130 1000]);
imagesc(s*xB,s*xB,log10(IBp),log10(max(IBc(:))) + [-6 0]);
axis equal tight;
xlabel(['x' u]);
ylabel(['y' u]);
set(gca,'FontSize',F,'XLim',[-3 3],'XTick',-2:1:2,'YLim',[-3 3],'YTick',-2:1:2,'YDir','normal');
cb = colorbar;
%set(cb,'YTick',1:6,'YTickLabel',{'10^1','10^2','10^3','10^4','10^5','10^6'});
set(cb,'YTick',2:2:6,'YTickLabel',{'10^2','10^4','10^6'});
if sav == 1
    print([ps 'BFP_Partial.png'],'-dpng','-r600');
end


