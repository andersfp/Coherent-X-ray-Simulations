% Initialization
clear;
close all;
clc;


%% Load data
% Load the lens exit field data
load('N_Fig_3_data.mat');

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Make plots
% Plotting parameters
F = 41;
sr = 1e3;
u = '[mm]';
% r1t = (-1.5:0.5:1.5).';
% r2t = (-1.5:0.5:1.5).';
r1t = (-1:1:1).';
r2t = (-1:1:1).';
dcm = '0';

% Plot the pupil
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(sr*r1,sr*r2,P,[0 1]);
xlabel(['r_2 ' u]);
ylabel(['r_1 ' u]);
axis equal tight;
colorbar(gca,'Ticks',0:0.2:1);
set(gca,'FontSize',F,'YDir','normal','XTick',r1t,'XTickLabel',num2str(r1t,['%.' dcm 'f']),'YTick',r2t,'YTickLabel',num2str(r2t,['%.' dcm 'f']));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Pupil_Function.pdf'],'-dpdf','-r0');
    print([ps 'Pupil_Function.png'],'-dpng','-r600');
end

% Plot the phase error
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(sr*r1,sr*r2,2*pi*cos(sqrt(r1.'.^2 + r2.^2)*2*pi/5e-4),[-2*pi 2*pi]);
xlabel(['r_2 ' u]);
ylabel(['r_1 ' u]);
axis equal tight;
%colorbar(gca,'Ticks',linspace(-2*pi,2*pi,9),'TickLabels',{'-2\pi','-3/2\pi','-\pi','-1/2\pi','0','1/2\pi','\pi','3/2\pi','2\pi'});
colorbar(gca,'Ticks',linspace(-2*pi,2*pi,5),'TickLabels',{'-2\pi','-\pi','0','\pi','2\pi'});
%colorbar(gca,'Ticks',linspace(-2*pi,2*pi,9),'TickLabels',{'$-2\pi$','$-\frac{3}{2}\pi$','$-\pi$','$-\frac{1}{2}\pi$','$0$','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$','$2\pi$'},'TickLabelInterpreter','latex');
set(gca,'FontSize',F,'YDir','normal','XTick',r1t,'XTickLabel',num2str(r1t,['%.' dcm 'f']),'YTick',r2t,'YTickLabel',num2str(r2t,['%.' dcm 'f']));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Phase_Errors.pdf'],'-dpdf','-r0');
    print([ps 'Phase_Errors.png'],'-dpng','-r600');
end

% Plot the lens field
Ef = Ef./max(abs(Ef(:)));
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(sr*r1,sr*r2,log10(abs(Ef)),[-3 0]);
colormap gray;
xlabel(['r_2 ' u]);
ylabel(['r_1 ' u]);
axis equal tight;
colorbar(gca,'Ticks',-3:0,'TickLabels',{'10^0','10^1','10^2','10^3'});
set(gca,'FontSize',F,'YDir','normal','XTick',r1t,'XTickLabel',num2str(r1t,['%.' dcm 'f']),'YTick',r2t,'YTickLabel',num2str(r2t,['%.' dcm 'f']));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Lens_Field.pdf'],'-dpdf','-r0');
    print([ps 'Lens_Field.png'],'-dpng','-r600');
end

% Plot the lens field after pupil/phase errors
B = B./max(abs(B(:)));
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(sr*r1,sr*r2,log10(abs(B)),[-3 0]);
colormap gray;
xlabel(['r_2 ' u]);
ylabel(['r_1 ' u]);
axis equal tight;
colorbar(gca,'Ticks',-3:0,'TickLabels',{'10^0','10^1','10^2','10^3'});
set(gca,'FontSize',F,'YDir','normal','XTick',r1t,'XTickLabel',num2str(r1t,['%.' dcm 'f']),'YTick',r2t,'YTickLabel',num2str(r2t,['%.' dcm 'f']));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Lens_Field_Pupil.pdf'],'-dpdf','-r0');
    print([ps 'Lens_Field_Pupil.png'],'-dpng','-r600');
end

% Save two insets
i11 = 471;
i12 = 595;
i21 = 504;
i22 = 627;
C = flip(log10(abs(Ef(i11:i12,i21:i22)))/3 + 1,1);
D = flip(log10(abs(B(i11:i12,i21:i22)))/3 + 1,1);
figure;imagesc(C,[0 1]);axis equal tight;colormap gray;
figure;imagesc(D,[0 1]);axis equal tight;colormap gray;
if sav == 1
    imwrite(C,[ps 'Lens_Field_Inset.png']);
    imwrite(D,[ps 'Lens_Field_Pupil_Inset.png']);
end


