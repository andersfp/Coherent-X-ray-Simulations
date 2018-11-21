% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the data sets
tic;
Ef = load_binary([p 'Ef.bin'],[ny nx n_omega]);
toc;
tic;
El = load_binary([p 'El.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Process the data
% Generate the axes
q1 = ((-ny/2):(ny/2 - 1)).'.*delta_q1;
q2 = ((-nx/2):(nx/2 - 1)).'.*delta_q2;
q3 = ((-n_omega/2):(n_omega/2 - 1)).'.*delta_q3;

% Get the intensities from non-rounded electric field
If = abs(Ef).^2;
Il = abs(El).^2;

% Normalize the intensities
Il = 1e6*Il./max(If(:));
If = 1e6*If./max(If(:));


%% Plot the data
% Extract slices to plot
If21 = log10(If(:,:,n_omega/2 + 1));
If23 = log10(squeeze(If(ny/2 + 1,:,:))).';
Il21 = log10(Il(:,:,n_omega/2 + 1));
Il23 = log10(squeeze(Il(ny/2 + 1,:,:))).';

% Set plotting font size
F = 41;
s = 1e-9;
u = '[nm^{-1}]';
%xt = -0.6:0.3:0.6;
%yt = -0.6:0.3:0.6;
%zt = -0.3:0.3:0.3;
xt = -0.5:0.5:0.5;
yt = -0.5:0.5:0.5;
zt = -0.3:0.3:0.3;
cl = [0 6];
ct1 = 0:6;
ctl1 = {'10^0','10^1','10^2','10^3','10^4','10^5','10^6'};
ct2 = 0:2:6;
ctl2 = {'10^0','10^2','10^4','10^6'};
ofs = 70;

% Plot q12 planes
figure;
set(gcf,'Position',[460 100 1000 1000-ofs]);
imagesc(s*q2,s*q1,If21,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
colorbar(gca,'Ticks',ct1,'TickLabels',ctl1);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'If21.pdf'],'-dpdf','-r0');
    print([ps 'If21.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000-ofs]);
imagesc(s*q2,s*q1,Il21,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
colorbar(gca,'Ticks',ct1,'TickLabels',ctl1);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Il21.pdf'],'-dpdf','-r0');
    print([ps 'Il21.png'],'-dpng','-r600');
end

% Plot q23 planes
figure;
set(gcf,'Position',[460 100 1000 500+ofs]);
imagesc(s*q2,s*q3,If23,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_3 ' u]);
colorbar(gca,'Ticks',ct2,'TickLabels',ctl2);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt,'OuterPosition',[0 0.05 1 1]);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'If23.pdf'],'-dpdf','-r0');
    print([ps 'If23.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 500+ofs]);
imagesc(s*q2,s*q3,Il23,cl);
axis equal tight;
colormap gray;
xlabel(['q_2 ' u]);
ylabel(['q_3 ' u]);
colorbar(gca,'Ticks',ct2,'TickLabels',ctl2);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt,'OuterPosition',[0 0.05 1 1]);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Il23.pdf'],'-dpdf','-r0');
    print([ps 'Il23.png'],'-dpng','-r600');
end

% Save two insets
[~,i11] = min(abs(q1 + 0.1e9));
[~,i12] = min(abs(q1 - 0.1e9));
[~,i21] = min(abs(q2 + 0.02e9));
[~,i22] = min(abs(q2 - 0.18e9));
A = flip(log10(If(i11:i12,i21:i22,n_omega/2 + 1))/6,1);
B = flip(log10(Il(i11:i12,i21:i22,n_omega/2 + 1))/6,1);
figure;imagesc(A,[0 1]);axis equal tight;colormap gray;
figure;imagesc(B,[0 1]);axis equal tight;colormap gray;
if sav == 1
    imwrite(A,[ps 'If21_Inset.png']);
    imwrite(B,[ps 'Il21_Inset.png']);
end

% Save two more insets
[~,i11] = min(abs(q3 + 0.1e9));
[~,i12] = min(abs(q3 - 0.1e9));
[~,i21] = min(abs(q2 + 0.02e9));
[~,i22] = min(abs(q2 - 0.18e9));
A = flip(If23(i11:i12,i21:i22)/6,1);
B = flip(Il23(i11:i12,i21:i22)/6,1);
A = interp2(A,1:124,linspace(1,139,124).');
B = interp2(B,1:124,linspace(1,139,124).');
figure;imagesc(A,[0 1]);axis equal tight;colormap gray;
figure;imagesc(B,[0 1]);axis equal tight;colormap gray;
if sav == 1
    imwrite(A,[ps 'If23_Inset.png']);
    imwrite(B,[ps 'Il23_Inset.png']);
end


