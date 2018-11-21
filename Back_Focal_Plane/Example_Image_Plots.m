% Initialization
clear;
close all;
clc;


%% Load data
% Load the simulation settings
load('Example_Image_Settings.mat');

% Load the coherent simulation
C = load('Example_Image_Coherent.mat');

% Load the partial transversal coherence simulation
T = load('Example_Image_Transversal.mat');

% Load the partial longitudinal coherence simulation
L = load('Example_Image_Longitudinal.mat');

% Load the partial coherence simulation
LT = load('Example_Image_LT.mat');

% Get the plotting coordinates
x = C.x;
xB = C.xB;


%% Make plots
% Make limits for the plots
blim = max(log10(C.B(:))) + [-6 0];
ilim = [0 0.01];

% Set plotting scales
ps = 1e6;
u = ' [\mum]';
fnt = 20;

% Plot the BFP in log scale
figure;
subplot(2,2,1);
imagesc(ps*xB,ps*xB,log10(C.B),blim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0 0.5 0.5 0.5]);
title('Coherent');
xlabel(['x' u]);
ylabel(['y' u]);
subplot(2,2,2);
imagesc(ps*xB,ps*xB,log10(L.B),blim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0.5 0.5 0.5 0.5]);
title('Partial longitudinal coherence');
xlabel(['x' u]);
ylabel(['y' u]);
subplot(2,2,3);
imagesc(ps*xB,ps*xB,log10(T.B),blim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0 0 0.5 0.5]);
title('Partial transversal coherence');
xlabel(['x' u]);
ylabel(['y' u]);
subplot(2,2,4);
imagesc(ps*xB,ps*xB,log10(LT.B),blim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0.5 0 0.5 0.5]);
title('Partial coherence');
xlabel(['x' u]);
ylabel(['y' u]);

% Make line plot through center of BFP
figure;
semilogy(ps.*xB,[C.B(:,m/2+1) L.B(:,m/2+1) T.B(:,m/2+1) LT.B(:,m/2+1)]);
set(gca,'FontSize',fnt,'XLim',[-300 300]);
title('Line through center of BFP');
xlabel(['x' u]);
ylabel('Intensity [a.u.]');
legend('Coherent','Longitudinal partial','Transverse partial','Partial');

% Make line plot through center of BFP
figure;
semilogy(ps.*xB,[C.B(:,m/2) L.B(:,m/2) T.B(:,m/2) LT.B(:,m/2)]);
set(gca,'FontSize',fnt,'XLim',[-300 300]);
title('Line through 1 pixel off center of BFP');
xlabel(['x' u]);
ylabel('Intensity [a.u.]');
legend('Coherent','Longitudinal partial','Transverse partial','Partial');

% Make line plot of 1D sum of BFP
figure;
semilogy(ps.*xB,[sum(C.B,2,'omitnan') sum(L.B,2,'omitnan') sum(T.B,2,'omitnan') sum(LT.B,2,'omitnan')]);
set(gca,'FontSize',fnt,'XLim',[-300 300]);
title('1D integrated intensity of BFP');
xlabel(['x' u]);
ylabel('Intensity [a.u.]');
legend('Coherent','Longitudinal partial','Transverse partial','Partial');

% Plot the image plane
figure;
subplot(2,2,1);
imagesc(ps*x,ps*x,C.I,ilim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0 0.5 0.5 0.5]);
title('Coherent');
xlabel(['x' u]);
ylabel(['y' u]);
subplot(2,2,2);
imagesc(ps*x,ps*x,L.I,ilim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0.5 0.5 0.5 0.5]);
title('Partial longitudinal coherence');
xlabel(['x' u]);
ylabel(['y' u]);
subplot(2,2,3);
imagesc(ps*x,ps*x,T.I,ilim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0 0 0.5 0.5]);
title('Partial transversal coherence');
xlabel(['x' u]);
ylabel(['y' u]);
subplot(2,2,4);
imagesc(ps*x,ps*x,LT.I,ilim);
axis equal tight;
set(gca,'YDir','normal','FontSize',fnt,'OuterPosition',[0.5 0 0.5 0.5]);
title('Partial coherence');
xlabel(['x' u]);
ylabel(['y' u]);


