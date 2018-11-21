% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters
load('Sim_Parameters.mat');

% Save the figure?
s = 0;


%% Make the object
% Set up object parameters
m = 1000;
dx = 0.1e-3;
w = 1e-6;

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
%E0 = rect(x0/(2*w)).*rect(x0.'/(2*w)).*(x0 > 0);
S = load('Star.mat');
E0 = padarray(S.S,([m m]-size(S.S))/2);

% Make line and intensity from object
L0 = E0(:,m/2+1);
I0 = abs(E0).^2;


%% Calculate the wave propagation
% Make the propagation distance and focal length arrays
D1 = [d0+N*T/2;N*T/2+dd];
F1 = f/N;
DN = [d0+T/2;T*ones(N-1,1);T/2+dd];
FN = f*ones(N,1);

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);
[a1,Rm1,Rp1,sm1,sp1,gm1,gp1] = Lens_Stack(D1,F1,lambda,R0,s0);
[aN,RmN,RpN,smN,spN,gmN,gpN] = Lens_Stack(DN,FN,lambda,R0,s0);

% Initialize GPU
t = initGPU();

% Propagate the single lens
tic;
[E1,~,x1] = propFrFFT(E0,x0,x0,Rm1(1),Rp1(end),sm1(1),sp1(end),sum(a1),lambda,'gpu');
toc;
L1 = E1(:,m/2+1);
I1 = abs(E1).^2;

% Propagate the N lenses
tic;
[EN,~,xN] = propFrFFT(E0,x0,x0,RmN(1),RpN(end),smN(1),spN(end),sum(aN),lambda,'gpu');
toc;
LN = EN(:,m/2+1);
IN = abs(EN).^2;


%% Make plots
% Set plot limits
xl = 1e6*dx/2*[-1 1];
yl = 1e6*dx/2*[-1 1];
cm = max(I1(:));
ti = -500:250:500;
tx = -40;
ty = -40;
fs = 28;

% Object
figure;
imagesc(1e6*x0,1e6*x0,I0,[0 M^2*cm]);
axis equal tight;
title('Object');
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',xl,'YLim',yl,'XTick',ti/M,'YTick',ti/M);

% Single lens
figure;
imagesc(1e6*x1,1e6*x1,I1,[0 cm]);
axis equal tight;
title('Single lens');
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',M*xl,'YLim',M*yl,'XTick',ti,'YTick',ti);

% N lenses
figure;
imagesc(1e6*xN,1e6*xN,IN,[0 cm]);
axis equal tight;
title('N lenses');
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',M*xl,'YLim',M*yl,'XTick',ti,'YTick',ti);

% Plot the object vs N lenses difference
figure;
imagesc(1e6*x0,1e6*x0,I0 - IN*M^2,[-2e-4 2e-4]);
axis equal tight;
title('Object vs N lenses: Difference');
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',xl,'YLim',yl,'XTick',ti/M,'YTick',ti/M);

% All plots
figure;
set(gcf,'Position',[150 370 1600 400]);
subplot(1,4,1);
imagesc(1e6*x0,1e6*x0,I0,[0 M^2*cm]);
axis equal tight;
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',xl,'YLim',yl,'XTick',ti/M,'YTick',ti/M,'OuterPosition',[0 0 0.25 1]);
text(tx,ty,'(a)','FontSize',fs,'Color','white');
subplot(1,4,2);
imagesc(1e6*x1,1e6*x1,I1,[0 cm]);
axis equal tight;
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',M*xl,'YLim',M*yl,'XTick',ti,'YTick',ti,'OuterPosition',[0.25 0 0.25 1]);
text(M*tx,M*ty,'(b)','FontSize',fs,'Color','white');
subplot(1,4,3);
imagesc(1e6*xN,1e6*xN,IN,[0 cm]);
axis equal tight;
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',M*xl,'YLim',M*yl,'XTick',ti,'YTick',ti,'OuterPosition',[0.5 0 0.25 1]);
text(M*tx,M*ty,'(c)','FontSize',fs,'Color','white');
subplot(1,4,4);
imagesc(1e6*x0,1e6*x0,I0 - IN*M^2,[-2e-4 2e-4]);
axis equal tight;
xlabel('Position (\mum)');
ylabel('Position (\mum)');
set(gca,'XLim',xl,'YLim',yl,'XTick',ti/M,'YTick',ti/M,'OuterPosition',[0.75 0 0.25 1]);
text(tx,ty,'(d)','FontSize',fs,'Color','white');

% Save all plots
if s == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('Single_Thin_vs_N_Thin.pdf','-dpdf','-r0');
    print('Single_Thin_vs_N_Thin.png','-dpng','-r600');
end

