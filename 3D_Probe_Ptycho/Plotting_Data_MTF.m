% Initialization
clear;
close all;
clc;


%% Load data
% Set the full file locations
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\';
f = 'Plotting_Data.mat';

% Load the data
load([p f]);


%% Extract MTFs
% Perform FFT
% fobj = fftshift(fftn(fftshift(obj)));
% frec = fftshift(fftn(fftshift(rec)));
% fna4 = fftshift(fftn(fftshift(na4)));
fobj = fftshift(fftn(fftshift(abs(obj))));
frec = fftshift(fftn(fftshift(abs(rec))));
fna4 = fftshift(fftn(fftshift(abs(na4))));

% Extract line plots
ox = abs(fobj(129,:,129)).';
oy = abs(fobj(:,129,129));
oz = squeeze(abs(fobj(129,129,:)));
rx = abs(frec(129,:,129)).';
ry = abs(frec(:,129,129));
rz = squeeze(abs(frec(129,129,:)));
nx = abs(fna4(129,:,129)).';
ny = abs(fna4(:,129,129));
nz = squeeze(abs(fna4(129,129,:)));

% Plot the line plots
figure;
semilogy(y,[ox rx nx]);
figure;
semilogy(y,[oy ry ny]);
figure;
semilogy(y,[oz rz nz]);

% Fit pupil function
fun1 = @(ref,w,y) ref.*gaussRMS(y,w);
fun2 = @(ref,test,w,y) sum((fun1(ref,w,y) - test).^2);
fun3 = @(ref,test,w,y) sum((log10(fun1(ref,w,y)) - log10(test)).^2);
wrx = fminsearch(@(w) fun2(ox,rx,w,y),31);
wry = fminsearch(@(w) fun2(oy,ry,w,y),29);
%wrz = fminsearch(@(w) fun2(oz,rz,w,y),3);
%wrz = fminsearch(@(w) fun2(oz(abs(y) < 50),rz(abs(y) < 50),w,y(abs(y) < 50)),10);
%wrz = fminsearch(@(w) fun3(oz,rz,w,y),10);
wrz = fminsearch(@(w) fun3(oz(abs(y) < 20),rz(abs(y) < 20),w,y(abs(y) < 20)),6);
wnx = fminsearch(@(w) fun2(ox,nx,w,y),51);
wny = fminsearch(@(w) fun2(oy,ny,w,y),94);
%wnz = fminsearch(@(w) fun2(oz,nz,w,y),3);
%wnz = fminsearch(@(w) fun2(oz(abs(y) < 50),nz(abs(y) < 50),w,y(abs(y) < 50)),10);
wnz = fminsearch(@(w) fun3(oz(abs(y) < 20),nz(abs(y) < 20),w,y(abs(y) < 20)),6);

% Plot the fits
figure;
semilogy(y,[ox rx fun1(ox,wrx,y)]);
figure;
semilogy(y,[oy ry fun1(oy,wry,y)]);
figure;
semilogy(y,[oz rz fun1(oz,wrz,y)]);
figure;
semilogy(y,[ox nx fun1(ox,wnx,y)]);
figure;
semilogy(y,[oy ny fun1(oy,wny,y)]);
figure;
semilogy(y,[oz nz fun1(oz,wnz,y)]);

% Reciprocal space pixel size
dq = 1./5e-6;

% Reciprocal space width
W = 2.*pi.*dq.*[wrx wry wrz;wnx wny wnz];

% Estimate real space resolution
R = 1./W;


%% Make plot
% Generate the axis
q = 2.*pi.*dq.*y;

% Make the plot
sc = 1e-6;
xlim = [-150 150];
ylim = [7 5e6];
fnt = 20;
lw = 2;
figure('Position',[250 310 1500 550]);
axes('OuterPosition',[0 0 1/3 1]);
semilogy(sc.*q,[ox rx nx],'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim,'FontSize',fnt);
xlabel('q_x [\mum^{-1}]');
ylabel('Amplitude (a.u.)');
text(-140,2.3e6,'(a)','FontSize',28);
axes('OuterPosition',[1/3 0 1/3 1]);
semilogy(sc.*q,[oy ry ny],'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim,'FontSize',fnt);
xlabel('q_y [\mum^{-1}]');
ylabel('Amplitude (a.u.)');
text(-140,2.3e6,'(b)','FontSize',28);
axes('OuterPosition',[2/3 0 1/3 1]);
semilogy(sc.*q,[oz rz nz],'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim,'FontSize',fnt);
xlabel('q_z [\mum^{-1}]');
ylabel('Amplitude (a.u.)');
text(-140,2.3e6,'(c)','FontSize',28);
%legend('Ground truth','Reconstruction','Reconstruction, 4\timesNA');
%legend('Ground truth','Rec.','Rec. 4\timesNA');
legend('Gr. truth','Rec.','4\timesNA');

% Save the plot
%print(gcf,'MTF_Figure.png','-dpng','-r480');


