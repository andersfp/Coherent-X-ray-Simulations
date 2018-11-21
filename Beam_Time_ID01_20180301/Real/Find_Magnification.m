% Initialization
clear;
close all;
clc;


%% Load data
% Load the diffraction patterns
tic;
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';
I = h5read([p 'Scan0008.cxi'],'/entry_1/instrument_1/detector_1/data');
L = h5read([p 'Scan0267.cxi'],'/entry_1/instrument_1/detector_1/data');
toc;

% Convert data
I = double(permute(I,[2 1 3]));
L = double(permute(L,[2 1 3]));


%% Compare data
% Plot 3D image of the two data sets
Slicer(log10(I),'displayRange',[0 5]);
Slicer(log10(L),'displayRange',[0 5]);

% Plot the central image
figure;
subplot(2,2,1);
imagesc(I(:,:,141));
axis equal tight;
subplot(2,2,2);
imagesc(L(:,:,51));
axis equal tight;
subplot(2,2,3);
imagesc(log10(I(:,:,141)));
axis equal tight;
subplot(2,2,4);
imagesc(log10(L(:,:,51)));
axis equal tight;

% Get lines
x = (1:516) - 182;
y = (1:516).' - 204;
% ix = I(206,:,141);
% iy = I(:,181,141);
% lx = L(202,:,141);
% ly = L(:,183,141);
w = 3;
ix = sum(I(204-w:204+w,:,141),1);
iy = sum(I(:,182-w:182+w,141),2);
lx = sum(L(204-w:204+w,:,51),1);
ly = sum(L(:,182-w:182+w,51),2);

% Plot the lines
figure;
plot(x,log10(ix),x,log10(lx));
figure;
plot(y,log10(iy),y,log10(ly));

% Fit the magnification
iy = iy/max(iy);
ly = ly/max(ly);
fun = @(x) sum(abs(iy - interp1(y,ly,x(1).*(y + x(2)),'linear',4)).^2);
m = fminsearch(fun,[1 0]);
M = 1./m(1);

% Plot the fit
ly2 = interp1(y,ly,m(1).*(y + m(2)),'linear',4);
figure;
plot(y,log10(iy),y,log10(ly2));

% % Fit the pupil function
% fun2 = @(x) sum(abs((gaussRMS(y,x).*iy + 3e-6) - ly2).^2);
% sigma = fminsearch(fun2,30);
% 
% % Plot the fit
% figure;
% plot(y,log10(gaussRMS(y,sigma).*iy + 3e-6),y,log10(ly2));



% [iye1,iye2] = envelope(log10(iy + 1e-6));
% [lye1,lye2] = envelope(log10(ly2 + 1e-6));
% figure;
% plot(y,[log10(iy) iye1 iye2 log10(ly2) lye1 lye2]);

figure;
subplot(2,2,1);
imagesc(sum(I,3));
axis equal tight;
subplot(2,2,2);
imagesc(sum(L,3));
axis equal tight;
subplot(2,2,3);
imagesc(log10(sum(I,3)));
axis equal tight;
subplot(2,2,4);
imagesc(log10(sum(L,3)));
axis equal tight;

figure;
plot(y,log10(gaussRMS(y,24.2).*iy),y,log10(ly2));


