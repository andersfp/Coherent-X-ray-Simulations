% Initialization
clear;
close all;
clc;


%% Load the corrected images
% Set the data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the measurement
load([dp 'Focus_Scan_1.mat']);

% Load the simulation
sim = load([dp 'Star_Simulation_Series_1_Partial_Avg.mat']);
S = sim.I;


%% Save the images
% Select the images to save
ix = (525:1524) - 100;
iy = 525:1524;
i0 = 24;
s0 = 4;
d = 6;
I1 = I(iy,ix,i0 - d);
I2 = I(iy,ix,i0);
I3 = I(iy,ix,i0 + d);
S1 = S(iy,ix,s0 - d/3);
S2 = S(iy,ix,s0);
S3 = S(iy,ix,s0 + d/3);

% Normalize the images
imin = 0.7;
imax = 1.3;
I1 = (I1 - imin)./(imax - imin);
I2 = (I2 - imin)./(imax - imin);
I3 = (I3 - imin)./(imax - imin);
smin = 20000;
smax = 50000;
S1 = (S1 - smin)./(smax - smin);
S2 = (S2 - smin)./(smax - smin);
S3 = (S3 - smin)./(smax - smin);

% Plot the 3 images to save
cl = [0 1];
figure;
subplot(2,3,1);
imagesc(I1,cl);
colormap gray;
axis equal tight;
subplot(2,3,2);
imagesc(I2,cl);
colormap gray;
axis equal tight;
subplot(2,3,3);
imagesc(I3,cl);
colormap gray;
axis equal tight;
subplot(2,3,4);
imagesc(S1,cl);
colormap gray;
axis equal tight;
subplot(2,3,5);
imagesc(S2,cl);
colormap gray;
axis equal tight;
subplot(2,3,6);
imagesc(S3,cl);
colormap gray;
axis equal tight;

% Set the save path
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\MLLs\Figures\';
imwrite(I1,[ps 'I1.png']);
imwrite(I2,[ps 'I2.png']);
imwrite(I3,[ps 'I3.png']);
imwrite(S1,[ps 'S1.png']);
imwrite(S2,[ps 'S2.png']);
imwrite(S3,[ps 'S3.png']);


