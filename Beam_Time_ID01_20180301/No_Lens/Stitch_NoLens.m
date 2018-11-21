% Initialization
clear;
close all;
clc;


%% Set the file names
% Set the data path
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';

% File names
f0 = 'Scan0008.cxi';
f1 = 'Scan0009.cxi';
f2 = 'Scan0010.cxi';


%% Read the data
% Set the data path
dp = '/entry_1/instrument_1/detector_1/data';

% Read the data
A0 = h5read([p f0],dp);
A1 = h5read([p f1],dp);
A2 = h5read([p f2],dp);


%% Process the data
% Convert data to double precision and transpose
A0 = permute(double(A0),[2 1 3]);
A1 = permute(double(A1),[2 1 3]);
A2 = permute(double(A2),[2 1 3]);

% Set the central slice index
ii = 141;

% Plot the middle sections
figure;
subplot(1,3,1);
imagesc(log10(A0(:,:,ii)),[1 6]);
axis equal tight;
subplot(1,3,2);
imagesc(log10(A1(:,:,ii)),[1 6]);
axis equal tight;
subplot(1,3,3);
imagesc(log10(A2(:,:,ii)),[1 6]);
axis equal tight;

% Extract central slices for further processing
a0 = A0(:,:,ii);
a1 = A1(:,:,ii);
a2 = A2(:,:,ii);

% Add zeropadding
z = zeros(516);
a0 = cat(1,z,a0,z);
a1 = cat(1,z,a1,z);
a2 = cat(1,z,a2,z);

% Calculate 2D cross correlations
b1 = xcorr2fft(a0,a1);
b2 = xcorr2fft(a0,a2);

% Plot the cross correlations
x = -258:257;
y = -774:773;
figure;
subplot(1,2,1);
imagesc(x,y,abs(b1));
axis equal tight;
subplot(1,2,2);
imagesc(x,y,abs(b2));
axis equal tight;

% Generate an RGB image with the 3 components
ii1 = 257;
ii2 = -262;
jj1 = 0;
jj2 = -1;
rgb = cat(3,a0,circshift(a1,[ii1 jj1]),circshift(a2,[ii2 jj2]));
rgb = log10(rgb)/6;

% Plot the RGB image
figure;
image(x,y,rgb);
axis equal tight;
title(num2str([ii1 jj1 ii2 jj2]));

figure;
subplot(1,3,1);
imagesc(x,y,rgb(:,:,1));
axis equal tight;
subplot(1,3,2);
imagesc(x,y,rgb(:,:,2));
axis equal tight;
subplot(1,3,3);
imagesc(x,y,rgb(:,:,3));
axis equal tight;

% Load the data mask
load('MaxiPix_Mask.mat');

% Add zeropadding to the mask
mask = cat(1,z,double(mask),z);

% Convert 0s to NaN
mask(mask == 0) = NaN;

% Generate a mask for each data set
m0 = mask;
m1 = circshift(mask,[ii1 jj1]);
m2 = circshift(mask,[ii2 jj2]);

% Combine the masks to get the final mask
mask = mean(cat(3,m0,m1,m2),3,'omitnan');

% Set NaNs back to 0
mask(isnan(mask)) = 0;

% Plot the total mask
figure;
imagesc(mask);
axis equal tight;

% Combine the intensity measurements
Z = zeros(516,516,281);
I = mean(cat(4,cat(1,Z,A0,Z).*m0,circshift(cat(1,Z,A1,Z),[ii1 jj1 0]).*m1,circshift(cat(1,Z,A2,Z),[ii2 jj2 0]).*m2),4,'omitnan');

% Cut down the size of the data array and mask to optimize for FFT
I = I(255:1278,1:512,1:280);
mask = mask(255:1278,1:512);

% Convert data NaN to 0
I(isnan(I)) = 0;


%% Save the data
% Save the mask and data as doubles
save('NoLens_Stitched.mat','I','mask');



