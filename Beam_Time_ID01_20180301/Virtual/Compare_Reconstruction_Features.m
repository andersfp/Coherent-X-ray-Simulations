% Initialization
clear;
close all;
clc;


%% Load data
% Load the no lens data
%dn = load('C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\Beam_Time_ID01_20180301\Data_Processing\Grain_G5\No_Lens\Corrected_NoLens_Center.mat');

% Load the virtual geometry center data
dv = load('Corrected_Virtual_Center.mat');

% Load the virtual geometry ptycho data
dp = load('Corrected_Virtual_Ptycho_1.mat');


%% Extract the two objects
% Get the virtual image reconstruction
V = double(dv.object_avg/100);

% Get the ptycho image reconstruction
P = double(dp.rho);

% Cut the size of the ptycho reconstruction
x0 = 109 - 29; % 109 - 29
xw = size(V,2);
y0 = 119 - 23; % 120 - 23
yw = size(V,1);
z0 = 122 - 25; % 123 - 25
zw = size(V,3);
P = P(y0 + 1:y0 + yw,x0 + 1:x0 + xw,z0 + 1:z0 + zw);

% Compare the xy-plane and zx-plane
figure;
subplot(1,2,1);
imshowpair(sum(abs(V),3),sum(abs(P),3),'falsecolor');
xlabel('x');
ylabel('y');
subplot(1,2,2);
imshowpair(squeeze(sum(abs(V),1)),squeeze(sum(abs(P),1)),'falsecolor');
xlabel('z');
ylabel('x');


%% Compare 2D slices
% Set slice number
i0 = 38;

% Make the 2D slices
v = V(:,:,i0);
p = P(:,:,i0);

% Normalize the amplitude of the two slices
v = v./max(abs(v(:)));
p = p./max(abs(p(:)));

% Plot the two slices with phase and amplitude
cmap = hsv(256);
rgbv = complex2rgb(v,cmap);
rgbp = complex2rgb(p,cmap);
figure;
subplot(1,2,1);
image(rgbv);
axis equal tight;
subplot(1,2,2);
image(rgbp);
axis equal tight;

% Plot the amplitude comparison
figure;
imshowpair(abs(v),abs(p),'falsecolor');

% Plot the amplitude difference
figure;
imagesc(abs(p) - abs(v));
axis equal tight;

% Plot the phase comparison
figure;
imshowpair(angle(v),angle(p),'falsecolor');

% Plot the phase difference
figure;
imagesc(angle(p) - angle(v));
axis equal tight;



