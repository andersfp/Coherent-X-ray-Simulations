% Initialization
clear;
close all;
clc;


%% Load data
% Load the no lens reconstruction
nl = load('C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\Beam_Time_ID01_20180301\Data_Processing\Grain_G5\No_Lens\Reconstruction_Center_2018_03_13_23_38_300.mat');

% Load the method 1 reconstruction
p1 = load('Reconstruction_Ptycho_1_2018_07_23_14_12_129.mat');

% Load the method 2 reconstruction
p5 = load('Reconstruction_Ptycho_5_2018_07_20_18_49_3809.mat');


%% Calculate the fields
% Get the objects
rn = nl.object_avg;
r1 = p1.rho;
r2 = p5.rho;

% Flip the no-lens reconstruction
rn = flip(flip(flip(conj(rn),1),2),3);

% Set the pixel sizes
delta_rn = 9.5583e-9;
delta_r1 = 1.3809e-8;

% Calculate the number of pixels for no-lens reconstruction
n1 = 256;
nn = delta_r1./delta_rn.*n1;
nn = 2*round(nn./2);

% Zeropad the no-lens reconstruction
rn = padarray(rn,([nn nn 280] - size(rn))/2);

% Calculate the fields
fn = fftshift(fftn(fftshift(rn)));
f1 = fftshift(fftn(fftshift(r1)));
f2 = fftshift(fftn(fftshift(r2)));

% De-pad the no-lens field
[~,ix] = max(max(abs(fn(:,:,141))));
[~,iy] = max(max(abs(fn(:,:,141)).'));
fn = fn(iy-128:iy+127,ix-128:ix+127,:);

% Calculate the intensities
In = abs(fn).^2;
I1 = abs(f1).^2;
I2 = abs(f2).^2;


%% Calculate and apply the mask
% Calculate the mask
thr = 100;
mask = single(In./max(In(:)).*500000 > thr);

% Calculate the truncated objects
r1t = fftshift(ifftn(fftshift(f1.*mask)));
r2t = fftshift(ifftn(fftshift(f2.*mask)));

% Correct the phase slope on the no-lens object
x = (-nn/2):(nn/2 - 1);
y = x.';
dx = nn./2 - ix;
dy = nn./2 - iy;
rn = rn.*exp(1i.*2.*pi.*dx./nn.*x).*exp(1i.*2.*pi.*dy./nn.*y);

% Set the slices to plot
iir = 134;
iif = 141;

% Change the global phase
rn = rn.*exp(-1i.*angle(rn(nn/2+1,nn/2+1,iir)));
r1 = r1.*exp(-1i.*angle(r1(129,129,iir)));
r2 = r2.*exp(-1i.*angle(r2(129,129,iir)));
r1t = r1t.*exp(-1i.*angle(r1t(129,129,iir)));
r2t = r2t.*exp(-1i.*angle(r2t(129,129,iir)));

% Plot all the 3D objects
Slicer(abs(rn),'displayRange',[0 2.5]);
Slicer(abs(r1),'displayRange',[0 0.05]);
Slicer(abs(r2),'displayRange',[0 0.05]);
Slicer(abs(r1t),'displayRange',[0 0.05]);
Slicer(abs(r2t),'displayRange',[0 0.05]);

% Make slice plots
figure;
imagesc(abs(rn(:,:,iir)),[0 2.5]);
axis equal tight;
figure;
imagesc(abs(r1(:,:,iir)),[0 0.05]);
axis equal tight;
figure;
imagesc(abs(r2(:,:,iir)),[0 0.05]);
axis equal tight;
figure;
imagesc(abs(r1t(:,:,iir)),[0 0.05]);
axis equal tight;
figure;
imagesc(abs(r2t(:,:,iir)),[0 0.05]);
axis equal tight;
figure;
imagesc(angle(rn(:,:,iir)),[-pi pi]);
axis equal tight;
figure;
imagesc(angle(r1(:,:,iir)),[-pi pi]);
axis equal tight;
figure;
imagesc(angle(r2(:,:,iir)),[-pi pi]);
axis equal tight;
figure;
imagesc(angle(r1t(:,:,iir)),[-pi pi]);
axis equal tight;
figure;
imagesc(angle(r2t(:,:,iir)),[-pi pi]);
axis equal tight;

% Make reciprocal space plots
figure;
imagesc(log10(abs(fn(:,:,iif))),[2 5]);
axis equal tight;
figure;
imagesc(log10(abs(f1(:,:,iif))),[0 3]);
axis equal tight;
figure;
imagesc(log10(abs(f2(:,:,iif))),[0 3]);
axis equal tight;
figure;
imagesc(mask(:,:,iif),[0 1]);
axis equal tight;
figure;
imagesc(angle(fn(:,:,iif)),[-pi pi]);
axis equal tight;
figure;
imagesc(angle(f1(:,:,iif)),[-pi pi]);
axis equal tight;
figure;
imagesc(angle(f2(:,:,iif)),[-pi pi]);
axis equal tight;



