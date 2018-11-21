% Initialization
clear;
close all;
clc;


%% Load data
% Load the diffraction patterns
fn = 39:47;
n = length(fn);
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';
nx = 516;
ny = 516;
n_omega = 281;
I = zeros(ny,nx,n_omega,n);
for i = 1:n
    I(:,:,:,i) = h5read([p 'Scan00' num2str(fn(i)) '.cxi'],'/entry_1/instrument_1/detector_1/data');
end

% Convert data
I = permute(I,[2 1 3 4]);
I = single(I);

% Load the mask
load('MaxiPix_Mask.mat');
mask = single(mask);

% Add additional masking
mask2 = imread('Mask_Lens.png');
mask = single(mask > 0 & mask2 > 127);


%% Process the data
% Pad the data with zeros for shifting the data
pd = 200;
I = padarray(I,[pd pd 0 0]);

% Generate matching masks
mask = padarray(mask,[pd pd]);
mask = repmat(mask,1,1,1,n);

% Generate coordinates
x = ((-nx/2):(nx/2-1));
y = ((-ny/2):(ny/2-1)).';
z = permute(((-(n_omega - 1)/2):((n_omega - 1)/2)),[3 1 2]);
x = single(x);
y = single(y);
z = single(z);

% Generate an effective pupil
sigma_det = 24.2;
G = gaussRMS(x - (383 - pd - nx/2),sigma_det).*gaussRMS(y - (402 - pd - ny/2),sigma_det);
G = padarray(G,[pd pd]);

% % Remove the effect of the pupil
% Is = I./G;
% Is(isnan(Is)) = 0;
% Is(isinf(Is)) = 0;

% Generate full size pupil
G = repmat(G,1,1,1,n);


%% Shift the data
%I2 = Is;

% Set the shifting amount
sx = [0 0 0 -28 -28 -28 28 28 28];
sy = [0 -38 38 0 -38 38 0 -38 38];% + [0 0 0 -2 -2 -2 2 2 2];

% Shift the data
for i = 1:n
    %I2(:,:,:,i) = circshift(I2(:,:,:,i),[sy(i) sx(i) 0]);
    I(:,:,:,i) = circshift(I(:,:,:,i),[sy(i) sx(i) 0]);
    mask(:,:,:,i) = circshift(mask(:,:,:,i),[sy(i) sx(i) 0]);
    G(:,:,:,i) = circshift(G(:,:,:,i),[sy(i) sx(i) 0]);
end

% % Check the overlap in the center
% A = log10(squeeze(I2(:,:,141,1:3)))/6;
% B = log10(squeeze(I2(:,:,141,4:6)))/6;
% C = log10(squeeze(I2(:,:,141,7:9)))/6;
% figure;
% subplot(1,3,1);
% image(A);
% axis equal tight;
% subplot(1,3,2);
% image(B);
% axis equal tight;
% subplot(1,3,3);
% image(C);
% axis equal tight;
% 
% % Check the overlap in the center
% A = log10(squeeze(I2(:,:,141,1:3:9)))/6;
% B = log10(squeeze(I2(:,:,141,2:3:9)))/6;
% C = log10(squeeze(I2(:,:,141,3:3:9)))/6;
% figure;
% subplot(1,3,1);
% image(A);
% axis equal tight;
% subplot(1,3,2);
% image(B);
% axis equal tight;
% subplot(1,3,3);
% image(C);
% axis equal tight;

% Stitch the mask
mask = mask.*G;

% Cut out the relevant part
ix0 = 384;
iy0 = 403;
dxy = 256;
I = I(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1,1:280,:);
mask = mask(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1,:,:);

% Save the data
save('Virtual_Shifted.mat','I','mask');


