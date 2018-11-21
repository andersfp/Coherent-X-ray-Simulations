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

% Cut the data
ix0 = 183;
iy0 = 202;
dxy = 256;
I = I(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1,1:280,:);
mask = mask(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1);

% Generate the Gaussian filter
x = (-dxy/2):(dxy/2 - 1);
y = x.';
sigma_det = 25.2;
G = gaussRMS(x,sigma_det).*gaussRMS(y,sigma_det);
G = single(G);

% Save the data
save('Ptycho_Data.mat','I','mask','G');


