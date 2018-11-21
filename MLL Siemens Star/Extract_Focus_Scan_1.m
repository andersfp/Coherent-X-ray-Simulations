% Initialization
clear;
close all;
clc;


%% Data location
% Set the directory
p = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_2D_17.3keV\pco_focus_30s\';

% Set the base name
fp = '*.edf';

% Make the file lists
[~,f] = fileNames(p,fp);


%% Load the data
% Get number of images
n = length(f);

% Pre-allocate info structure
info(n) = edf_info(f{n});

% Set the image dimensions
mx = info(n).dim_1;
my = info(n).dim_2;

% Read the images
img = zeros(my,mx,n);
fprintf([num2str(n) ': ']);
for i = 1:n
    img(:,:,i) = edf_read(f{i});
    info(i) = edf_info(f{i});
    fprintf('.');
end
fprintf('\n');


%% Process the data
% Plot the raw data
Slicer(img,'displayRange',[0 1200]);

% Load the background
load('Background.mat');

% Load the flatfield
load('Flatfield.mat');

% Subtract the background from the data
A = abs(img - bg);

% Find the center of mass of each image
x = 1:mx;
y = (1:my).';
cx = squeeze(sum(sum(A.*x))./sum(sum(A)));
cy = squeeze(sum(sum(A.*y))./sum(sum(A)));

% Shift the images to align the lens FOV to the flatfield
i0 = 31;
sx = round(cx - cx(i0));
sy = round(cy - cy(i0));
B = A;
for i = 1:n
    B(:,:,i) = circshift(B(:,:,i),-sx(i),2);
    B(:,:,i) = circshift(B(:,:,i),-sy(i),1);
    fprintf('.');
end
fprintf('\n');

% Perform the flatfield correction
I = abs(B./(ff2 - bg));

% Remove NaNs and Infs
I(isnan(I)) = 0;
I(isinf(I)) = 0;


%% Shift the data to align the star in each slice
% Select a small ROI without noise
C = I(745:1256,451:962,:);

% Calculate the slice-by-slice cross-correlation
D = zeros(size(C)-[0 0 1]);
for i = 2:n
    D(:,:,i-1) = xcorr2fft(C(:,:,i-1),C(:,:,i));
end

% Determine the maximum of each cross-correlation
[~,ix] = max(max(D,[],1),[],2);
[~,iy] = max(max(D,[],2),[],1);

% Find the amount that should be switched relative to the first slice
ix = [0;cumsum(squeeze(ix) - 257)];
iy = [0;cumsum(squeeze(iy) - 257)];

% Shift the slices
i0x = 225;
i0y = 30;
for i = 1:n
    I(:,:,i) = circshift(I(:,:,i),i0x-ix(i),2);
    I(:,:,i) = circshift(I(:,:,i),i0y-iy(i),1);
    fprintf('.');
end
fprintf('\n');

% Plot the corrected images
Slicer(I,'displayRange',[0.7 1.3]);


%% Save the data
% Save the corrected image stack
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';
save([dp 'Focus_Scan_1.mat'],'I','info','-v7.3');


