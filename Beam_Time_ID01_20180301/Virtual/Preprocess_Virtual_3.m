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
I = double(I);

% Load the mask
load('MaxiPix_Mask.mat');
mask = double(mask);

% Add additional masking
mask2 = imread('Mask_Lens.png');
mask = double(mask > 0 & mask2 > 127);


%% Process the data
% Pad the data with zeros for shifting the data
pd = 0;
I = padarray(I,[pd pd 0 0]);

% Generate matching masks
mask = padarray(mask,[pd pd]);
mask = repmat(mask,1,1,1,n);

% Generate coordinates
x = ((-nx/2):(nx/2-1));
y = ((-ny/2):(ny/2-1)).';
z = permute(((-(n_omega - 1)/2):((n_omega - 1)/2)),[3 1 2]);

% Generate an effective pupil
sigma_det = 24.2;
G = gaussRMS(x - (183 - nx/2),sigma_det).*gaussRMS(y - (202 - ny/2),sigma_det);
G = padarray(G,[pd pd]);

% Generate full size pupil
G = repmat(G,1,1,1,n);


%% Shift the data in the rocking direction
% Extract the rocking curves
A = I./G.*(G > 0.05);
A = squeeze(sum(sum(A,1,'omitnan'),2,'omitnan'));

% Find the 10 highest values of the rocking curves
[~,ii] = sort(A,1,'descend');

% Calculate COM
z0 = zeros(n,1);
w = z0;
for i = 1:n
    w(i) = sum(A(ii(1:10,i),i));
    z0(i) = sum(squeeze(z(ii(1:10,i))).*A(ii(1:10,i),i))./w(i);
end

% Shift and normalize the data
zs = z0 - min(z0) - 1;
for i = 1:n
    I(:,:,:,i) = imshift3(I(:,:,:,i),[0 0 -zs(i)])./w(i);
end


%% Make a rough shift in the xy direction
% Set the shifting amount
sx = [0 0 0 -27 -27 -27 27 27 27]; % 29
sy = [0 -35 35 0 -35 35 0 -35 35]; % 38

% Shift the data
for i = 1:n
    I(:,:,:,i) = circshift(I(:,:,:,i),[sy(i) sx(i) 0]);
    mask(:,:,:,i) = circshift(mask(:,:,:,i),[sy(i) sx(i) 0]);
    G(:,:,:,i) = circshift(G(:,:,:,i),[sy(i) sx(i) 0]);
end


%% Make a fine adjustment of the xy shift
% Extract the central slices
A = squeeze(I(:,:,143,:));

% Cut out an ROI
A = A(208-64+1:208+64,184-64+1:184+64,:);

% Sort the images
A = reshape(A,128,128,3,3);
A = A(:,:,[3 1 2],[2 1 3]);

% Run the registration
data.intensity = log10(A);
scale = 1;
m = 100;
dyx1 = zeros(2,2,m);
for i = 1:m
    dyx0 = rand(2)*2 - 1;
    [dyx1(:,:,i)] = dxm2d_findRegistration(data,scale,dyx0);
end

% Determine the shift
dyx = median(dyx1,3);

% Calculate fine shift amounts
u = [1 1 1;2 2 2;3 3 3] - 3/2;
v = [1 2 3;1 2 3;1 2 3] - 3/2;
fsx = dyx(1,2)*u + dyx(2,2)*v;
fsy = dyx(1,1)*u + dyx(2,1)*v;

% Reshape the shift amounts
fsx = reshape(fsx([2 3 1],[2 1 3]),9,1);
fsy = reshape(fsy([2 3 1],[2 1 3]),9,1);

% Shift the intensity
for i = 1:n
    I(:,:,:,i) = imshift3(I(:,:,:,i),[-fsx(i) -fsy(i) 0]);
    mask(:,:,:,i) = imshift2(mask(:,:,:,i),[-fsx(i) -fsy(i)]);
    G(:,:,:,i) = imshift2(G(:,:,:,i),[-fsx(i) -fsy(i)]);
end


%% Save the data
% Cut out the relevant part
ix0 = 183;
iy0 = 202;
dxy = 256;
I = I(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1,5:280,:);
mask = mask(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1,:,:);
G = G(iy0 - dxy/2:iy0 + dxy/2 - 1,ix0 - dxy/2:ix0 + dxy/2 - 1,:,:);

% Scale the intensity and convert to single
I = I./max(I(:)).*1e6;
I = single(I);
mask = single(mask);
G = single(G);

% Check the overlap in the center
A = log10(squeeze(I(:,:,143,1:3)))/6;
B = log10(squeeze(I(:,:,143,4:6)))/6;
C = log10(squeeze(I(:,:,143,7:9)))/6;
figure;
subplot(1,3,1);
image(A);
axis equal tight off;
set(gca,'Position',[0 0 1/3 1]);
subplot(1,3,2);
image(B);
axis equal tight off;
set(gca,'Position',[1/3 0 1/3 1]);
subplot(1,3,3);
image(C);
axis equal tight off;
set(gca,'Position',[2/3 0 1/3 1]);

% Check the overlap in the center
A = log10(squeeze(I(:,:,143,1:3:9)))/6;
B = log10(squeeze(I(:,:,143,2:3:9)))/6;
C = log10(squeeze(I(:,:,143,3:3:9)))/6;
figure;
subplot(1,3,1);
image(A);
axis equal tight off;
set(gca,'Position',[0 0 1/3 1]);
subplot(1,3,2);
image(B);
axis equal tight off;
set(gca,'Position',[1/3 0 1/3 1]);
subplot(1,3,3);
image(C);
axis equal tight off;
set(gca,'Position',[2/3 0 1/3 1]);
return
% Save the shifted data
save('Processed_Virtual_Shifted_3.mat','I','mask','G');

% Stitch the data
I = sum(I,4);
mask = sum(mask,4);
mask = mask./max(mask(:));
G = sum(G,4);
G = G./max(G(:));

% Save the stitched data
save('Processed_Virtual_Stitched_3.mat','I','mask','G');


