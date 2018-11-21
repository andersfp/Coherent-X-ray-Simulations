% Initialization
clear;
close all;
clc;


%% Load data
% Load the shifted data
load('Processed_Virtual_Shifted_2.mat');

% Extract the central images
i0 = size(I,3)/2 + 1;
im = squeeze(I(:,:,i0,:));
%im = log10(squeeze(I(:,:,i0,:)));
%im(isinf(im)) = 0;

% Sort the images
im2 = reshape(im,256,256,3,3);
im2 = im2(:,:,[3 1 2],[2 1 3]);


%% Image registration
% Try Hugh's function
data.intensity = im2;
scale = 1;
n = 50;
dyx1 = zeros(2,2,n);
for i = 1:n
    dyx0 = rand(2)*2 - 1;
    [dyx1(:,:,i)] = dxm2d_findRegistration(data,scale,dyx0);
end


figure;
plot(1:n,squeeze([dyx1(1,1,:) dyx1(2,1,:) dyx1(1,2,:) dyx1(2,2,:)]));

