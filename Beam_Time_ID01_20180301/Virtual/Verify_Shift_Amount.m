% Initialization
clear;
close all;
clc;


%% Load the data
% Load the .mat file
load('Virtual_Shifted.mat');

% Perform intensity correction
C = [1.0000 2.3865 1.3163 1.5545 5.3719 3.5688 4.2160 14.4266 7.6019];
C = permute(C,[3 4 1 2]);
I = C.*I;

% Shift the data along the rocking direction
S = [-1 -1 -2 -2 -2 -3 -3 -3 -4];
n = length(S);
for i = 1:n
    I(:,:,:,i) = circshift(I(:,:,:,i),S(i),3);
end

% Plot the central slices
Slicer(squeeze(log10(I(:,:,141,:))),'displayRange',[0 6]);

% Plot the mask
Slicer(squeeze(mask),'displayRange',[0 1]);


%% Process the data
% Check the total intensity
figure;
plot(squeeze(sum(sum(sum(I)))));
figure;
plot(squeeze(sum(sum(I(:,:,141,:)))));

% Plot the images side by side
i0 = 141;
cl = [0 6];
ii = [5 2 8 4 1 7 6 3 9];
figure;
set(gcf,'Position',[2200 42 1074 1074]);
for i = 1:n
    subplot(3,3,i);
    imagesc(log10(I(:,:,i0,ii(i))),cl);
    axis equal tight off;
    x0 = mod(i-1,3)/3;
    y0 = (2 - floor((i-1)/3))/3;
    set(gca,'Position',[x0 y0 1/3 1/3]);
end

% Check the intensity
i1 = 125;
i2 = 133;
Ii = squeeze(sum(sum(I(i1:i2,i1:i2,:,:)./mask(i1:i2,i1:i2,:,:))));
figure;
plot(Ii);

% Calculate corrections
A = sort(Ii,'descend');
c = mean(A(1:10,:));
c = c(1)./c;

% Calculate shifts
[~,s] = max(Ii);
s = 141 - s;

% Compute cross-correlations
A = squeeze(I(:,:,i0,:)./sum(I,3));
A(isnan(A)) = 0;
A(isinf(A)) = 0;
A = A(109:147,109:147,:);
B = zeros(size(A));
for i = 1:n
    B(:,:,i) = xcorr2fft(A(:,:,1),A(:,:,i));
end
B = B./max(max(B));

Slicer(A);
Slicer(B);


