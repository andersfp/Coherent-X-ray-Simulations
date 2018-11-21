% Initialization
clear;
close all;
clc;


%% Set parameters and load data
% Set the resolution
res = 100e-9;

% Load the simulated intensities
f = 'MLL_Displace_Aligned_Lenses.mat';
A = load(f);
I = A.I;
d = A.d;
xd = A.xd;

% Get the image size
m = length(xd);
n = size(I,3);

% Renormalize the images
for i = 1:n
    [N,bin] = histcounts(I(:,:,i),10001);
    N(1) = NaN;
    [~,ii] = max(N);
    I(:,:,i) = I(:,:,i)/bin(ii);
end

% Plot the images
ii = [1:5 n:-1:5];
Ip = I(:,:,ii);
dp = d(ii);
maxval = max(Ip(:));
figure;
for i = 1:10
    subplot(2,5,i);
    imagesc(1e3*xd,1e3*xd,Ip(:,:,i),[0 maxval]);
    axis equal tight;
    set(gca,'YDir','normal');
    title(['Distance: ' num2str(1e3*dp(i),'%+i') ' mm']);
    xlabel('x (mm)');
    ylabel('y (mm)');
end


%% Convolute with resolution function
% Calculate magnification
f = 38e-3;
M = (2 - f)/f;

% Get the resolution width in the detector plane
mres = M*res;

% Generate full coordinate grid
[X,Y] = meshgrid(xd,xd);

% Calculate the point spread function
psf = gauss(X,mres).*gauss(Y,mres)*mean(diff(xd))^2;

% Perform the convolution
Ic = zeros(size(I));
for i = 1:n
    Ic(:,:,i) = fconv(I(:,:,i),psf);
end

% Plot the images again
Ipc = Ic(:,:,ii);
figure;
for i = 1:10
    subplot(2,5,i);
    imagesc(1e3*xd,1e3*xd,Ipc(:,:,i),[0 maxval]);
    axis equal tight;
    set(gca,'YDir','normal');
    title(['Distance: ' num2str(1e3*dp(i),'%+i') ' mm']);
    xlabel('x (mm)');
    ylabel('y (mm)');
end


%% Make line scans
% Set the coordinates
% p1 = [-0.3 0]*1e-3;
% p2 = [0 0.3]*1e-3;
p1 = [-0.20 0.28]*1e-3;
p2 = [-0.31 0.16]*1e-3;

% Make a plot showing the position of the line profile
figure;
imagesc(1e3*xd,1e3*xd,Ic(:,:,5),[0 maxval]);
line(1e3*[p1(1) p2(1)],1e3*[p1(2) p2(2)],'Color','black','LineWidth',2);
axis equal tight;
set(gca,'YDir','normal');
title('Line profile');
xlabel('x (mm)');
ylabel('y (mm)');

% Calculate the coordinate index
[~,ix1] = min(abs(xd - p1(1)));
[~,iy1] = min(abs(xd - p1(2)));
[~,ix2] = min(abs(xd - p2(1)));
[~,iy2] = min(abs(xd - p2(2)));

% Get number of pixels in resulting line
nl = ceil(sqrt((ix1 - ix2)^2 + (iy1 - iy2)^2));

% Calculate coordinate length for the line
xlmax = sqrt(dot(p2 - p1,p2 - p1));
xl = linspace(0,xlmax,nl)';

% Get the line profiles
l = zeros(nl,n);
for i = 1:n
    l(:,i) = improfile(Ic(:,:,i),[ix1 ix2],[iy1 iy2],nl,'bilinear');
end

% Set plotting parameters
ylim = [0.5 1.5];

% Plot the line profiles of positive shifted lens
figure;
plot(1e3*xl,l(:,5:end));
set(gca,'YLim',ylim);
title('Positively shifted line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('0 mm','+1 mm','+2 mm','+3 mm','+4 mm');

% Plot the line profiles of negative shifted lens
figure;
plot(1e3*xl,l(:,5:-1:1));
set(gca,'YLim',ylim);
title('Negtively shifted line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('0 mm','-1 mm','-2 mm','-3 mm','-4 mm');

% Plot the line profiles of 1 mm pair
figure;
plot(1e3*xl,l(:,[4 6]));
set(gca,'YLim',ylim);
title('Pair of 1 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-1 mm','+1 mm');

% Plot the line profiles of 2 mm pair
figure;
plot(1e3*xl,l(:,[3 7]));
set(gca,'YLim',ylim);
title('Pair of 2 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-2 mm','+2 mm');

% Plot the line profiles of 3 mm pair
figure;
plot(1e3*xl,l(:,[2 8]));
set(gca,'YLim',ylim);
title('Pair of 3 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-3 mm','+3 mm');

% Plot the line profiles of 4 mm pair
figure;
plot(1e3*xl,l(:,[1 9]));
set(gca,'YLim',ylim);
title('Pair of 4 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-4 mm','+4 mm');

% Save the line profiles for further analysis
%save('Simulated_Line_Profiles.mat','xl','l','d','Ic');



