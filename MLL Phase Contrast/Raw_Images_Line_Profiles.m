% Initialization
clear;
close all;
clc;


%% Load the images
% Get a list of the files
p = [pwd '\phase_data_MLL\'];
lst = dir([p '*.png']);
lst = struct2cell(lst);
lst = lst(1,:)';

% Load the files
n = length(lst);
nx = 775;
ny = 800;
img = zeros(ny,nx,n);
for i = 1:n
    A = double(imread([p lst{i}]));
    img(:,:,i) = sqrt(A(:,:,1).^2 + A(:,:,2).^2 + A(:,:,3).^2);
end

% Get the lens motor position
d = cellfun(@(x) str2double(x(10:15)),lst);


%% Make line profiles
% Establish moving zero-point
x0 = linspace(298,209,n)';
%y0 = linspace(247,281,n)';
y0 = linspace(247,281,n)';

% Make full coordinates
x = (1:nx)';
y = (1:ny)';
[X,Y] = meshgrid(x,y);

% Plot first image
figure;
imagesc(x-x0(1),y-y0(1),img(:,:,1),[150 320]);
axis equal tight;
title('First image');

% Make direction vectors for query grid
dxq = [137 + 118;76 + 50];
dxq = dxq/sqrt(dot(dxq,dxq));
dyq = [-dxq(2);dxq(1)];

% Set origin for the query grid
%q0 = [119;-20];
q0 = [100;-20];

% Make query grid
%nxq = 40;
nxq = 80;
nyq = 150;
x0q = (0:(nxq - 1))';
y0q = (0:(nyq - 1))';
[X0q,Y0q] = meshgrid(x0q,y0q);
Xq = dxq(1)*X0q + dyq(1)*Y0q + q0(1);
Yq = dxq(2)*X0q + dyq(2)*Y0q + q0(2);

% Interpolate the area
L = zeros(nyq,nxq,n);
for i = 1:n
    L(:,:,i) = interp2(X-x0(i),Y-y0(i),img(:,:,i),Xq,Yq);
end

% Plot the first image section
figure;
imagesc(x0q,y0q,L(:,:,i));
axis equal tight;
title('Line profiles to be averaged');

% Make the line profiles
l = squeeze(mean(L,2));

% Plot selected line profiles
xl2 = y0q*1e-6;
figure;
plot(1e3*xl2,l(:,2:10:end));
title('Line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');

% Plot all the line profiles
ii = (2:1:n)';
m = length(ii);
ofs = 15*(0:(m-1));
m2 = round(m/2);
cmap1 = [linspace(1,0,m2)' zeros(m2,2)];
cmap2 = [zeros(m2,2) linspace(0,1,m2)'];
if m2 == m/2
    cmap = [cmap1;cmap2];
else
    cmap = [cmap1(1:end-1,:);cmap2];
end
figure;
set(gca,'ColorOrder',cmap,'NextPlot','replacechildren');
plot(1e3*xl2,l(:,ii) + ofs);
title('Line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');

% Extract comparable line profiles
i0 = 37;
d2 = d - d(i0);
[~,i1] = min(abs(d2 + 4));
[~,i2] = min(abs(d2 + 3));
[~,i3] = min(abs(d2 + 2));
[~,i4] = min(abs(d2 + 1));
[~,i5] = min(abs(d2));
[~,i6] = min(abs(d2 - 1));
[~,i7] = min(abs(d2 - 2));
ii = [i1 i2 i3 i4 i5 i6 i7];
% img2 = img(:,:,ii);
% L2 = L(:,:,ii);
% l2 = l(:,ii);
% d2 = d2(ii);
img2 = (img(:,:,ii) + img(:,:,ii+1) + img(:,:,ii-1) + img(:,:,ii+2) + img(:,:,ii-2))/5;
L2 = (L(:,:,ii) + L(:,:,ii+1) + L(:,:,ii-1) + L(:,:,ii+2) + L(:,:,ii-2))/5;
l2 = (l(:,ii) + l(:,ii+1) + l(:,ii-1) + l(:,ii+2) + l(:,ii-2))/5;
d2 = (d2(ii) + d2(ii+1) + d2(ii-1) + d2(ii+2) + d2(ii-2))/5;

% k = mean(l2([1:10 nyq-10:nyq],:));
% %k = mean(k);
% l2 = l2./k;

% Plot all 7 selected line profiles
figure;
plot(1e3*xl2,l2);
title('Selected line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-4 mm','-3 mm','-2 mm','-1 mm','0 mm','+1 mm','+2 mm');

% Plot the 1 mm pair
figure;
plot(1e3*xl2,l2(:,[4 6]));
title('Pair of 1 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-1 mm','+1 mm');

% Plot the 2 mm pair
figure;
plot(1e3*xl2,l2(:,[3 7]));
title('Pair of 2 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-2 mm','+2 mm');

% Save the measured line profiles
%save('Measured_Line_Profiles_2.mat','xl2','l2','d2','img2');


