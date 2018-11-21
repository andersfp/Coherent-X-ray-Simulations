% Initialization
clear;
close all;
clc;


%% Load data
% Load the diffraction patterns
tic;
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';
I = h5read([p 'Scan0039.cxi'],'/entry_1/instrument_1/detector_1/data');
load('MaxiPix_Mask.mat');
toc;

% Convert data
I = permute(I,[2 1 3]);
I = single(I);
mask = single(mask);
%I = double(I);
%mask = double(mask);

% Add additional masking
mask2 = imread('Mask_Lens.png');
mask = single(mask > 0 & mask2 > 127);
%mask = double(mask > 0 & mask2 > 127);

% Determine center of mass
x = 1:516;
y = x.';
cx = sum(sum(I(:,:,141).*x))./sum(sum(I(:,:,141)));
cy = sum(sum(I(:,:,141).*y))./sum(sum(I(:,:,141)));
ix = round(cx);
iy = round(cy);

% Crop the data
w = 128;
I = I(iy-w:iy+w-1,ix-w:ix+w-1,1:280);
mask = mask(iy-w:iy+w-1,ix-w:ix+w-1);


%% Reconstruction
% Use GPU
gpu = 1;

% Set up reconstruction parameters
iterER = 1;
iterHIO = 20;
cyc = 100;
beta = 0.9;
tol = 0.10;
sig0 = 2.5;
sigmin = 1.7;
sinit = 0.0005;
sig = sig0;
ER = 4;
HIO = 3;

% Coordinates
nx = size(I,2);
ny = size(I,1);
n_omega = size(I,3);
x = ((-nx/2):(nx/2-1));
y = ((-ny/2):(ny/2-1)).';
z = permute(((-n_omega/2):(n_omega/2-1)),[3 1 2]);
x = single(x);
y = single(y);
z = single(z);

% Perform the reconstruction
data_input = sqrt(I);
support = ifftshift(fftn(fftshift(data_input)));
%support = abs(support) > sinit*max(abs(support(:)));
support = single(abs(support) > sinit*max(abs(support(:))));
%support = double(abs(support) > sinit*max(abs(support(:))));
field = data_input.*exp(1i*2*pi*rand(ny,nx,n_omega));
object = ifftshift(ifftn(fftshift(field)));
chi_square_real = [];
chi_square_reci = [];

% Make the mask variable intensity
sigma_det = 24.2;
G = sqrt(gaussRMS(x,sigma_det).*gaussRMS(y,sigma_det));
data_input = data_input./G;
mask = mask.*(G > 0.4);

% Calculate the mean value of unmasked data next to masked pixels
if isempty(gcp('nocreate'))
    parpool;
end
tic;
[h,dm,ds] = maskNeighbor2(mask,data_input,100,0.1);
%[h,dm,ds] = maskNeighbor3(mask,data_input,100);
toc;

% Set the amplitude limit
lim = dm + 0.5*ds;

% Group the masked pixel indices, mean and deviation of the neighbors
epsilon.h = h;
epsilon.data_mean = dm;
epsilon.lim = lim;
%epsilon.lim = 0;


%% Start the reconstruction
% Send to GPU if used
if gpu == 1
    data_input = gpuArray(data_input);
    support = gpuArray(support);
    field = gpuArray(field);
    object = gpuArray(object);
    mask = gpuArray(mask);
    x = gpuArray(x);
    y = gpuArray(y);
    z = gpuArray(z);
    epsilon = structfun(@gpuArray,epsilon,'UniformOutput',false);
end

% Perform the reconstruction
tic;
for i = 1:cyc
    % HIO
    [field,object,chi_square_real,chi_square_reci] = cycle(HIO,iterHIO,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    % ER
    [field,object,chi_square_real,chi_square_reci] = cycle(ER,iterER,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    % Shrinkwrap
    [support,sig] = shrinkwrap(object,x,y,z,sig,sigmin,tol,support);
    % Plots
    figure(10);
    semilogy(chi_square_reci);
    figure(100);
    subplot(1,2,1);
    imagesc(x,y,support(:,:,n_omega/2+1));
    axis equal tight;
    subplot(1,2,2);
    imagesc(x,y,abs(object(:,:,n_omega/2+1)));
    axis equal tight;
    drawnow;
end
toc;

% Obtain averaging of reconstructions
%object_avg = gather(object);
object_avg = object;
tic;
for i = 1:cyc
    % HIO
    [field,object,chi_square_real,chi_square_reci] = cycle(HIO,iterHIO,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    % ER
    [field,object,chi_square_real,chi_square_reci] = cycle(ER,iterER,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    % Averaging
    %object_avg = object_avg + gather(object);
    object_avg = object_avg + object;
    % Plots
    figure(10);
    semilogy(chi_square_reci);
    figure(100);
    subplot(1,2,2);
    imagesc(x,y,abs(object(:,:,n_omega/2+1)));
    axis equal tight;
    drawnow;
end
toc;

% Collect data from GPU if used
if gpu == 1
    data_input = gather(data_input);
    support = gather(support);
    field = gather(field);
    object = gather(object);
    mask = gather(mask);
    chi_square_real = gather(chi_square_real);
    chi_square_reci = gather(chi_square_reci);
    object_avg = gather(object_avg);
    x = gather(x);
    y = gather(y);
    z = gather(z);
    epsilon = structfun(@gather,epsilon,'UniformOutput',false);
end

% Save the error metric figure for quick reconstruction estimation
figure(10);
saveas(gcf,['Error_Metric_Center_' timestr(3) '_' num2str(round(chi_square_reci(end)),'%.0f') '.png']);

% Plot the support
figure(20);
imagesc(support(:,:,n_omega/2+1));
axis equal tight;
title('Support');

% Plot the reciprocal space reconstruction
figure(40);
subplot(1,2,1);
imagesc(log10(abs(field(:,:,n_omega/2+1))));
axis equal tight;
title('Reciprocal space amplitude');
subplot(1,2,2);
imagesc(angle(field(:,:,n_omega/2+1)));
axis equal tight;
title('Reciprocal space phase');

% Plot the reconstruction
figure(50);
subplot(1,2,1);
imagesc(abs(object(:,:,n_omega/2+1)));
axis equal tight;
title('Object amplitude');
subplot(1,2,2);
imagesc(angle(object(:,:,n_omega/2+1)));
axis equal tight;
title('Object phase');

% Make 3D plot of the raw reconstructed object
figure;
pp = patch(isosurface(abs(object),0.1*max(abs(object(:)))));
isonormals(abs(object),pp);
pp.FaceColor = 'red';
pp.EdgeColor = 'none';
daspect([1 1 1]);
view(3);
camlight(45,45);
camlight(-135,-45);
lighting gouraud;
xlabel('y');
ylabel('z');
zlabel('x');

% Crop the output object
xs = x.*support;
ys = y.*support;
zs = z.*support;
xmin = min(xs(:));
xmax = max(xs(:));
ymin = min(ys(:));
ymax = max(ys(:));
zmin = min(zs(:));
zmax = max(zs(:));
x0 = round((xmax + xmin)/2) + nx/2;
xw = round(xmax - xmin);
y0 = round((ymax + ymin)/2) + ny/2;
yw = round(ymax - ymin);
z0 = round((zmax + zmin)/2) + n_omega/2;
zw = round(zmax - zmin);
object = object(y0 - yw:y0 + yw - 1,x0 - xw:x0 + xw - 1,z0 - zw:z0 + zw - 1);
object_avg = object_avg(y0 - yw:y0 + yw - 1,x0 - xw:x0 + xw - 1,z0 - zw:z0 + zw - 1);

% Save the resulting objects
save(['Reconstruction_Center_2_' timestr(3) '_' num2str(round(chi_square_reci(end)),'%.0f') '.mat'],...
    'object','object_avg','chi_square_real','chi_square_reci','ER','HIO','iterER','iterHIO','cyc','beta','tol','sig0','sigmin','sig','sinit','nx','ny','n_omega');

