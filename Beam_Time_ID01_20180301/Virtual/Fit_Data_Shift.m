% Initialization
clear;
close all;
clc;


%% Load data
% Load the diffraction patterns
load('Ptycho_Data.mat');

% Load an initial guess
dat = load('Reconstruction_Center_2018_04_17_16_58_1161.mat');
rho = dat.object_avg./dat.cyc;

% Set the correct data size
nx = size(I,2);
ny = size(I,1);
n_omega = size(I,3);
n = size(I,4);
rho = padarray(rho,([ny nx n_omega] - size(rho))/2,0);

% Convert to double precision
I = double(I);
rho = double(rho);
mask = double(mask);
G = double(G);


%% Fit the displacements
% Generate pixel coordinates
x = (-nx/2):(nx/2 - 1);
y = ((-ny/2):(ny/2 - 1)).';
z = permute((-n_omega/2):(n_omega/2 - 1),[3 1 2]);

% Set initial guess of xy-position alignment of diffraction data
bx0 = [0 0 0 -1 -1 -1 1 1 1].'.*29; % 29
by0 = [0 -1 1 0 -1 1 0 -1 1].'.*38 + [0 0 0 -1 -1 -1 1 1 1].'.*2; % 38, 2
% bx = single(permute(bx0,[2 3 4 1]));
% by = single(permute(by0,[2 3 4 1]));

% Calculate guess for z-position alignment of diffraction data
A = sum(sum(I));
bz = -sum(A.*z)./sum(A);
bz = bz - bz(1);
bz0 = squeeze(bz);

% Send data to GPU
rho = gpuArray(rho);
I = gpuArray(I);
mask = gpuArray(mask);
G = gpuArray(G);

% Optimization options
opt = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1e10,'TolX',1e-3);

% Initial guess
%x0 = cat(1,bx0.',by0.',bz0.',squeeze(sum(sum(sum(gather(I))))).'./2.6349e8);
% x0 = [    0.0006    0.0034    0.0013  -28.5979  -28.8212  -28.7341   27.0028   26.9619   27.3943
%     0.0009  -38.1134   32.4225   -3.9212  -40.6277   30.9473   -4.6579  -41.8149   30.8222
%     0.0004   -0.2720    0.2061   -0.1339   -0.8718   -0.0756   -0.8533   -1.3481   -1.5224
%     1.3176    0.6479    1.0804    0.9225    0.3223    0.5287    0.3259    0.1272    0.2403];
x0 = [    0.0006   -0.2205   -0.2921  -28.5979  -28.8216  -28.7341   27.0028   26.9619   27.3943
    0.0008  -38.0828   32.4441   -3.9212  -40.6277   30.9477   -4.6579  -41.8149   30.8222
    0.0004   -0.2330    0.2651   -0.1339   -0.8720   -0.0762   -0.8533   -1.3481   -1.5224
    1.3184    0.6457    1.0792    0.9225    0.3223    0.5288    0.3259    0.1272    0.2403];

% Optimize each dataset
xmin = zeros(4,n);
fval = zeros(1,n);
for i = 1:n
    fun = @(x) optshft(rho,I(:,:,:,i),mask,G,x);
    [xmin(:,i),fval(i)] = fminsearch(fun,x0(:,i),opt);
end

% Generate the simulated field
field = zeros(ny,nx,n_omega,n,'like',rho);
for i = 1:n
    px = -1i.*2.*pi.*xmin(1,i)./nx.*x;
    py = -1i.*2.*pi.*xmin(2,i)./ny.*y;
    pz = -1i.*2.*pi.*xmin(3,i)./n_omega.*z;
    tmp = rho.*exp(px).*exp(py).*exp(pz);
    field(:,:,:,i) = sqrt(xmin(4,i)).*fftshift(fftn(fftshift(tmp)));
end

% Convert simulated field to simulated intensity
field = abs(field).^2;

% Collect data from GPU
rho = gather(rho);
I = gather(I);
mask = gather(mask);
G = gather(G);
field = gather(field);

% Make RGB comparison plots
j0 = 141;
imin = 0;
rng = 6;
rgb = cat(5,field,I,field);
rgb = permute(rgb,[1 2 5 4 3]);
rgb = (log10(rgb) - imin)./rng;
rgb = rgb(:,:,:,[5 2 8 4 1 7 6 3 9],:);
figure;
set(gcf,'Position',[550 250 768 768],'Name',num2str(j0),'WindowScrollWheelFcn',{@interactive5D,rgb});
for i = 1:n
    subplot(3,3,i);
    image(rgb(:,:,:,i,j0));
    axis equal tight off;
    set(gca,'Position',[rem(i - 1,3)/3 (2 - floor((i - 1)/3))/3 1/3 1/3]);
end



% function fom = optshft(rho,I,mask,G,var)
% % Separate the shift vectors
% % bx = var(1:9);
% % by = var(10:18);
% % bz = var(19:27);
% bx = var(1);
% by = var(2);
% bz = var(3);
% a = var(4);
% 
% % Change the dimension of the shift vectors
% bx = permute(bx,[2 3 4 1]);
% by = permute(by,[2 3 4 1]);
% bz = permute(bz,[2 3 4 1]);
% 
% % Get the data size
% nx = size(I,2);
% ny = size(I,1);
% n_omega = size(I,3);
% 
% % Generate pixel coordinates
% x = (-nx/2):(nx/2 - 1);
% y = ((-ny/2):(ny/2 - 1)).';
% z = permute((-n_omega/2):(n_omega/2 - 1),[3 1 2]);
% 
% % Calculate the phase shifts
% px = -1i.*2.*pi.*bx./nx.*x;
% py = -1i.*2.*pi.*by./ny.*y;
% pz = -1i.*2.*pi.*bz./n_omega.*z;
% 
% % Send phase shifts to GPU
% px = gpuArray(px);
% py = gpuArray(py);
% pz = gpuArray(pz);
% 
% % Apply the phase factors
% tmp = rho.*exp(px).*exp(py).*exp(pz);
% 
% % Fourier transform to the detector space
% tmp = fftshift(fftshift(fftshift(tmp,1),2),3);
% tmp = fft(tmp);
% tmp = permute(tmp,[2 3 1 4]);
% tmp = fft(tmp);
% tmp = permute(tmp,[2 3 1 4]);
% tmp = fft(tmp);
% tmp = permute(tmp,[2 3 1 4]);
% tmp = fftshift(fftshift(fftshift(tmp,1),2),3);
% 
% % Calculate intensity
% tmp = abs(tmp).^2;
% 
% % Apply detector mask and pupil function
% tmp = tmp.*mask.*G;
% 
% % Calculate RMS difference between simulated and measured intensity
% fom = gather(sum(sum(sum(sum(abs(I - a.*tmp).^2)))));
% 
% end
% 
