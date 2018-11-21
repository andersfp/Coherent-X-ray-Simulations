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

% Reduce data size
% G = G(65:192,65:192);
% I = I(65:192,65:192,:,:);
% mask = mask(65:192,65:192);
% rho = bin3(rho,[2 2 1]);
% rho(end+1,end+1,:) = 0;
% rho = 4*rho;

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

% Load the intensity corrections
cor = load('Ptycho_Corrections.mat');


%% Fit the displacements
% Generate pixel coordinates
x = (-nx/2):(nx/2 - 1);
y = ((-ny/2):(ny/2 - 1)).';
z = permute((-n_omega/2):(n_omega/2 - 1),[3 1 2]);

% Set initial guess of xy-position alignment of diffraction data
bx0 = [0 0 0 -1 -1 -1 1 1 1].'.*28; % 29
by0 = [0 -1 1 0 -1 1 0 -1 1].'.*38 + [0 0 0 -1 -1 -1 1 1 1].'.*2; % 38, 2

% Correct the intensity
%I = I./permute(cor.a,[2 3 4 1]);
I = I./permute(cor.af,[2 3 4 1]);

% Correct the slope of the initial reconstruction
bz = cor.z0;
pz = 1i.*2.*pi.*bz(1)./n_omega.*z;
rho = rho.*exp(pz);

% Set the initial guess of the intensity matching
a0 = ones(n,1);

% Starting guess after running the optimization
bx0 = [0.0275   -0.2525   -0.3311  -28.6346  -28.8591  -28.7290   26.9370   26.8944   27.3283].';
by0 = [-0.9228  -38.1500   31.9645   -4.3790  -40.8335   30.4793   -4.9145  -41.7507   30.4394].';
%a0 = [1.3025    1.4038    1.2582    1.3122    1.8048    1.3569    1.3033    2.3554    1.2643].';
a0 = [1.3053    0.6485    1.2335    1.3532    0.5728    1.3166    1.3219    0.6535    1.2664].';

% Send data to GPU
rho = gpuArray(rho);
I = gpuArray(I);
mask = gpuArray(mask);
G = gpuArray(G);

% Optimization options
opt = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1e8,'TolX',1e-3);

% Optimize each dataset
xmin = zeros(3,n);
fval = zeros(1,n);
for i = 1:n
    fun = @(x) optshft(rho,I(:,:,:,i),mask,G,x(1),x(2),bz(i),x(3));
    [xmin(:,i),fval(i)] = fminsearch(fun,[bx0(i) by0(i) a0(i)],opt);
end

% Assign the results
bx = xmin(1,:).';
by = xmin(2,:).';
a = xmin(3,:).';

% Generate the simulated field
field = zeros(ny,nx,n_omega,n,'like',rho);
for i = 1:n
    px = -1i.*2.*pi.*bx(i)./nx.*x;
    py = -1i.*2.*pi.*by(i)./ny.*y;
    pz = -1i.*2.*pi.*bz(i)./n_omega.*z;
    tmp = rho.*exp(px).*exp(py).*exp(pz);
    field(:,:,:,i) = sqrt(a(i)).*fftshift(fftn(fftshift(tmp)));
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

