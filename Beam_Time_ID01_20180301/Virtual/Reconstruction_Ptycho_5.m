% Initialization
clear;
close all;
clc;

% Use GPU?
gpu = 1;


%% Load data
% Load the diffraction patterns
load('Ptycho_Data.mat');

% Load an initial guess
dat = load('Reconstruction_Center_2018_04_17_16_58_1161.mat');
rho0 = dat.object_avg./dat.cyc;

% Load data corrections
cor = load('Ptycho_Corrections.mat');

% Set the correct data size
nx = size(I,2);
ny = size(I,1);
n_omega = size(I,3);
n = size(I,4);
rho0 = padarray(rho0,([ny nx n_omega] - size(rho0))/2,0);
field0 = fftshift(fftn(ifftshift(rho0)));


%% Reconstruction
% Generate pixel coordinates
x = (-nx/2):(nx/2 - 1);
y = ((-ny/2):(ny/2 - 1)).';
z = permute((-n_omega/2):(n_omega/2 - 1),[3 1 2]);

% Correct the intensity
I = I./permute(cor.a,[2 3 4 1]);
%I = I./permute(cor.af,[2 3 4 1]);

% Set initial guess of positions
bx = [0.0275   -0.2525   -0.3311  -28.6346  -28.8591  -28.7290   26.9370   26.8944   27.3283].';
by = [-0.9228  -38.1500   31.9645   -4.3790  -40.8335   30.4793   -4.9145  -41.7507   30.4394].';
bz = cor.z0;
a = [1.3025    1.4038    1.2582    1.3122    1.8048    1.3569    1.3033    2.3554    1.2643].';
%a = [1.3053    0.6485    1.2335    1.3532    0.5728    1.3166    1.3219    0.6535    1.2664].';

% Correct the slope of the initial reconstruction
pz = 1i.*2.*pi.*bz(1)./n_omega.*z;
rho0 = rho0.*exp(pz);

% Generate support
support = single(abs(rho0) > 0);

% Shrink the support
%nhood = [0 1 0;1 1 1;0 1 0];
%nhood = [0 0 1 0 0;0 0 1 0 0;1 1 1 1 1;0 0 1 0 0;0 0 1 0 0];
% nhood = [0 0 0 1 0 0 0;0 0 0 1 0 0 0;0 0 0 1 0 0 0;1 1 1 1 1 1 1;0 0 0 1 0 0 0;0 0 0 1 0 0 0;0 0 0 1 0 0 0];
% se = strel('arbitrary',nhood);
% support = imerode(support,se);
% nhood = [1 0 1;0 1 0;1 0 1];
% se = strel('arbitrary',nhood);
% support = imdilate(support,se);

% Load support
load('Ptycho_5_Support.mat');

% Generate amplitude data
data_input = sqrt(I);

% Generate an effective mask
mask = mask.*sqrt(G);

% Generate the initial object
rho = rho0;

% Send to GPU if used
if gpu == 1
    data_input = gpuArray(data_input);
    support = gpuArray(support);
    rho = gpuArray(rho);
    mask = gpuArray(mask);
    x = gpuArray(x);
    y = gpuArray(y);
    z = gpuArray(z);
end

% Set parameters
beta = 0.9;
epsilon = [];
chi_square_real = [];
chi_square_reci = [];
tol = 0.15;
sig0 = 2.0;
sigmin = 1.5;
sig = sig0;
ufo = 0.75;

cyci = 20;
cyco = 20;

bxp = NaN(n,cyco);
byp = bxp;
bzp = bxp;
ap = bxp;
bxp(:,1) = bx;
byp(:,1) = by;
bzp(:,1) = bz;
ap(:,1) = a;

co = parula(n);
set(groot,'DefaultAxesColorOrder',co);

n2 = 7;

% Run cycles
for iii = 1:cyco
    % Update the data positions
    if iii <= 4
        [bx,by,bz,a] = optshft_fit(bx,by,bz,a,rho,data_input,mask,1);
    end
    for ii = 1:cyci
        % Run ER on each dataset individually
        for i = 1:n2
            object = rho.*exp(-1i.*2.*pi.*bx(i)./nx.*x).*exp(-1i.*2.*pi.*by(i)./ny.*y).*exp(-1i.*2.*pi.*bz(i)./n_omega.*z);
            field = fftshift(fftn(fftshift(object)));
            %field_new = (data_input(:,:,:,i)./sqrt(a(i))).*exp(1i.*angle(field)).*mask + field.*(1 - mask);
            field_new = data_input(:,:,:,i).*exp(1i.*angle(field)).*mask + field.*(1 - mask);
            object_new = fftshift(ifftn(fftshift(field_new)));
            object = object_new.*support;
            rho = object.*exp(1i.*2.*pi.*bx(i)./nx.*x).*exp(1i.*2.*pi.*by(i)./ny.*y).*exp(1i.*2.*pi.*bz(i)./n_omega.*z);
            chi_square_real = error_metric_real(object_new,object,chi_square_real);
            chi_square_reci = error_metric_reci(field_new,field,chi_square_reci);
        end
        % Plots
        chi = reshape(chi_square_reci,n2,ii + (iii - 1)*cyci).';
        figure(10);
        semilogy(chi);
        figure(100);
        subplot(1,2,1);
        imagesc(x,y,support(:,:,n_omega/2+1));
        axis equal tight;
        subplot(1,2,2);
        imagesc(x,y,abs(rho(:,:,n_omega/2+1)));
        axis equal tight;
        drawnow;
    end
    % Shrinkwrap
    %[support,sig] = shrinkwrap(rho,x,y,z,sig,sigmin,tol,support);
    % Smoothness filter
    %rho = imgaussfilt3(abs(rho),ufo).*exp(1i.*angle(rho));
    % Plot alignment
    bxp(:,iii) = bx;
    byp(:,iii) = by;
    bzp(:,iii) = bz;
    ap(:,iii) = a;
    figure(20);
    subplot(1,4,1);
    plot(bxp.');
    subplot(1,4,2);
    plot(byp.');
    subplot(1,4,3);
    plot(bzp.');
    subplot(1,4,4);
    plot(ap.');
    drawnow;
end

% Apply the support at the end
%rho = rho.*support;

% Generate separate fields
fieldi = zeros(ny,nx,n_omega,n,'like',rho);
for i = 1:n
    px = -1i.*2.*pi.*bx(i)./nx.*x;
    py = -1i.*2.*pi.*by(i)./ny.*y;
    pz = -1i.*2.*pi.*bz(i)./n_omega.*z;
    tmp = rho.*exp(px).*exp(py).*exp(pz);
    fieldi(:,:,:,i) = fftshift(fftn(fftshift(tmp)));
end

% Genereate the field for the reconstructed object
field = fftshift(fftn(fftshift(rho)));

% Gather from GPU if used
if gpu == 1
    data_input = gather(data_input);
    support = gather(support);
    rho = gather(rho);
    mask = gather(mask);
    x = gather(x);
    y = gather(y);
    z = gather(z);
    chi_square_real = gather(chi_square_real);
    chi_square_reci = gather(chi_square_reci);
    field = gather(field);
    field_new = gather(field_new);
    object = gather(object);
    object_new = gather(object_new);
    fieldi = gather(fieldi);
end

% Generate an effective pupil function
Ge = zeros(ny,nx,n2);
for i = 1:n2
    Ge(:,:,i) = imshift2(sqrt(G),[bx(i) by(i)]);
end
Ge = abs(Ge);
Ge2 = sum(Ge,3);
Ge2(Ge2 > 1) = 1;
se2 = strel('disk',14,0);
Ge3 = imdilate(Ge2,se2);
field = field.*Ge3;
fieldi = fieldi.*Ge3;
rho = fftshift(ifftn(fftshift(field)));
rho = rho.*support;

% Make RGB comparison plots
j0 = 141;
imin = 0;
rng = 3;
rgb = cat(5,abs(fieldi),data_input,abs(fieldi));
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

% Make 3D plots
Slicer(log10(abs(field)),'displayRange',[0 3]);
Slicer(abs(rho));

% Save the reconstruction
save(['Reconstruction_Ptycho_5_' timestr(3) '_' num2str(round(sum(chi_square_reci(end-n:end))),'%.0f') '.mat'],...
    'rho','chi_square_real','chi_square_reci','cyci','cyco','ufo','tol','sig0','sigmin','sig','nx','ny','n_omega','n','n2','bx','by','bz','a','Ge');


%% Test
% s2 = single(abs(rho) > 0.012);
% se = strel('sphere',3);
% s3 = imerode(s2,se);
% s4 = bwareaopen(s3,100);
% s5 = imdilate(s4,se);
% s6 = zeros(size(s5));
% for i = 1:n_omega
%     s6(:,:,i) = imfill(s5(:,:,i),'holes');
% end
% Slicer(s6);
