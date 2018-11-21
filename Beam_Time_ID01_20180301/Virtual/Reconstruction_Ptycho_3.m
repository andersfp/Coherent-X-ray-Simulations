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

% Set initial guess of xy-position alignment of diffraction data
bx0 = [0 0 0 -1 -1 -1 1 1 1].'.*29; % 29
by0 = [0 -1 1 0 -1 1 0 -1 1].'.*38 + [0 0 0 -1 -1 -1 1 1 1].'.*2; % 38, 2
bx = single(permute(bx0,[2 3 4 1]));
by = single(permute(by0,[2 3 4 1]));

% Calculate guess for z-position alignment of diffraction data
A = sum(sum(I));
bz = -sum(A.*z)./sum(A);
bz = bz - bz(1);
bz0 = squeeze(bz);

% Generate support
support = single(abs(rho0) > 0);

% Generate amplitude data
data_input = sqrt(I);

% Generate an effective mask
mask = mask.*G;

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

% Generate the individual objects
rhoi = rho.*exp(-1i.*2.*pi.*bx./nx.*x).*exp(-1i.*2.*pi.*by./ny.*y).*exp(-1i.*2.*pi.*bz./n_omega.*z);

% Generate the initial field
field = fftshift(fftn(ifftshift(rho)));

% Set parameters
beta = 0.9;
epsilon = [];
chi_square_real = [];
chi_square_reci = [];
tol = 0.12;
sig0 = 2.0;
sigmin = 1.5;
sig = sig0;

ER = 2;
iterER1 = 1;
cyci = 1;
cyco = 1;

bxp = NaN(n,cyco+1);
byp = bxp;
bzp = bxp;
bxp(:,1) = bx0;
byp(:,1) = by0;
bzp(:,1) = bz0;

co = hsv(n);
set(groot,'DefaultAxesColorOrder',co);

% Run cycles
for iii = 1:cyco
    for ii = 1:cyci
        % Update the individual estimates
        rhoi = rho.*exp(-1i.*2.*pi.*bx./nx.*x).*exp(-1i.*2.*pi.*by./ny.*y).*exp(-1i.*2.*pi.*bz./n_omega.*z);
        % Run ER on each dataset individually
        for i = 1:n
            % ER
            [field,rhoi(:,:,:,i),chi_square_real,chi_square_reci] = cycle(ER,iterER1,beta,epsilon,data_input(:,:,:,i),support,field,rhoi(:,:,:,i),mask,chi_square_real,chi_square_reci);
        end
        % Make a new overall estimate
        rho = mean(rhoi.*exp(1i.*2.*pi.*bx./nx.*x).*exp(1i.*2.*pi.*by./ny.*y).*exp(1i.*2.*pi.*bz./n_omega.*z),4);
        % Plots
        a = reshape(chi_square_reci,iterER1,n,ii + (iii - 1)*cyci);
        a = permute(a,[1 3 2]);
        a = reshape(a,iterER1*(ii + (iii - 1)*cyci),n);
        figure(10);
        semilogy(a);
        figure(100);
        subplot(1,2,1);
        imagesc(x,y,support(:,:,n_omega/2+1));
        axis equal tight;
        subplot(1,2,2);
        imagesc(x,y,abs(rho(:,:,n_omega/2+1)));
        axis equal tight;
        drawnow;
    end
    % Update the phase slopes
    [bx,by,bz] = phaseSlopeFit(rhoi);
    % Shrinkwrap
    %[support,sig] = shrinkwrap(rho,x,y,z,sig,sigmin,tol,support);
    % Plot alignment
    bxp(:,iii+1) = squeeze(bx);
    byp(:,iii+1) = squeeze(by);
    bzp(:,iii+1) = squeeze(bz);
    figure(20);
    subplot(1,3,1);
    plot(bxp.');
    subplot(1,3,2);
    plot(byp.');
    subplot(1,3,3);
    plot(bzp.');
    drawnow;
end

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
    rhoi = gather(rhoi);
end

% Make RGB comparison plots
j0 = 141;
imin = 0;
rng = 3;
rgb = cat(5,repmat(abs(field),1,1,1,n),data_input,zeros(ny,nx,n_omega,n));
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
Slicer(log10(abs(field)),'displayRange',[-3 3]);
Slicer(abs(rho));


return
% Save the reconstruction
save(['Reconstruction_Ptycho_3_' timestr(3) '_' num2str(round(chi_square_reci(end)),'%.0f') '.mat'],...
    'rho','chi_square_real','chi_square_reci','ER','iterER1','cyci','cyco','beta','tol','sig0','sigmin','sig','nx','ny','n_omega');


