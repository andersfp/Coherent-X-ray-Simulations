% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental details
load('Exp_Param.mat');

% Load the diffraction patterns
%load([p 'Diffraction_Intensity.mat'],'If');
If1 = load_binary([p 'If1.bin'],[ny nx n_omega]);


%% Support reconstruction
% Downsample the measured intensity
i1 = 404;
i2 = 387;
i3 = 122;
I = If1(i1+1:end-i1,i2+1:end-i2,i3+1:end-i3);

% Use GPU
gpu = 1;

% Set up reconstruction parameters
beta = 0.9;
epsilon = 0.0001;
sigmin = 2.0;
tol = 0.25;
sig = 3.0;

% Coordinates
%[X,Y,Z] = meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1,-n_omega/2:n_omega/2-1);
X = ((-size(I,2)/2):(size(I,2)/2-1));
Y = ((-size(I,1)/2):(size(I,1)/2-1)).';
Z = permute(((-size(I,3)/2):(size(I,3)/2-1)),[3 1 2]);

% Perform the reconstruction
data_input = sqrt(I);
support = ifftshift(fftn(fftshift(data_input)));
support = double(abs(support) > 0.002*max(abs(support(:))));
field = data_input.*exp(1i*2*pi*rand(size(I,1),size(I,2),size(I,3)));
object = ifftshift(ifftn(fftshift(field)));
mask = 1;
chi_square_real = [];
chi_square_reci = [];

% Send to GPU if used
if gpu == 1
    data_input = gpuArray(data_input);
    support = gpuArray(support);
    field = gpuArray(field);
    object = gpuArray(object);
    mask = gpuArray(mask);
end

% Perform the support reconstruction
tic;
% for i = 1:30
%     [field,object,chi_square_real,chi_square_reci] = cycle(8,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
%     [field,object,chi_square_real,chi_square_reci] = cycle(18,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
%     figure(10);
%     semilogy(chi_square_reci);
%     drawnow;
%     [support,sig] = shrinkwrap(object,X,Y,Z,sig,sigmin,tol);
%     figure(100 + i);
%     imagesc(x,y,support(:,:,size(I,3)/2+1));
%     axis equal tight;
%     drawnow;
% end
i = 0;
while sig > sigmin
    i = i + 1;
    [field,object,chi_square_real,chi_square_reci] = cycle(8,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    [field,object,chi_square_real,chi_square_reci] = cycle(18,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    figure(10);
    semilogy(chi_square_reci);
    drawnow;
    [support,sig] = shrinkwrap(object,X,Y,Z,sig,sigmin,tol);
    figure(100 + i);
    imagesc(x,y,support(:,:,size(I,3)/2+1));
    axis equal tight;
    drawnow;
end
% sig2 = linspace(sig,sigmin,round((sig - sigmin)/(0.01*sig))).';
% for i = 1:length(sig2)
%     [field,object,chi_square_real,chi_square_reci] = cycle(8,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
%     [field,object,chi_square_real,chi_square_reci] = cycle(18,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
%     figure(10);
%     semilogy(chi_square_reci);
%     drawnow;
%     [support,sig] = shrinkwrap(object,X,Y,Z,sig2(i),sigmin,tol);
%     figure(100 + i);
%     imagesc(x,y,support(:,:,size(I,3)/2+1));
%     axis equal tight;
%     drawnow;
% end
% i = 0;
% while sig > sigmin
%     i = i + 1;
%     [field,object,chi_square_real,chi_square_reci] = cycle(8,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
%     figure(10);
%     semilogy(chi_square_reci);
%     drawnow;
%     [support,sig] = shrinkwrap(object,X,Y,Z,sig,sigmin,tol);
%     figure(100 + 2*i - 1);
%     imagesc(x,y,support(:,:,size(I,3)/2+1));
%     axis equal tight;
%     drawnow;
%     [field,object,chi_square_real,chi_square_reci] = cycle(18,10,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
%     figure(10);
%     semilogy(chi_square_reci);
%     drawnow;
%     [support,sig] = shrinkwrap(object,X,Y,Z,sig,sigmin,tol);
%     figure(100 + 2*i);
%     imagesc(x,y,support(:,:,size(I,3)/2+1));
%     axis equal tight;
%     drawnow;
% end
toc;

% Collect data from GPU if used
if gpu == 1
    data_input = gather(data_input);
    support = gather(support);
    field = gather(field);
    object = gather(object);
    mask = gather(mask);
end

% Plot the support
figure(20);
imagesc(1e6*x,1e6*y,support(:,:,n_omega/2+1));
axis equal tight;
title('Support');
xlabel('x [\mum]');
ylabel('y [\mum]');

% Save the resulting support
save([p 'Lowres_Free_Space_Shrinkwrap_1.mat'],'support','field','i1','i2','i3','-v7.3');


