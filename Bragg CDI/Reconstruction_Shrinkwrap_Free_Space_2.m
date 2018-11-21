% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental details
load('Exp_Param.mat');

% Load the diffraction patterns
If2 = load_binary([p 'If2.bin'],[ny nx n_omega]);

% Load the low resolution initial results
load([p 'Lowres_Free_Space_Shrinkwrap_2.mat']);


%% Reconstruction
% Use GPU
gpu = 0;

% Set up reconstruction parameters
iterER = 10;
iterHIO = 10;
cyc = 50;
beta = 0.9;
epsilon = 0.0001;
sigmin = 0.0;
tol = 0.15; % 0.25
sig = 3.0; % 3.0

% Coordinates
%[X,Y,Z] = meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1,-n_omega/2:n_omega/2-1);
X = ((-nx/2):(nx/2-1));
Y = ((-ny/2):(ny/2-1)).';
Z = permute(((-n_omega/2):(n_omega/2-1)),[3 1 2]);

% Perform the reconstruction
data_input = sqrt(If2);
support = ifftshift(fftn(fftshift(support)));
support = padarray(support,[i1 i2 i3],0,'both');
support = abs(fftshift(ifftn(ifftshift(support))));
support = double(support > max(support(:))/3);
field = padarray(field,[i1 i2 i3],0,'both');
object = fftshift(ifftn(fftshift(field)));
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

% Perform the reconstruction
tic;
for i = 1:cyc
    [field,object,chi_square_real,chi_square_reci] = cycle(8,iterHIO,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    figure(10);
    semilogy(chi_square_reci);
    drawnow;
    [field,object,chi_square_real,chi_square_reci] = cycle(18,iterER,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    figure(10);
    semilogy(chi_square_reci);
    drawnow;
    [support,sig] = shrinkwrap(object,X,Y,Z,sig,sigmin,tol);
    figure(100 + i);
    imagesc(x,y,support(:,:,n_omega/2+1));
    axis equal tight;
    drawnow;
end
toc;

% % Pause and save the reconstruction halfway to run over 2 nights
% save([p 'temp.mat'],'-v7.3');
% return
% clear
% load('Exp_Param.mat');
% load([p 'temp.mat']);

% Obtain averaging of reconstructions
object_avg = object;
tic;
for i = 1:cyc
    [field,object,chi_square_real,chi_square_reci] = cycle(8,iterHIO,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    figure(10);
    semilogy(chi_square_reci);
    drawnow;
    [field,object,chi_square_real,chi_square_reci] = cycle(18,iterER,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
    figure(10);
    semilogy(chi_square_reci);
    drawnow;
    object_avg = object_avg + object;
end
toc;

% Collect data from GPU if used
if gpu == 1
    data_input = gather(data_input);
    support = gather(support);
    field = gather(field);
    object = gather(object);
    mask = gather(mask);
end

% Save the error metric figure for quick reconstruction estimation
figure(10);
saveas(gcf,'Error_Metric_Shrinkwrap_Free_Space_2.png');

% Plot the support
figure(20);
imagesc(1e6*x,1e6*y,support(:,:,n_omega/2+1));
axis equal tight;
title('Support');
xlabel('x [\mum]');
ylabel('y [\mum]');

% % Plot the mask
% figure(30);
% imagesc(1e-9*qx,1e-9*qy,mask(:,:,n_omega/2+1));
% axis equal tight;
% title('Mask');
% xlabel('q_x [nm^{-1}]');
% ylabel('q_y [nm^{-1}]');

% Plot the reciprocal space reconstruction
figure(40);
subplot(1,2,1);
imagesc(1e-9*qx,1e-9*qy,log10(abs(field(:,:,n_omega/2+1))));
axis equal tight;
title('Reciprocal space amplitude');
xlabel('q_x [nm^{-1}]');
ylabel('q_y [nm^{-1}]');
subplot(1,2,2);
imagesc(1e-9*qx,1e-9*qy,angle(field(:,:,n_omega/2+1)));
axis equal tight;
title('Reciprocal space phase');
xlabel('q_x [nm^{-1}]');
ylabel('q_y [nm^{-1}]');

% Plot the reconstruction
figure(50);
subplot(1,2,1);
imagesc(1e6*x,1e6*y,abs(object(:,:,n_omega/2+1)));
axis equal tight;
title('Object amplitude');
xlabel('x [\mum]');
ylabel('y [\mum]');
subplot(1,2,2);
imagesc(1e6*x,1e6*y,angle(object(:,:,n_omega/2+1)));
axis equal tight;
title('Object phase');
xlabel('x [\mum]');
ylabel('y [\mum]');

% % Make 3D plot of the raw reconstructed object
% figure;
% pp = patch(isosurface(abs(object),0.1*max(abs(object(:)))));
% isonormals(abs(object),pp);
% pp.FaceColor = 'red';
% pp.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight(45,45);
% camlight(-135,-45);
% lighting gouraud;
% xlabel('y');
% ylabel('z');
% zlabel('x');

% Shift the object
object_xyz = object;
object_avg_xyz = object_avg;
for i = 1:n_omega
    object_xyz(:,:,i) = circshift(object(:,:,i),-round((i - n_omega/2)*shft),1);
    object_avg_xyz(:,:,i) = circshift(object_avg(:,:,i),-round((i - n_omega/2)*shft),1);
end

% % Make 3D plot of the orthogonal space reconstructed object
% figure;
% pp = patch(isosurface(abs(object_xyz),0.1*max(abs(object_xyz(:)))));
% isonormals(abs(object_xyz),pp);
% pp.FaceColor = 'red';
% pp.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight(45,45);
% camlight(-135,-45);
% lighting gouraud;
% xlabel('y');
% ylabel('z');
% zlabel('x');

% Save the resulting objects
save_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_xyz.bin'],object_xyz);
save_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object.bin'],object);
save_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_avg_xyz.bin'],object_avg_xyz);
save_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_object_avg.bin'],object_avg);
save_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_field.bin'],field);
save_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_support.bin'],support);

% Save the error metrics
save([p 'Reconstruction_Free_Space_Shrinkwrap_2_Error_Metrics.mat'],'chi_square_real','chi_square_reci');

