% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the electric field from the free space simulation
tic;
Ef = load_binary([p 'Ef.bin'],[ny nx n_omega]);
toc;

% Load the electric field from the lens-based simulation
tic;
El = load_binary([p 'El.bin'],[ny nx n_omega]);
toc;


%% Recover objects
% Set the phase curvatures
Rpxf = 4.0000;
Rpyf = 4.0000;
Rpxl = 6.0318;
Rpyl = 6.0318;

% Compensate the phase curvatures
tic;
xd = ((-nx/2):(nx/2 - 1)).*det_dx;
yd = ((-ny/2):(ny/2 - 1)).'.*det_dx;
Ef = Ef.*exp(-1i.*pi.*xd.^2./(lambda.*Rpxf)).*exp(-1i.*pi.*yd.^2./(lambda.*Rpyf));
El = El.*exp(-1i.*pi.*xd.^2./(lambda.*Rpxl)).*exp(-1i.*pi.*yd.^2./(lambda.*Rpyl));
toc;

% Generate the hybrid electric field
tic;
Eh = abs(El).*exp(1i.*angle(Ef));
toc;

% Get the free space object
tic;
F = ifftshift(ifftn(fftshift(Ef)));
toc;

% Get the lens-based object
tic;
L = ifftshift(ifftn(fftshift(El)));
% L = zeros(ny,nx,n_omega);
% for i = 1:n_omega
%     L(:,:,i) = frfft2gpu(El(:,:,i),-1-[4.0881e-4 3.9575e-4]);
% end
% L = ifftshift(ifft(fftshift(L,3),[],3),3);
toc;

% Get the hybrid object
tic;
H = ifftshift(ifftn(fftshift(Eh)));
toc;

% Shift the objects
tic;
for i = 1:n_omega
    k = -round((i - n_omega/2)*shft);
    F(:,:,i) = circshift(F(:,:,i),k,1);
    L(:,:,i) = circshift(L(:,:,i),k,1);
    H(:,:,i) = circshift(H(:,:,i),k,1);
end
toc;

% Get the absolute values
tic;
aF = abs(F);
aL = abs(L);
aH = abs(H);
toc;

% Plot the free space object
sf = Slicer(aF,'displayRange',[0 max(aF(:))]);

% Plot the lens-based object
sl = Slicer(aL,'displayRange',[0 max(aL(:))]);

% Plot the hybrid object
sh = Slicer(aH,'displayRange',[0 max(aH(:))]);

% Save the inversions
tic;
save_binary([p 'Inversion_Object_F.bin'],F);
save_binary([p 'Inversion_Object_L.bin'],L);
save_binary([p 'Inversion_Object_H.bin'],H);
toc;


