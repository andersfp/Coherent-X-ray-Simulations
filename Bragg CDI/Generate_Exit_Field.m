% Initialization
clear;
close all;
clc;


%% Load data
% Load the data path
load('Exp_Param.mat');

% Load the object
%load([p 'Object.mat'],'O123');
O123 = load_binary([p 'O123.bin'],[ny nx n_omega]);


%% Generate the exit field
% 3D Fourier transform the object
F = ifftshift(fftn(fftshift(O123)));

% 2D Inverse Fourier transform the slices
n = size(F,3);
E0 = zeros(size(F));
for i = 1:n
    E0(:,:,i) = ifftshift(ifft2(fftshift(F(:,:,i))));
end


%% Save the result
% Save the exit field
%save([p 'Exit_Field.mat'],'E0','-v7.3');
save_binary([p 'E0.bin'],E0);


