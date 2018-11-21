% Initialization
clear;
close all;
clc;


%% Input parameters
% Setup type
s = 'Lens';

% Save the line profiles?
sav = 0;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Load 1
tic;
O1 = load_binary([p 'Reconstruction_' s '_Shrinkwrap_1_object_xyz.bin'],[ny nx n_omega]);
toc;

% Load 1 averaged
tic;
OA1 = load_binary([p 'Reconstruction_' s '_Shrinkwrap_1_object_avg_xyz.bin'],[ny nx n_omega]);
toc;

% Load 2
tic;
O2 = load_binary([p 'Reconstruction_' s '_Shrinkwrap_2_object_xyz.bin'],[ny nx n_omega]);
toc;

% Load 2 averaged
tic;
OA2 = load_binary([p 'Reconstruction_' s '_Shrinkwrap_2_object_avg_xyz.bin'],[ny nx n_omega]);
toc;


%% Process the data
% Set the indices
% Free space 1: 442 498 267
% Free space 2: 525 524 224
% Lens 1: 552 498 258
% Lens 2: 456 483 269
if strcmp(s,'Free_Space')
    i11 = 442;
    i12 = 498;
    i13 = 267;
    i21 = 525;
    i22 = 524;
    i23 = 224;
elseif strcmp(s,'Lens')
    i11 = 552;
    i12 = 498;
    i13 = 258;
    i21 = 456;
    i22 = 483;
    i23 = 269;
else
    error('Wrong file name.');
end

% Get the lines
lx1 = abs(O1(i11,:,i13)).';
ly1 = abs(O1(:,i12,i13));
lz1 = squeeze(abs(O1(i11,i12,:)));
lxa1 = abs(OA1(i11,:,i13)).';
lya1 = abs(OA1(:,i12,i13));
lza1 = squeeze(abs(OA1(i11,i12,:)));
lx2 = abs(O2(i21,:,i23)).';
ly2 = abs(O2(:,i22,i23));
lz2 = squeeze(abs(O2(i21,i22,:)));
lxa2 = abs(OA2(i21,:,i23)).';
lya2 = abs(OA2(:,i22,i23));
lza2 = squeeze(abs(OA2(i21,i22,:)));

% Plot the lines
figure;
plot(x,[lx1 lxa1 lx2 lxa2]);
figure;
plot(y,[ly1 lya1 ly2 lya2]);
figure;
plot(z,[lz1 lza1 lz2 lza2]);

% Save the line profiles
if sav == 1
    save([s '_Line_Profiles.mat'],'lx1','lx2','lxa1','lxa2','ly1','ly2','lya1','lya2','lz1','lz2','lza1','lza2','x','y','z');
end


