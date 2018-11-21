% Initialization
clear;
close all;
clc;


%% Define the file
% Set the data path
p1 = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set1_2D_17.3keV\PCO_dark\';
p2 = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set3_2D_34.5keV\pco_bg_3s\';
p3 = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_17.3keV_darkfield\pco_dark_30s\';

% Set the base name
fp = '*.edf';

% Make the file lists
[~,f1] = fileNames(p1,fp);
[~,f2] = fileNames(p2,fp);
[~,f3] = fileNames(p3,fp);
f = cat(1,f1,f2,f3);


%% Load the data
% Get number of images
n = length(f);

% Pre-allocate info structure
info(n) = edf_info(f{n});

% Set the image dimensions
mx = info(n).dim_1;
my = info(n).dim_2;

% Read the images
img = zeros(my,mx,n);
fprintf([num2str(n) ': ']);
for i = 1:n
    img(:,:,i) = edf_read(f{i});
    info(i) = edf_info(f{i});
    fprintf('.');
    if rem(i,100) == 0
        fprintf('\n');
    end
end
fprintf('\n');

% Plot the background images
Slicer(img,'displayRange',[80 120]);


%% Create a mean background
% Average background without rectangle contribution
bg = median(img(:,:,[1:20 26:30]),3);

% Plot the background
figure;
imagesc(bg,[80 120]);
axis equal tight;
colormap gray;

% Save the background
save('Background.mat','bg','-v7.3');


