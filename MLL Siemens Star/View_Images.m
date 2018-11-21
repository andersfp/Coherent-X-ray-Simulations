% Initialization
clear;
close all;
clc;


%% Data location
% Set the directory
%p = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_2D_17.3keV\pco_focus_30s\';
%p = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_2D_17.3keV\smx_star_scan_1\';
p = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set3_2D_34.5keV\pco_focus_3s\';

% Set the base name
fp = '*.edf';

% Make the file lists
[~,f] = fileNames(p,fp);


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


%% Display the data
% Cut out an ROI
ix1 = 1;
ix2 = mx;
iy1 = 1;
iy2 = my;

% Cut the image
img = img(iy1:iy2,ix1:ix2,:);

% Plot the data
Slicer(img);


