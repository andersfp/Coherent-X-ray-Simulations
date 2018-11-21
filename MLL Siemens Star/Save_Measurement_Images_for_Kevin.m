% Initialization
clear;
close all;
clc;


%% Load the corrected images
% Set the data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the measurement
load([dp 'Focus_Scan_1.mat']);


%% Save the images
% Convert the image data to 16 bit integer
I = uint16(round(32000*I));

% Set the base name and extension
f = 'Siemens_Star_';
ext = '.tiff';

% Save the images
n = size(I,3);
for i = 1:n
    imwrite(I(:,:,i),[f sprintf('%02.0f',i) ext]);
end

