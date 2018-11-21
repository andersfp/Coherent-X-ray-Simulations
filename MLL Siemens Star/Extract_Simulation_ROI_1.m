% Initialization
clear;
close all;
clc;


%% Load the data
% Set the data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the simulation
tic;
B = load([dp 'Star_Simulation_Series_1.mat']);
toc;
S = B.I;
d = B.d.';
x0 = B.x;
y0 = B.y;

% Clear unused data
clear B;

% Plot the simulation
Slicer(S,'displayRange',[16000 52000]);


%% Get the line profiles
% Set the positions
ix1 = 1000;
ix2 = 1049;
iy1 = 681;
iy2 = 780;

% Get the ROI
St = S(iy1:iy2,ix1:ix2,:);
xt = x0(ix1:ix2);
yt = y0(iy1:iy2);

% Set positions
ix1 = 581;
ix2 = 680;
iy1 = 1000;
iy2 = 1049;

% Get the ROI
Sl = S(iy1:iy2,ix1:ix2,:);
xl = x0(ix1:ix2);
yl = y0(iy1:iy2);

% Save the ROIs
save('ROI_Simulation_1.mat','d','St','xt','yt','Sl','xl','yl','-v7.3');



