% Initialization
clear;
close all;
clc;


%% Load the data
% Set the data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the measurement
tic;
load([dp 'Focus_Scan_1.mat']);
toc;

% Plot the data
Slicer(I,'displayRange',[0.7 1.3]);


%% Extract the lines
% Get the horizontal line top
t = I(678:858,1000:1048,28);

% Get the x- and y-axes
xt = ((1:size(t,2)).' - 21)*740e-9;
yt = ((1:size(t,1)).' - 0)*740e-9;

% Plot the top ROI
figure;
imagesc(xt,yt,t);
axis equal tight;

% Get the horizontal line bottom
b = I(1198:1370,1015:1052,29);

% Get the x- and y-axes
xb = ((1:size(b,2)).' - 19)*740e-9;
yb = ((1:size(b,1)).' - 0)*740e-9;

% Plot the bottom ROI
figure;
imagesc(xb,yb,b);
axis equal tight;

% Get the vertical line left
l = I(1009:1038,555:792,20);

% Get the x- and y-axes
xl = ((1:size(l,2)).' - 10)*740e-9;
yl = ((1:size(l,1)).' - 0)*740e-9;

% Plot the left ROI
figure;
imagesc(xl,yl,l);
axis equal tight;


%% Save the data
% Save the profiles
save('ROI_Resolution_1.mat','t','xt','yt','b','xb','yb','l','xl','yl','-v7.3');

