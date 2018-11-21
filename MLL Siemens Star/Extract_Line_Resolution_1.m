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
t = squeeze(mean(I(678:858,1000:1048,:),1));

% Get the x- and y-axes
xt = ((1:size(t,1)).' - 21)*740e-9;
yt = t(:,28);

% Plot the top lines
figure;
plotInt(xt,t);
figure;
plot(xt,yt);

% Get the horizontal line bottom
b = squeeze(mean(I(1198:1370,1015:1052,:),1));

% Get the x- and y-axes
xb = ((1:size(b,1)).' - 19)*740e-9;
yb = b(:,29);

% Plot the bottom lines
figure;
plotInt(xb,b);
figure;
plot(xb,yb);

% Get the vertical line left
l = squeeze(mean(I(1009:1038,555:792,:),2));

% Get the x- and y-axes
xl = ((1:size(l,1)).' - 10)*740e-9;
yl = l(:,20);

% Plot the left lines
figure;
plotInt(xl,l);
figure;
plot(xl,yl);


%% Save the data
% Save the profiles
save('Line_Resolution_1.mat','t','xt','yt','b','xb','yb','l','xl','yl','-v7.3');

