% Initialization
clear;
close all;
clc;


%% Load the data
% Set the data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the measurement
tic;
A = load([dp 'Focus_Scan_1.mat']);
toc;
D = A.I;
info = A.info;

% Load the simulation
tic;
B = load([dp 'Star_Simulation_Series_1.mat']);
toc;
S = B.I;
d = B.d.';
x0 = B.x;
y0 = B.y;

% Clear unused data
clear A B;


%% Plot the slice data
% Plot the measurement
Slicer(D,'displayRange',[0.7 1.3]);

% Plot the simulation
Slicer(S,'displayRange',[16000 52000]);


%% Get the line profiles
% Set the positions
w = 2;
ix = 600;
iy = 700;
ix1 = 865;
ix2 = 1100;
iy1 = 865;
iy2 = 1125;

% Get the data line profiles
dv = squeeze(mean(D(iy1:iy2,ix-w:ix+w,:),2));
dh = squeeze(mean(D(iy-w:iy+w,ix1:ix2,:),1));

% Get the data line profiles
sv = squeeze(mean(S(iy1:iy2,ix-w:ix+w,:),2));
sh = squeeze(mean(S(iy-w:iy+w,ix1:ix2,:),1));

% Get the axes
x = x0(ix1:ix2);
y = y0(iy1:iy2);

% Save the line profiles
%save('Line_Profiles_1.mat','d','dh','dv','sh','sv','x','y','-v7.3');


%% Plot the position of the lines
% Get the center image of the simulation as RGB
rgbs = repmat(S(:,:,41),1,1,3);
rgbs = (rgbs - min(rgbs(:)))./(max(rgbs(:)) - min(rgbs(:)));

% Highlight the line profiles
pw = w;
rgbs(iy1:iy2,ix-pw:ix+pw,2:3) = 0;
rgbs(iy-pw:iy+pw,ix1:ix2,1:2) = 0;

% Plot the position of the lines on the simulation
figure;
image(rgbs);
axis equal tight;
set(gca,'XLim',[400 1400],'YLim',[500 1500]);

% Get the center image of the data as RGB
rgbd = repmat(D(:,:,24),1,1,3);
rgbd = (rgbd - 0.7)./(1.3 - 0.7);

% Highlight the line profiles
rgbd(iy1:iy2,ix-pw:ix+pw,2:3) = 0;
rgbd(iy-pw:iy+pw,ix1:ix2,1:2) = 0;

% Plot the position of the lines on the simulation
figure;
image(rgbd);
axis equal tight;
set(gca,'XLim',[400 1400],'YLim',[500 1500]);

% Plot both images in one figure
figure;
subplot(1,2,1);
image(rgbs);
axis equal tight;
set(gca,'XLim',[400 1400],'YLim',[500 1500]);
subplot(1,2,2);
image(rgbd);
axis equal tight;
set(gca,'XLim',[400 1400],'YLim',[500 1500]);


