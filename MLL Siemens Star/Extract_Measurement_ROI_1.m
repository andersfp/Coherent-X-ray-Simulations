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

% Clear unused data
clear A;


%% Get the line profiles
% Plot the simulation
Slicer(D,'displayRange',[0.7 1.3]);

% Set the axes
x0 = (-1024:1023)*0.74e-6;
y0 = x0.';

% Fitting model
fun = @(a,b,d,theta,sig,x,y) a.*erfc((cosd(theta).*x - sind(theta).*y - b)./(sqrt(2).*sig)) + d;

% Fit the ROI
Z = D(675:858,998:1043,28);
[X,Y] = meshgrid(x0(998:1043),y0(675:858));
ft = fit([X(:) Y(:)],Z(:),fun,'StartPoint',[(mean(Z(:,1))-mean(Z(:,end)))/2 0 mean(Z(:,end)) 2.3 2.6e-6]);

% Plot the fit
figure;
plot(ft,[X(:) Y(:)],Z(:));

% Shift the measurement
shft = sind(ft.theta);
m = size(D,1);
n = size(D,3);
D0 = zeros(m,m,n);
for i = 1:m
    D0(i,:,:) = circshift(D(i,:,:),round((m/2 - i).*shft),2);
end
Slicer(D0,'displayRange',[0.7 1.3]);

% Set the positions
ix1 = 1000;
ix2 = 1049;
iy1 = 681;
iy2 = 780;

% Get the ROI
Dt = D0(iy1:iy2,ix1:ix2,:);
xt = x0(ix1:ix2);
yt = y0(iy1:iy2);

% Fit the ROI
Z = D(1200:1370,1013:1055,29);
[X,Y] = meshgrid(x0(1013:1055),y0(1200:1370));
fb = fit([X(:) Y(:)],Z(:),fun,'StartPoint',[(mean(Z(:,end))-mean(Z(:,1)))/2 0 mean(Z(:,1)) 181.1 1.8e-6]);

% Plot the fit
figure;
plot(fb,[X(:) Y(:)],Z(:));

% Shift the measurement
shft = sind(fb.theta - 180);
m = size(D,1);
n = size(D,3);
D0 = zeros(m,m,n);
for i = 1:m
    D0(i,:,:) = circshift(D(i,:,:),round((m/2 - i).*shft),2);
end
Slicer(D0,'displayRange',[0.7 1.3]);

% Set the positions
ix1 = 1000;
ix2 = 1049;
iy1 = 1261;
iy2 = 1360;

% Get the ROI
Db = D0(iy1:iy2,ix1:ix2,:);
xb = x0(ix1:ix2);
yb = y0(iy1:iy2);

% Fit the ROI
Z = D(1006:1041,562:792,20);
[X,Y] = meshgrid(x0(562:792),y0(1006:1041));
fl = fit([X(:) Y(:)],Z(:),fun,'StartPoint',[(mean(Z(end,:))-mean(Z(1,:)))/2 0 mean(Z(1,:)) 88.9 5.2e-6]);

% Plot the fit
figure;
plot(fl,[X(:) Y(:)],Z(:));

% Shift the measurement
shft = -sind(fl.theta - 90);
m = size(D,1);
n = size(D,3);
D0 = zeros(m,m,n);
for i = 1:m
    D0(:,i,:) = circshift(D(:,i,:),round((m/2 - i).*shft),1);
end
Slicer(D0,'displayRange',[0.7 1.3]);

% Set the positions
ix1 = 581;
ix2 = 680;
iy1 = 1000;
iy2 = 1049;

% Get the ROI
Dl = D0(iy1:iy2,ix1:ix2,:);
xl = x0(ix1:ix2);
yl = y0(iy1:iy2);

% Save the ROIs
save('ROI_Measurement_1.mat','Dt','xt','yt','Dl','xl','yl','Db','xb','yb','-v7.3');



