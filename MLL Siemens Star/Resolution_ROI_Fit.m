% Initialization
clear;
close all;
clc;


%% Load data
% Load the line profiles
load('ROI_Resolution_1.mat');


%% Fit the data
% Make a fitting model
fun = @(a,b,d,theta,sig,x,y) a.*erfc((cosd(theta).*x - sind(theta).*y - b)./(sqrt(2).*sig)) + d;

% Fit the top line
[Xt,Yt] = meshgrid(xt,yt);
ft = fit([Xt(:) Yt(:)],t(:),fun,'StartPoint',[(mean(t(:,1))-mean(t(:,end)))/2 0 min(t(:)) 2.4 2.5e-6]);

% Plot the top line fit
figure;
plot(ft,[Xt(:) Yt(:)],t(:));

% Fit the bottom line
[Xb,Yb] = meshgrid(xb,yb);
fb = fit([Xb(:) Yb(:)],b(:),fun,'StartPoint',[(mean(b(:,1))-mean(b(:,end)))/2 0 min(b(:)) 1 1.8e-6]);

% Plot the bottom line fit
figure;
plot(fb,[Xb(:) Yb(:)],b(:));

% Fit the left line
[Xl,Yl] = meshgrid(xl,yl);
fl = fit([Xl(:) Yl(:)],l(:),fun,'StartPoint',[(mean(l(1,:))-mean(l(end,:)))/2 0 min(l(:)) -91 0.5e-6]);

% Plot the left line fit
figure;
plot(fl,[Xl(:) Yl(:)],l(:));


%% Get the fitting results
% Set the confidence level
lvl = 0.682689492;
is = 5;

% Top line width
st = ft.sig;
ste = confint(ft,lvl);
ste = diff(ste(:,is));

% Bottom line width
sb = fb.sig;
sbe = confint(fb,lvl);
sbe = diff(sbe(:,is));

% Left line width
sl = fl.sig;
sle = confint(fl,lvl);
sle = diff(sle(:,is));

% Get the magnifications
Mv = (1404 - 658)*740e-9/(2*6e-6);
Mh = (1022 - 529)*740e-9/(1*6e-6);

% Get the real space resolutions
rt = st./Mh;
rte = ste./Mh;
rb = sb./Mh;
rbe = sbe./Mh;
rl = sl./Mv;
rle = sle./Mv;

disp([rt rte]*1e9);
disp([rb rbe]*1e9);
disp([rl rle]*1e9);


