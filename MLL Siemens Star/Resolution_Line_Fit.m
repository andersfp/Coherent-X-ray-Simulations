% Initialization
clear;
close all;
clc;


%% Load data
% Load the line profiles
load('Line_Resolution_1.mat');


%% Fit the data
% Make a fitting model
fun = @(a,b,d,sig,x) a*erfc((x - b)./(sqrt(2)*sig)) + d;

% Fit the top line
ft = fit(xt,yt,fun,'StartPoint',[(yt(1)-yt(end))/2 0 min(yt) 2e-6]);

% Plot the top line fit
figure;
plot(ft,xt,yt);

% Fit the bottom line
fb = fit(xb,yb,fun,'StartPoint',[(yb(1)-yb(end))/2 0 min(yb) 2e-6]);

% Plot the bottom line fit
figure;
plot(fb,xb,yb);

% Fit the left line
fl = fit(xl,yl,fun,'StartPoint',[(yl(1)-yl(end))/2 0 min(yl) 2e-6]);

% Plot the left line fit
figure;
plot(fl,xl,yl);


%% Get the fitting results
% Set the confidence level
lvl = 0.682689492;

% Top line width
st = ft.sig;
ste = confint(ft,lvl);
ste = diff(ste(:,4));

% Bottom line width
sb = fb.sig;
sbe = confint(fb,lvl);
sbe = diff(sbe(:,4));

% Left line width
sl = fl.sig;
sle = confint(fl,lvl);
sle = diff(sle(:,4));

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


