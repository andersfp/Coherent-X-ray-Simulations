% Initialization
clear;
close all;
clc;


%% Load data
% Load the line profiles
load('Line_Profiles.mat');


%% Process the lines
% Scale the real space axes
sc = 1e6;
rl = rl.*sc;
rd = rd.*sc;

% Number of lines
m = length(rl);
n = size(lx,2);

% Normalize the line profiles
lim = 1;
lx = lx./mean(abs(lx(abs(rl) < lim,:)));
ly = ly./mean(abs(ly(abs(rl) < lim,:)));
lz = lz./mean(abs(lz(abs(rl) < lim,:)));
dpp = dpp./mean(abs(dpp(abs(rd) < lim,:)));
dpm = dpm./mean(abs(dpm(abs(rd) < lim,:)));
dmm = dmm./mean(abs(dmm(abs(rd) < lim,:)));
dmp = dmp./mean(abs(dmp(abs(rd) < lim,:)));

% Plot the line profiles
figure;
subplot(2,3,1);
plot(rl,abs(lx));
subplot(2,3,2);
plot(rl,abs(ly));
subplot(2,3,3);
plot(rl,abs(lz));
subplot(2,4,5);
plot(rd,abs(dpp));
subplot(2,4,6);
plot(rd,abs(dpm));
subplot(2,4,7);
plot(rd,abs(dmm));
subplot(2,4,8);
plot(rd,abs(dmp));


%% Fit the lines
% Combine xyz
yl = abs([lx ly lz]);
nl = size(yl,2);

% Combine diagonals
yd = abs([dpp dpm dmm dmp]);
nd = size(yd,2);

% Fitting model
fun = @(x0,w,x) erfc((x - x0)./(sqrt(2).*w))./2;

% Fit xyz
lx0 = zeros(nl,1);
lw = lx0;
figure;
for i = 1:nl
    f = fit(rl,yl(:,i),fun,'StartPoint',[2 0.1],'Exclude',rl < 0);
    lx0(i) = f.x0;
    lw(i) = f.w;
    subplot(3,5,i);
    plot(f,rl,yl(:,i));
end

% Fit diagonals
dx0 = zeros(nd,1);
dw = dx0;
figure;
for i = 1:nd
    f = fit(rd,yd(:,i),fun,'StartPoint',[2 0.1],'Exclude',rd < 0);
    dx0(i) = f.x0;
    dw(i) = f.w;
    subplot(4,5,i);
    plot(f,rd,yd(:,i));
end



