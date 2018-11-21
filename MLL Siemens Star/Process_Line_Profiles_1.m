% Initialization
clear;
close all;
clc;


%% Load the data
% Load line profiles
load('Line_Profiles_1.mat');

% Get data sizes
mh = length(x);
mv = length(y);
nd = size(dh,2);
ns = size(sh,2);


%% Process
% Normalize the simulated horizontal profiles
sh = sh - mean(sh(:));
sh = sh./std(sh(:));

% Normalize the simulated vertical profiles
sv = sv - mean(sv(:));
sv = sv./std(sv(:));

% Normalize the measured profiles
% ih1 = 1;
% ih2 = mh;
% iv1 = 1;
% iv2 = mv;
% for i = 1:nd
%     dh(:,i) = dh(:,i) - mean(dh(ih1:ih2,i));
%     dv(:,i) = dv(:,i) - mean(dv(iv1:iv2,i));
%     dh(:,i) = dh(:,i)./std(dh(ih1:ih2,i));
%     dv(:,i) = dv(:,i)./std(dv(iv1:iv2,i));
% end
dh = dh - mean(dh(:));
dv = dv - mean(dv(:));
dh = dh./std(dh(:));
dv = dv./std(dv(:));

% Plot the horizontal data
figure;
plotInt2(x,dh,(x + 14.066e-6)/1.1172,sh(:,11:51));
axis tight;
set(gca,'YLim',[-4 4]);

% Plot the vertical data
figure;
plotInt2(y,dv,(y - 2.1291e-5)*1.1066,sv(:,23:63));
axis tight;
set(gca,'YLim',[-4 4]);


%% Plot the line profiles
% Plot settings
ofs = 3;
lw = 2;
r = [0.8500 0.3250 0.0980]; % [215,25,28]/255;
b = [0 0.4470 0.7410]; % [44,123,182]/255;
z = ((1:41).' - 19)*0.1;
F = 14;

% Plot the horizontal line profiles
ih = 6:5:41;
ofsh = ofs*repmat(1:length(ih),mh,1);
figure;
plot(1e6*x,dh(:,ih) + ofsh,'-','LineWidth',lw,'Color',b);
hold on;
plot(1e6*(x + 14.066e-6)/1.1172,sh(:,16:5:51) + ofsh,'-','LineWidth',lw,'Color',r);
set(gca,'XLim',[-90 50],'YLim',[0 9*ofs],'XTick',-80:20:40,'YTick',(1:8)*ofs,'YTickLabel',num2str(z(ih),'%.1f'),'FontSize',F);
xlabel('x [\mum]');
ylabel('z [mm]');

% Plot the vertical line profiles
iv = 4:5:41;
ofsv = ofs*repmat(1:length(iv),mv,1);
figure;
plot(1e6*y,dv(:,iv) + ofsv,'-','LineWidth',lw,'Color',b);
hold on;
plot(1e6*(y - 2.1291e-5)*1.1066,sv(:,26:5:63) + ofsv,'-','LineWidth',lw,'Color',r);
set(gca,'XLim',[-110 50],'YLim',[0 9*ofs],'XTick',-100:20:40,'YTick',(1:8)*ofs,'YTickLabel',num2str(z(iv),'%.1f'),'FontSize',F);
xlabel('y [\mum]');
ylabel('z [mm]');


