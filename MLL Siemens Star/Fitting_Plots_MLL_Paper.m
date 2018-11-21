% Initialization
clear;
close all;
clc;

% Save figures?
sav = 0;


%% Load the data
% Load the measurement ROI
load('ROI_Measurement_1.mat');

% Load the simulation ROI
load('ROI_Simulation_1_Partial_Avg.mat');


%% Generate line profiles
% Make line profiles
st = squeeze(mean(St,1));
sl = squeeze(mean(Sl,2));
sb = flip(st,1);
dt = squeeze(mean(Dt,1));
dl = squeeze(mean(Dl,2));
db = squeeze(mean(Db,1));

% Get signal levels of the line profiles
ist = 4;
isl = 4;
isb = 4;
idt = 28;
idl = 19;
idb = 29;
hst = mean(st(xt < -10e-6,ist));
lst = mean(st(xt > 10e-6,ist));
hsl = mean(sl(yl > 10e-6,isl));
lsl = mean(sl(yl < -10e-6,isl));
hsb = mean(sb(xb > 10e-6,isb));
lsb = mean(sb(xb < -10e-6,isb));
hdt = mean(dt(xt < -10e-6,idt));
ldt = mean(dt(xt > 10e-6,idt));
hdl = mean(dl(yl > 10e-6,idl));
ldl = mean(dl(yl < -10e-6,idl));
hdb = mean(db(xb > 12e-6,idb));
ldb = mean(db(xb < -10e-6,idb));

% Normalize the line profiles
st = (st - lst)./(hst - lst);
sl = (sl - lsl)./(hsl - lsl);
sb = (sb - lsb)./(hsb - lsb);
dt = (dt - ldt)./(hdt - ldt);
dl = (dl - ldl)./(hdl - ldl);
db = (db - ldb)./(hdb - ldb);


%% Fit the line profiles
% Set up the fitting model
fun = @(a,b,d,sig,x) a.*0.5.*erfc((x - b)./(sqrt(2).*sig)) + d;

% Set the magnifications
Mx = 61.4;
My = 46.0;

% Fit the data
xt = 1e9*xt;
ft = fit(xt.'/Mx,dt(:,idt),fun,'StartPoint',[1 0 0 30]);
et = diff(confint(ft));
yl = 1e9*yl;
fl = fit(-yl/My,dl(:,idl),fun,'StartPoint',[1 0 0 100]);
el = diff(confint(fl));

% Plot the data
F = 20;
lw = 3;
ms = 20;
figure;
plot(xt.'/Mx,ft(xt.'/Mx),'-',xt.'/Mx,dt(:,idt),'.','LineWidth',lw,'MarkerSize',ms);
set(gca,'XLim',[-400 400],'YLim',[-0.1 1.3],'FontSize',F);
text(-350,0.2,{['\sigma = ' num2str(ft.sig,'%.1f') ' nm'],['\pm' num2str(et(end)/2,'%.1f') ' nm']},'FontSize',F);
xlabel('x [nm]');
ylabel('Normalized intensity');
if sav == 1
    print('C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\MLLs\Figures\Edge_Fit_Top.png','-dpng','-r600');
end

figure;
plot(yl/My,fl(yl/My),'-',-yl/My,dl(:,idl),'.','LineWidth',lw,'MarkerSize',ms);
set(gca,'XLim',[-400 400],'YLim',[-0.1 1.3],'FontSize',F);
text(-350,0.2,{['\sigma = ' num2str(fl.sig,'%.1f') ' nm'],['\pm' num2str(el(end)/2,'%.1f') ' nm']},'FontSize',F);
xlabel('y [nm]');
ylabel('Normalized intensity');
if sav == 1
    print('C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\MLLs\Figures\Edge_Fit_Left.png','-dpng','-r600');
end



