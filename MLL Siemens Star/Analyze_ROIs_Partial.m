% Initialization
clear;
close all;
clc;


%% Load the data
% Load the measurement ROI
load('ROI_Measurement_1.mat');

% Load the simulation ROI
load('ROI_Simulation_1_Partial_Avg.mat');

d(4) = 0;


%% Generate line profiles
% Plot the ROIs
Slicer(St);
Slicer(Sl);
Slicer(Dt);
Slicer(Dl);
Slicer(Db);

% Make line profiles
st = squeeze(mean(St,1));
sl = squeeze(mean(Sl,2));
sb = flip(st,1);
dt = squeeze(mean(Dt,1));
dl = squeeze(mean(Dl,2));
db = squeeze(mean(Db,1));

% Plot the line profiles
figure;
plotInt(xt,st);
figure;
plotInt(yl,sl);
figure;
plotInt(xb,sb);
figure;
plotInt(xt,dt);
figure;
plotInt(yl,dl);
figure;
plotInt(xb,db);

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

% Plot the corresponding line profiles
m = 9;
ylim = [-2 3];
figure;
plotInt2(xt,st,xt,dt(:,idt-m:3:idt+m));
set(gca,'YLim',ylim);
figure;
plotInt2(yl,sl,yl,dl(:,idl-m:3:idl+m));
set(gca,'YLim',ylim);
figure;
plotInt2(xb,sb,xb,db(:,idb-m:3:idb+m));
set(gca,'YLim',ylim);


%% Process and plot the line profiles
% Align the line profiles
n = size(sb,2);
sdt = zeros(n,1);
sdl = sdt;
sdb = sdt;
for i = 1:n
    sdt(i) = finddelay(st(:,i),dt(:,idt + (i - 4)*3));
    sdl(i) = finddelay(sl(:,i),dl(:,idl + (i - 4)*3));
    sdb(i) = finddelay(sb(:,i),db(:,idb + (i - 4)*3));
end

% Plot settings
dxt = mean(diff(xt));
dyl = mean(diff(yl));
dxb = mean(diff(xb));
k = -9:3:9;
s = 3;
ofs = 0:s:(s*(length(k) - 1));
lw = 2;
col1 = [0 0.4470 0.7410];
col2 = [0.8500 0.3250 0.0980];
col3 = [0.9290 0.6940 0.1250];
col4 = [0.4940 0.1840 0.5560];
xlim = [-15 15];
ylim2 = [ylim(1) max(ofs)+4];
fnt = 16;
xs = 1e6;
xu = '[\mum]';
xtick = -10:10:10;
ys = 1e3;
yu = '[mm]';

% Plot the aligned profiles
figure;
plot(xs*(xt.'),st + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(xt.' - sdt.'.*dxt),dt(:,idt + k) + ofs,'Color',col2,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);
figure;
plot(xs*(yl),sl + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(yl - sdl.'.*dyl),dl(:,idl + k) + ofs,'Color',col2,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);
figure;
plot(xs*(xb.'),sb + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(xb.' - sdb.'.*dxb),db(:,idb + k) + ofs,'Color',col2,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);

% Plot the aligned profiles in a single figure
figure;
subplot(1,3,1);
plot(xs*(xt.'),st + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(xt.' - sdt.'.*dxt),dt(:,idt + k) + ofs,'Color',col2,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
set(gca,'OuterPosition',[0 0 1/3 1]);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);
title('Top');
subplot(1,3,2);
plot(xs*(yl),sl + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(yl - sdl.'.*dyl),dl(:,idl + k) + ofs,'Color',col3,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
set(gca,'OuterPosition',[1/3 0 1/3 1]);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);
title('Left');
subplot(1,3,3);
plot(xs*(xb.'),sb + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(xb.' - sdb.'.*dxb),db(:,idb + k) + ofs,'Color',col4,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
set(gca,'OuterPosition',[2/3 0 1/3 1]);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);
title('Bottom');


%% Fit the line profiles
% Set up the fitting model
sigs = 1e-7;
bs = 1e-7;
fun = @(a,b,d,sig,x) a.*0.5.*erfc((x - bs.*b)./(sqrt(2).*sigs.*sig)) + d;

% Set the magnifications
Mx = 61.4;
My = 46.0;

% Fit the top data
m2 = 2;
sig = zeros(2*m2+1,3);
sige = sig;
figure;
for i = -m2:m2
    ft = fit(xt.'/Mx,dt(:,idt+i),fun,'StartPoint',[1 0 0 0.5]);
    et = diff(confint(ft));
    subplot(3,2*m2+1,sub2ind([2*m2+1 3],i+m2+1,1));
    plot(ft,xt.'/Mx,dt(:,idt+i));
    set(gca,'XLim',[-4e-7 4e-7],'YLim',[-0.5 1.6]);
    text(-3e-7,0,{['\sigma = ' num2str(ft.sig*sigs*1e9,'%.1f') ' nm'],['\pm' num2str(et(end)*sigs*1e9/2,'%.1f') ' nm']});
    title(['Top: ' num2str(i/10,'%+.1f') ' mm']);
    fl = fit(-yl/My,dl(:,idl+i),fun,'StartPoint',[1 0 0 0.5]);
    el = diff(confint(fl));
    subplot(3,2*m2+1,sub2ind([2*m2+1 3],i+m2+1,2));
    plot(fl,-yl/My,dl(:,idl+i));
    set(gca,'XLim',[-4e-7 4e-7],'YLim',[-0.5 1.6]);
    text(-3e-7,0,{['\sigma = ' num2str(fl.sig*sigs*1e9,'%.1f') ' nm'],['\pm' num2str(el(end)*sigs*1e9/2,'%.1f') ' nm']});
    title(['Left: ' num2str(i/10,'%+.1f') ' mm']);
    fb = fit(-xb.'/Mx,db(:,idb+i),fun,'StartPoint',[1 0 0 0.5]);
    eb = diff(confint(fb));
    subplot(3,2*m2+1,sub2ind([2*m2+1 3],i+m2+1,3));
    plot(fb,-xb.'/Mx,db(:,idb+i));
    set(gca,'XLim',[-4e-7 4e-7],'YLim',[-0.5 1.6]);
    text(-3e-7,0,{['\sigma = ' num2str(fb.sig*sigs*1e9,'%.1f') ' nm'],['\pm' num2str(eb(end)*sigs*1e9/2,'%.1f') ' nm']});
    title(['Bottom: ' num2str(i/10,'%+.1f') ' mm']);
    sig(i+m2+1,:) = [ft.sig fl.sig fb.sig];
    sige(i+m2+1,:) = [et(end) el(end) eb(end)]/2;
end
sig = sigs.*sig;
sige = sigs.*sige;


