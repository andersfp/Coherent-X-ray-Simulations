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

d(4) = 0;

% Load the images
p = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\MLLs\Figures\';
I1 = imread([p 'I1.png']);
I2 = imread([p 'I2.png']);
I3 = imread([p 'I3.png']);
S1 = imread([p 'S1.png']);
S2 = imread([p 'S2.png']);
S3 = imread([p 'S3.png']);


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

% Plot the corresponding line profiles
ylim = [-2 3];



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
set(gcf,'Position',[300 500 1350 600]);

axes('Position',[0 0.5 2/9 0.5]);
imagesc(I1,[0 255]);
colormap gray;
axis equal tight off;
text(50,120,'a)','FontSize',40,'Color','white');

axes('Position',[2/9 0.5 2/9 0.5]);
imagesc(I2,[0 255]);
colormap gray;
hold on;
rectangle('Position',[145 470 100 50],'EdgeColor',col3,'Linewidth',3);
axis equal tight off;
text(50,120,'b)','FontSize',40,'Color','white');
line(900-50-[180 0],900-[0 0],'LineWidth',5,'Color',col2);
line(900-[0 0],900-50-[124 0],'LineWidth',5,'Color',col2);

axes('Position',[4/9 0.5 2/9 0.5]);
imagesc(I3,[0 255]);
colormap gray;
axis equal tight off;
text(50,120,'c)','FontSize',40,'Color','white');

axes('Position',[0 0 2/9 0.5]);
imagesc(S1,[0 255]);
colormap gray;
axis equal tight off;
text(50,120,'d)','FontSize',40,'Color','white');

axes('Position',[2/9 0 2/9 0.5]);
imagesc(S2,[0 255]);
colormap gray;
axis equal tight off;
text(50,120,'e)','FontSize',40,'Color','white');

axes('Position',[4/9 0 2/9 0.5]);
imagesc(S3,[0 255]);
colormap gray;
axis equal tight off;
text(50,120,'f)','FontSize',40,'Color','white');

axes('OuterPosition',[6/9 0 3/9 1.05]);
plot(xs*(yl),sl + ofs,'Color',col1,'LineWidth',lw);
hold on;
plot(xs*(yl - sdl.'.*dyl),dl(:,idl + k) + ofs,'Color',col3,'LineWidth',lw);
set(gca,'XLim',xlim,'YLim',ylim2,'FontSize',fnt,'XTick',xtick,'YTick',ofs,'YTickLabel',ys*d);
xlabel(['x ' xu]);
ylabel(['\Delta ' yu]);
text(-14,20.6,'g)','FontSize',40);

if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)],'InvertHardcopy','off');
    print('C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\MLLs\Figures\Figure_09.pdf','-dpdf','-r0');
end

