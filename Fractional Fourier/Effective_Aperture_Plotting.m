% Initialization
clear;
close all;
clc;


%% Load data
% Run the simulation scripts?
%run('Effective_Aperture_Small');
%run('Effective_Aperture_Large');

% Save the figures?
s = 0;

% Load the simulation parameters
load('Sim_Parameters.mat');

% Load the small object results
S = load('Effective_Aperture_Small_Profiles.mat');

% Load the large object results
L = load('Effective_Aperture_Large_Profiles.mat');


%% Process data
% Get the coordinate
xs = S.x1;
xl = L.x1;

% Extract the line profiles
ls1 = S.L1;
lsN = S.LN;
ll1 = L.L1;
llN = L.LN;

% Get the intensity
is1 = abs(ls1).^2;
isN = abs(lsN).^2;
il1 = abs(ll1).^2;
ilN = abs(llN).^2;

% Get the phase
R = -0.2946;
ps1 = unwrap(angle(ls1.*exp(-1i*pi*xs.^2/(lambda*R)))) + pi*xs.^2/(lambda*R);
psN = unwrap(angle(lsN.*exp(-1i*pi*xs.^2/(lambda*R)))) + pi*xs.^2/(lambda*R);
pl1 = unwrap(angle(ll1.*exp(-1i*pi*xl.^2/(lambda*R)))) + pi*xl.^2/(lambda*R);
plN = unwrap(angle(llN.*exp(-1i*pi*xl.^2/(lambda*R)))) + pi*xl.^2/(lambda*R);
ps1 = ps1 - ps1(xs == 0);
psN = psN - psN(xs == 0);
pl1 = pl1 - pl1(xl == 0);
plN = plN - plN(xl == 0);

% Extract the effective attenuation
th = 1e-9;
as = isN./is1;
as(is1/max(is1) < th) = NaN;
al = ilN./il1;
al(il1/max(il1) < th) = NaN;

% Calculate the theoretical small object attenuation profile
y0 = sigA*d0;
a0 = sigA;
yn = y0*(cos((1:N)*phi) + phi/2*sin((1:N)*phi)) + a0*(f*phi*sin((1:N)*phi) - T/2*cos((1:N)*phi));
yN = yn(end);
at = exp(-xs.^2/(2*yN.^2));

% Calculate vignetting parameters
gamma = ((f + 2*d0*N)*sqrt(d0^2 + (f*phi)^2)*cos(2*N*phi + atan(d0/(f*phi))))/(2*(d0*f + N*(d0^2 + f^2)) + (d0^2 + (f*phi)^2)*sin(2*N*phi)/phi);
sigV = delta/(mu*sigA*sqrt((N*phi)^2 - sin(N*phi)^2));

% Small object parameter check
xmaxs = sqrt(2)*1.2e-7/2;
xmaxl = sqrt(2)*600e-6/2;
disp([abs(gamma)*xmaxs/sigA (xmaxs/sigV)^2]);
disp([abs(gamma)*xmaxl/sigA (xmaxl/sigV)^2]);

% Fit the effective aperture widths
fun = @(sig,x) exp(-x.^2./(2*sig.^2));
fs = fit(xs,as,fun,'StartPoint',yN);
fl = fit(xl,al,fun,'StartPoint',yN);


%% Make plots
% Plot the small object intensities
figure;
plot(1e6*xs,[is1 isN]);
title('Small object intensity');
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
legend('No attenuation','Individual lens attenuation');

% Plot the small object phases
figure;
plot(1e6*xs,[ps1 psN]);
title('Small object phase');
xlabel('Position (\mum)');
ylabel('Phase (rad)');
legend('No attenuation','Individual lens attenuation');

% Plot the small object phase differences
figure;
plot(1e6*xs,ps1 - psN);
title('Small object phase difference');
xlabel('Position (\mum)');
ylabel('Phase (rad)');

% Plot the large object intensities
figure;
plot(1e6*xl,[il1 ilN]);
title('Large object intensity');
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
legend('No attenuation','Individual lens attenuation');

% Plot the large object phases
figure;
plot(1e6*xl,[pl1 plN]);
title('Large object phase');
xlabel('Position (\mum)');
ylabel('Phase (rad)');
legend('No attenuation','Individual lens attenuation');

% Plot the large object phase differences
figure;
plot(1e6*xl,pl1 - plN);
title('Large object phase difference');
xlabel('Position (\mum)');
ylabel('Phase (rad)');

% Plot the attenuations
figure;
plot(1e6*xs,as,1e6*xl,al,1e6*xs,at);
title('Attenuation profiles');
xlabel('Position (\mum)');
ylabel('Attenuation');
legend('Small object','Large object','Theoretical small object');

% Plot the image plane for small object
figure;
plot(1e6*S.x2,abs(S.L2).^2,1e6*S.xM,abs(S.LM).^2);
title('Small object image intensity');
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
legend('Effective aperture','Lens by lens');

% Plot the image phase for small object
figure;
plot(1e6*S.x2,unwrap(angle(S.L2)),1e6*S.xM,unwrap(angle(S.LM)));
title('Small object image phase');
xlabel('Position (\mum)');
ylabel('Phase (rad)');
legend('Effective aperture','Lens by lens');

% Plot the image plane for large object
figure;
plot(1e6*L.x2,abs(L.L2).^2,1e6*L.xM,abs(L.LM).^2);
title('Large object image intensity');
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
legend('Effective aperture','Lens by lens');

% Plot the image phase for large object
figure;
plot(1e6*L.x2,unwrap(angle(L.L2)),1e6*L.xM,unwrap(angle(L.LM)));
title('Large object image phase');
xlabel('Position (\mum)');
ylabel('Phase (rad)');
legend('Effective aperture','Lens by lens');


%% Plots for paper
% Plotting parameters
F = 14;
lw = 2;

% Plot the attenuations
figure;
plot(1e6*xs,as,'-',1e6*xs,at,'--','LineWidth',lw);
xlabel('Position (\mum)');
ylabel('Relative intensity');
legend('Small object','Theoretical');
set(gca,'FontSize',F,'XLim',[-300 300],'YLim',[0 1.1]);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Effective_Aperture.pdf','-dpdf','-r0');
    print('Effective_Aperture.png','-dpng','-r600');
end

% Plot the objects
figure;
set(gcf,'Position',[150 370 1200 500]);
subplot(1,2,1);
yyaxis left;
plot(1e6*S.x2,abs(S.L2).^2,'LineWidth',lw);
hold on;
plot(1e6*S.xM,abs(S.LM).^2,'--','LineWidth',lw,'Color',[0.9290 0.6940 0.1250]);
set(gca,'FontSize',F,'XLim',[-3 3],'OuterPosition',[0 0 0.5 1]);
xlabel('Position (\mum)');
ylabel('Intensity (a.u.)');
text(-2,0.009,'(a)','FontSize',22,'Color','black');
yyaxis right;
plot(1e6*S.x2,unwrap(angle(S.L2))-1015+0.2791,'LineWidth',lw);
hold on;
plot(1e6*S.xM,unwrap(angle(S.LM))-1065+0.0135,'--','LineWidth',lw,'Color',[0.9290 0.6940 0.1250]);
set(gca,'YLim',[-4 8]);
ylabel('Phase (rad)');
subplot(1,2,2);
yyaxis left;
plot(1e3*L.x2,abs(L.L2).^2,'LineWidth',lw);
hold on;
plot(1e3*L.xM,abs(L.LM).^2,'--','LineWidth',lw,'Color',[0.9290 0.6940 0.1250]);
set(gca,'FontSize',F,'XLim',[-4 4],'OuterPosition',[0.5 0 0.5 1]);
xlabel('Position (mm)');
ylabel('Intensity (a.u.)');
text(-3,0.009,'(b)','FontSize',22,'Color','black');
yyaxis right;
plot(1e3*L.x2,unwrap(angle(L.L2)),'LineWidth',lw);
hold on;
plot(1e3*L.xM,unwrap(angle(L.LM)),'--','LineWidth',lw,'Color',[0.9290 0.6940 0.1250]);
ylabel('Phase (rad)');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Object_Image_Effective.pdf','-dpdf','-r0');
    print('Object_Image_Effective.png','-dpng','-r600');
end

% Plot the differences
figure;
set(gcf,'Position',[150 370 1200 500]);
subplot(1,2,1);
yyaxis left;
plot(1e6*S.x2,100*(abs(S.L2).^2 - abs(S.LM).^2)./abs(S.L2).^2,'LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0 0 0.5 1],'XLim',[-3 3],'YLim',[-1 1]);
xlabel('Position (\mum)');
ylabel('Relative intensity difference (%)');
text(-2.8,0.9,'(a)','FontSize',22,'Color','black');
yyaxis right;
plot(1e6*S.x2,unwrap(angle(S.L2)) - unwrap(angle(S.LM)) + 50.27 - 0.004416,'LineWidth',lw);
set(gca,'YLim',[-0.06 0.06]);
ylabel('Phase difference (rad)');
subplot(1,2,2);
yyaxis left;
plot(1e3*L.x2,100*(abs(L.L2).^2 - abs(L.LM).^2)./abs(L.L2).^2,'LineWidth',lw);
set(gca,'FontSize',F,'OuterPosition',[0.5 0 0.5 1],'XLim',[-5 5],'YLim',[-1 1],'XTick',-4:2:4);
xlabel('Position (mm)');
ylabel('Relative intensity difference (%)');
text(-4.7,0.9,'(b)','FontSize',22,'Color','black');
yyaxis right;
plot(1e3*L.x2,unwrap(angle(L.L2)) - unwrap(angle(L.LM)),'LineWidth',lw);
set(gca,'YLim',[-0.06 0.06]);
ylabel('Phase difference (rad)');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if s == 1
    print('Object_Image_Effective_Difference.pdf','-dpdf','-r0');
    print('Object_Image_Effective_Difference.png','-dpng','-r600');
end


