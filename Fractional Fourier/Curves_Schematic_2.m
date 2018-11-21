% Initialization
clear;
close all;
clc;

% Save figure?
sav = 0;


%% Set parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42/E;

% Set CRL parameters
R_CRL = 50e-6;
delta = 1.17843774e-6;
f = 0.1*R_CRL/(2*delta);
T = 0.016;
N = 7;
mu = 1/24590.6e-6;
Tweb = 2e-6;

% Calculate focal length
phi = sqrt(T/f);
fN = f*phi*cot(N*phi);

% Set object and image positions (see CRL_Focal_Length_Sim)
M = 2;
d1 = fN*(1 + 1/(M*cos(N*phi)));
d2 = fN*(1 + M/cos(N*phi));

% Calculate analytical apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R_CRL,mu,f,phi,d1);

% Calculate the effective vignetting width
sigma_v = Vignetting(R_CRL,N,mu,d1,T,f,lambda,sigma_p);


%% Make curved surfaces
% Generate distances and focal lengths
D = [d1+T/2;T*ones(N-1,1);T/2+d2];
F = f*ones(N,1);

% Calculate the propagation parameters
m = 100000;
dx = 1e-3;
R0 = -1e10;
s0 = dx/sqrt(m);

% Make a tiny object
x0 = ((-m/2):(m/2 - 1)).'/m*dx;
w = 100e-9;
E0 = exp(-x0.^2./(2*w.^2));

% Optimize the d2 distance
shwplt = 0;
tic;
d2c = optimImagePos(E0,x0,D,F,lambda,R0,s0,sigma_v,sigma_p,shwplt);
toc;
D(end) = D(end) + d2c;

% Calculate the FrFT parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);

% Propagate the object
n = length(D);
Rm2 = zeros(n,1);
Rm2(1) = R0;
for i = 2:n
    [E,x] = propFrFT1(E0,x0,R0,Rm(i),s0,sm(i),sum(a(1:i-1)),lambda,0);
    Rm2(i) = wavefrontCurvature(x,E,lambda,shwplt);
end
Rp2 = zeros(n,1);
for i = 1:n
    [E,x] = propFrFT1(E0,x0,R0,Rp(i),s0,sp(i),sum(a(1:i)),lambda,0);
    Rp2(i) = wavefrontCurvature(x,E,lambda,shwplt);
end

% Generate combined radius and scale
%R = reshape([Rm2,Rp2].',2*length(a),1);
R = reshape([Rm,Rp].',2*length(a),1);
s = reshape([sm,sp].',2*length(a),1);

% Generate position
p = [0;cumsum(D)];
p = reshape([p(1:end-1),p(2:end)].',2*length(a),1);

% Get center of circles
x0 = p - R;

% Get maximum height of circle
ym = s/s0*dx*100;

% Get angles
th = asind(ym./R);

% Plot the curves
c1 = get(groot,'DefaultAxesColorOrder');
bl = c1(1,:);
re = c1(2,:);
col = [0,0,1;repmat([re;bl],length(a)-1,1);1,0,0];
figure;
set(gcf,'Position',[400 300 900 500]);
plot([-10 10],[0 0],'k','LineWidth',1);
hold on;
ys = 0.02;
plot([d1-fN d1-fN],[-ys ys],'k','LineWidth',1);
plot([d1+N*T+fN d1+N*T+fN],[-ys ys],'k','LineWidth',1);
for i = 1:length(R)
    t = linspace(th(i),-th(i),1000);
    x = x0(i) + cosd(t)*R(i);
    y = sind(t)*R(i);
    plot(x,y,'LineWidth',3,'Color',col(i,:));
    %plot(x,y,'LineStyle','none','Marker','o','MarkerSize',3,'MarkerEdgeColor',col(i,:),'MarkerFaceColor',col(i,:));
    %plot(x(1:2),y(1:2),'LineStyle','none','Marker','o','MarkerSize',3,'MarkerEdgeColor',col(i,:),'MarkerFaceColor',col(i,:));
    %plot(x(end-1:end),y(end-1:end),'LineStyle','none','Marker','o','MarkerSize',3,'MarkerEdgeColor',col(i,:),'MarkerFaceColor',col(i,:));
    %scatter(x,y,10,'MarkerEdgeColor',col(i,:),'MarkerFaceColor',col(i,:),'MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2);
    %scatter(x,y,10,'MarkerEdgeColor','none','MarkerFaceColor',col(i,:),'MarkerFaceAlpha',0.09);
end
axis equal;
set(gca,'XLim',[-0.05 1.5],'YLim',[-0.4 0.4],'FontSize',14,'XTick',0:0.2:1.4,'YTick',-0.4:0.2:0.4,'OuterPosition',[-0.15 0 1.3 1]);
xlabel('z [m]');
ylabel('x [m]');

% Save the figure
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('Curves_Schematic_2.pdf','-dpdf','-r0');
    print('Curves_Schematic_2.png','-dpng','-r600');
end



