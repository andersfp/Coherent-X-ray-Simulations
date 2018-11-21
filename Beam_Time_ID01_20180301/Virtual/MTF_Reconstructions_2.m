% Initialization
clear;
close all;
clc;


%% Load the data
% Load the no lens reconstruction
nl = load('Reconstruction_Center_2018_04_17_16_58_1161.mat');

% Load the ptycho 1 reconstruction
p1 = load('Reconstruction_Ptycho_1_2018_07_23_14_12_129.mat');

% Load the ptycho 5 reconstruction
p5 = load('Reconstruction_Ptycho_5_2018_07_20_18_49_3809.mat');


%% Process the data
% Generate the fields
F0 = fftshift(fftn(fftshift(padarray(nl.object_avg./nl.cyc,[97 79 95]))));
F1 = fftshift(fftn(fftshift(p1.rho)));
F5 = fftshift(fftn(fftshift(p5.rho)));

% Plot the fields
rng = [-1 3];
Slicer(log10(abs(F0)),'displayRange',rng);
Slicer(log10(abs(F1)),'displayRange',rng);
Slicer(log10(abs(F5)),'displayRange',rng);

% Set reciprocal space pixel size
delta_q1 = 1.8915e6;
delta_q2 = 1.8915e6;
delta_q3 = 1.4521e6;

% Generate reciprocal space axes
q1 = (127:-1:-128).'.*delta_q1;
q2 = (-128:127).*delta_q2;
q3 = permute((-140:139).*delta_q3,[3 1 2]);

% Calculate orthogonal reciprocal space coordinates
th = 20;
qx = q2;
qy = cosd(th).*q1;
qz = q3 - sind(th).*q1;

% Calculate the magnitude of Q
Q = sqrt(qx.^2 + qy.^2);

% Calculate the intensities
I0 = sum(abs(F0).^2,3);
I1 = sum(abs(F1).^2,3);
I5 = sum(abs(F5).^2,3);

% Normalize the intensities
I0 = I0./max(I0(:));
I1 = I1./max(I1(:));
I5 = I5./max(I5(:));

% Discretize the intensity
n = 100;
tic;
[y,e] = discretize(Q(:),n);
toc;
m0 = zeros(n,1);
m1 = m0;
m5 = m0;
tic;
for i = 1:n
    m0(i) = mean(I0(y == i));
    m1(i) = mean(I1(y == i));
    m5(i) = mean(I5(y == i));
end
toc;

% Plot the distribution
e = (e(1:end-1) + e(2:end))/2;
figure;
semilogy(e,m0,e,m1,e,m5);
legend('Virtual','Ptycho 1','Ptycho 5');

% Calculate fitting friendly parameters
X = 1e-9.*e.';
Y0 = log10(m0);
Y1 = log10(m1);
Y5 = log10(m5);

% Fit the distributions
fun = @(a,b,c,x) a.*exp(-b.*x) + c;
f0 = fit(X,Y0,fun,'StartPoint',[8.6 10.4 -8.7]);
f1 = fit(X,Y1,fun,'StartPoint',[7.6 10.6 -7.7]);
f5 = fit(X,Y5,fun,'StartPoint',[6.6 10.4 -6.9]);

% Plot the fits
figure;
plot(f0,X,Y0);
figure;
plot(f1,X,Y1);
figure;
plot(f5,X,Y5);

%% Save figure
% Figure for saving
s = 1e-9;
figure;
semilogy(s.*e,m0,s.*e,m1,s.*e,m5,'LineWidth',3);
legend('Single data set','Parallel synthesis','Serial synthesis');
xlabel('|\bfq\rm_{\perp}| (nm^{-1})');
ylabel('Normalized intensity');
set(gca,'FontSize',18,'Xlim',[0 0.33],'XTick',0:0.05:0.3,'XTickLabel',{0,[],0.1,[],0.2,[],0.3},'YTick',10.^(-8:2:0),'YMinorTick','off');
line([0 0.07747],[1e-3 1e-3],'LineStyle','--','Color','k','LineWidth',2);
line([1 1]*0.047195,[1e-3 1e-8],'LineStyle','--','Color','k','LineWidth',2);
line([1 1]*0.06415,[1e-3 1e-8],'LineStyle','--','Color','k','LineWidth',2);
line([1 1]*0.07747,[1e-3 1e-8],'LineStyle','--','Color','k','LineWidth',2);
%print('MTF_Virtual_Synthesis.png','-dpng','-r192');


