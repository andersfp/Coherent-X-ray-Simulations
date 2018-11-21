% Initialization
clear;
close all;
clc;


%% Load data
% Load the no lens data
dn = load('C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\Beam_Time_ID01_20180301\Data_Processing\Grain_G5\No_Lens\Reconstruction_Center_2018_03_13_23_38_300.mat');

% Load the virtual geometry center data
dv = load('Reconstruction_Center_2018_04_17_16_58_1161.mat');

% Load the virtual geometry ptycho data
dp = load('Reconstruction_Ptycho_1_2018_06_21_16_07_20.mat');


%% Extract the objects and fields
% No lens object
on = padarray(double(dn.object_avg),[216 184 96])./dn.cyc;

% Virtual center object
ov = padarray(double(dv.object_avg),[97 79 95])./dv.cyc;

% Virtual ptycho object
op = double(dp.rho);

% No lens field
fn = fftshift(fftn(fftshift(on)));

% Virtual center field
fv = fftshift(fftn(fftshift(ov)));

% Virtual ptycho field
fp = fftshift(fftn(fftshift(op)));


%% Generate MTFs
% Set up the experimental parameters
Dn = 1.632;
Dv = 1.662/1.48;
det_dx = 55e-6;
E = 8e3;
lambda = E2lambda(E);
k = 2.*pi./lambda;
delta_omega = 0.003*pi/180;
a = 3.9242e-10;
q0 = 2*pi/a*sqrt(3);
th = asin(q0.*lambda./(4.*pi));

% Calculate pixel size
[delta_q1n,delta_q2n,delta_q3n] = q_range2(Dn,det_dx,th,lambda,delta_omega,512,512,280);
[delta_q1v,delta_q2v,delta_q3v] = q_range2(Dv,det_dx,th,lambda,delta_omega,256,256,280);

% Generate intrinsic reciprocal space coordinates
q1n = (-203:308).'.*delta_q1n;
q2n = (-180:331).*delta_q2n;
q3n = permute((-141:138).*delta_q3n,[3 1 2]);
q1v = (-127:128).'.*delta_q1v;
q2v = (-127:128).*delta_q2v;
q3v = permute((-141:138).*delta_q3v,[3 1 2]);

% Generate orthgonal reciprocal space coordinates
qxn = q2n;
qyn = q1n - sin(th).*q3n;
qzn = cos(th).*q3n;
qxv = q2v;
qyv = q1v - sin(th).*q3v;
qzv = cos(th).*q3v;

% Calculate the magnitude of q
qn = sqrt(qxn.^2 + qyn.^2 + qzn.^2);
qv = sqrt(qxv.^2 + qyv.^2 + qzv.^2);

% Plot the MTFs
n = 100;
tic;
[yn,en] = discretize(qn(:),n);
[yv,ev] = discretize(qv(:),n);
toc;
mn = zeros(n,1);
mv = mn;
mp = mn;
tic;
for i = 1:n
    mn(i) = mean(abs(fn(yn == i)).^2);
    mv(i) = mean(abs(fv(yv == i)).^2);
    mp(i) = mean(abs(fp(yv == i)).^2);
end
toc;

% mn = mn./mean(mn(1:6));
% mv = mv./mean(mv(1:8));
% mp = mp./mean(mp(1:8));

en = (en(1:end-1) + en(2:end))/2;
ev = (ev(1:end-1) + ev(2:end))/2;

figure;
semilogy(en,mn,ev,mv,ev,mp);
legend('No lens','Virtual center','Virtual ptycho');


