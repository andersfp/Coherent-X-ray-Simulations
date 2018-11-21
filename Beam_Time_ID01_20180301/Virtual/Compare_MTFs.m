% Initialization
clear;
close all;
clc;


%% Load data
% Load the no lens data
dn = load('C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\Beam_Time_ID01_20180301\Data_Processing\Grain_G5\No_Lens\Corrected_NoLens_Center.mat');

% Load the virtual geometry center data
dv = load('Corrected_Virtual_Center.mat');

% Load the virtual geometry ptycho data
dp = load('Corrected_Virtual_Ptycho_1.mat');


%% Extract the objects and fields
% No lens object
on = padarray(double(dn.object_avg),[216 184 96])./100;

% Virtual center object
ov = padarray(double(dv.object_avg),[97 79 95])./100;

% Virtual ptycho object
op = double(dp.rho);

% No lens field
fn = fftshift(fftn(fftshift(on)));

% Virtual center field
fv = fftshift(fftn(fftshift(ov)));

% Virtual ptycho field
fp = fftshift(fftn(fftshift(op)));


%% Generate MTFs
% Generate real space axes
rxn = (-256:255).*mean(diff(dn.rx));
ryn = (-256:255).'.*mean(diff(dn.ry));
rzn = permute((-140:139).*mean(diff(dn.rz)),[3 1 2]);
rxv = (-128:127).*mean(diff(dv.rx));
ryv = (-128:127).'.*mean(diff(dv.ry));
rzv = permute((-140:139).*mean(diff(dv.rz)),[3 1 2]);

% Generate reciprocal space axes
qxn = (-256:255)./(512.*mean(diff(dn.rx)));
qyn = (-256:255).'./(512.*mean(diff(dn.ry)));
qzn = permute((-140:139)./(280.*mean(diff(dn.rz))),[3 1 2]);
qxv = (-128:127)./(256.*mean(diff(dv.rx)));
qyv = (-128:127).'./(256.*mean(diff(dv.ry)));
qzv = permute((-140:139)./(280.*mean(diff(dv.rz))),[3 1 2]);

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
    mn(i) = mean(abs(fn(yn == i)));
    mv(i) = mean(abs(fv(yv == i)));
    mp(i) = mean(abs(fp(yv == i)));
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


