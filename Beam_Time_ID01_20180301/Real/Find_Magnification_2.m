% Initialization
clear;
close all;
clc;


%% Load data
% Load the diffraction patterns
tic;
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';
I = h5read([p 'Scan0008.cxi'],'/entry_1/instrument_1/detector_1/data');
L = h5read([p 'Scan0267.cxi'],'/entry_1/instrument_1/detector_1/data');
toc;

% Convert data
I = double(permute(I,[2 1 3]));
L = double(permute(L,[2 1 3]));


%% Compare data
% Plot 3D image of the two data sets
Slicer(log10(I),'displayRange',[0 5]);
Slicer(log10(L),'displayRange',[0 5]);

% Plot the central image
figure;
subplot(2,2,1);
imagesc(I(:,:,142));
axis equal tight;
subplot(2,2,2);
imagesc(L(:,:,51));
axis equal tight;
subplot(2,2,3);
imagesc(log10(I(:,:,142)));
axis equal tight;
subplot(2,2,4);
imagesc(log10(L(:,:,51)));
axis equal tight;

% Extract ROIs
A = I(156:255,133:232,131);
B = flip(flip(L(156:255,133:232,40),1),2);

% Plot the ROIs
figure;
subplot(1,2,1);
imagesc(log10(A),[0 4]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[0 4]);
axis equal tight;

% Fit the magnification and pupil
fun2 = @(x) sum(sum((fun(A,x(1),x(2),x(3),0,0,x(4),x(5)) - B).^2));
x0 = [1.15 1.3 0.2 5.7 0.9]; % [1.11 0.86 0.08 -1.37 -0.15 5.42]
opts = optimset('MaxFunEvals',10000,'MaxIter',10000,'PlotFcns',@optimplotfval,'TolFun',1e3);
[xmin,fval,exitflag,output] = fminsearch(fun2,x0,opts);

% Plot the fitting result
A3 = fun(A,xmin(1),xmin(2),xmin(3),0,0,xmin(4),xmin(5));
figure;
subplot(1,2,1);
imagesc(log10(A3),[0 4]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[0 4]);
axis equal tight;

% Plot an overlay
rgb = cat(3,log10(A3)./4,log10(B)./4,log10(A3)./4);
figure;
image(rgb);
axis equal tight;

% Calculate the width based on sigma_a
N = 70;
T = 2e-3;
R = 50e-6;
E = 8e3;
[delta,mu] = Be_Prop(E);
d1 = 0.121;
d2 = 1.632 - N*T - d1;
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);
[sigma_D,sigma_a,sigma_v,gamma,yN] = CRL_Parameters_2(N,R,mu,f,phi,d1);
M11 = cos(N.*phi);
M12 = f.*phi.*sin(N.*phi);
M21 = -sin(N.*phi)./(f.*phi);
M22 = cos(N.*phi);
K12 = M12 + d1.*(M11 + d2.*M21) + d2.*M22;
sigma_det = K12.*sigma_a;
sigma_pix = sigma_det./55e-6;


% Define optimization function
function A3 = fun(A,s,dx,dy,x0,y0,w,sf)
% Get the number of pixels
nx = size(A,2);
ny = size(A,1);

% Make axes
x = (-nx/2):(nx/2 - 1);
y = ((-ny/2):(ny/2 - 1)).';

% Interpolate the reference measurement A
A2 = interp2(x,y,A,s.*(x + dx),s.*(y + dy),'linear',0);

% Apply the pupil function
A3 = A2.*gaussRMS(x - x0,w).*gaussRMS(y - y0,w);

% Apply an overall scale factor
A3 = sf.*A3;

end


