% Initialization
clear;
close all;
clc;


%% Load data
% Load the diffraction patterns
tic;
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';
I = h5read([p 'Scan0008.cxi'],'/entry_1/instrument_1/detector_1/data');
L = h5read([p 'Scan0039.cxi'],'/entry_1/instrument_1/detector_1/data');
toc;

% Convert data
I = double(permute(I,[2 1 3]));
L = double(permute(L,[2 1 3]));

% Load the mask
load('MaxiPix_Mask.mat');
mask = double(mask);
mask(mask == 0) = NaN;

% Apply the mask
I = mask.*I;
L = mask.*L;

% Load the second mask
mask2 = imread('Mask_Lens.png');
mask2 = double(mask2);
mask2(mask2 < 127) = NaN;
mask2(mask2 > 127) = 1;

% Apply the second mask
L = mask2.*L;


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
imagesc(L(:,:,142));
axis equal tight;
subplot(2,2,3);
imagesc(log10(I(:,:,142)));
axis equal tight;
subplot(2,2,4);
imagesc(log10(L(:,:,142)));
axis equal tight;

% Extract ROIs
A = I(106:305,83:282,133);
B = L(106:305,83:282,133);

% Plot the ROIs
figure;
subplot(1,2,1);
imagesc(log10(A),[0 5]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[0 5]);
axis equal tight;

% Fit the magnification and pupil
%fun2 = @(x) sum(sum((log10(fun(A,x(1),x(2),x(3),x(4),x(5),x(6))+1) - log10(B+1)).^2,'omitnan'),'omitnan');
fun2 = @(x) sum(sum(((fun(A,x(1),x(2),x(3),0,0,x(4),x(5))) - (B)).^2,'omitnan'),'omitnan');
x0 = [1.38 -1.70 2.66 27 1.5]; % [1.11 0.86 0.08 -1.37 -0.15 5.42]
opts = optimset('MaxFunEvals',10000,'MaxIter',10000,'PlotFcns',@optimplotfval,'TolFun',1e6);
xs = linspace(1.38,1.38,1);
fval = zeros(length(xs),1);
xmin = zeros(length(xs),5);
for i = 1:length(xs)
    x0(1) = xs(i);
    [xmin(i,:),fval(i),exitflag,output] = fminsearch(fun2,x0,opts);
    fprintf('.');
end
fprintf('\n');

figure;
plot(fval);

figure;
plot(fval,xmin);



% Plot the fitting result
A3 = fun(A,xmin(1,1),xmin(1,2),xmin(1,3),0,0,xmin(1,4),xmin(1,5));
figure;
subplot(1,2,1);
imagesc(log10(A3),[0 4]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[0 4]);
axis equal tight;

% Plot an overlay
rgb = cat(3,log10(A3)./5,log10(B)./5,log10(A3)./5);
figure;
image(rgb);
axis equal tight;

% Calculate the width based on sigma_a
N = 20;
T = 2e-3;
R = 50e-6;
E = 8e3;
[delta,mu] = Be_Prop(E);
d1 = 0.0499;
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
A2 = interp2(x,y,A,s.*(x + dx),s.*(y + dy),'linear',NaN);

% Apply the pupil function
A3 = A2.*gaussRMS(x - x0,w).*gaussRMS(y - y0,w);

% Apply an overall scale factor
A3 = sf.*A3;

end


