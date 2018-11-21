% Initialization
clear;
close all;
clc;


%% Load the data
% Load the prepped intensity data
load('Ptycho_Data.mat');

% Get the data size
nx = size(I,2);
ny = size(I,1);
n_omega = size(I,3);
n = size(I,4);


%% Process the data
% Correct the intensity using the aperture
I = I./G;

% Set the pixel offsets
bx = [0.0006;-0.2205;-0.2921;-28.5979;-28.8216;-28.7341;27.0028;26.9619;27.3943];
by = [0.0008;-38.0828;32.4441;-3.9212;-40.6277;30.9477;-4.6579;-41.8149;30.8222];
bz = [0.0004;-0.2330;0.2651;-0.1339;-0.8720;-0.0762;-0.8533;-1.3481;-1.5224];

% Plot the central slice of the first dataset
figure;
imagesc(I(:,:,141,1),[0 1e6]);
axis equal tight;

% Set the ROI on the first image
ix = 127:131;
iy = 127:131;

% Expand the ROI to all 9 datasets
ix = round(ix - bx);
iy = round(iy - by);

% Extract the identical ROIs
A = zeros(n_omega,n);
for i = 1:n
    A(:,i) = squeeze(mean(mean(I(iy(i,:),ix(i,:),:,i))));
end

% Plot the normalized rocking curve
co = parula(n);
figure;
set(gca,'ColorOrder',co,'NextPlot','replacechildren');
semilogy(A);
l = num2cell(1:n);
l = cellfun(@(x) ['Dataset ' num2str(x)],l,'UniformOutput',false);
legend(l);


%% Fit the rocking curves
% Set the exclusions
ex = ones(n_omega,1);
ex(135:149) = 0;

% Set the fitting function
fun = @(a,w,x0,x) a.*exp(-(x - x0).^2./(2.*w.^2));

% Set a starting guess
sg = [max(A).' repmat([3.3 142],n,1)];

% Fit the central rocking curve
a = zeros(n,1);
w = a;
z0 = a;
for i = 1:n
    f = fit((1:n_omega).',A(:,i),fun,'Exclude',ex,'StartPoint',sg(i,:));
    a(i) = f.a;
    w(i) = f.w;
    z0(i) = 141 - f.x0;
    figure;
    plot(f,(1:n_omega).',A(:,i));
end

% Normalize the intensity
a = a./max(a);

% Plot the fitting results
figure;
subplot(1,3,1);
plot(a);
subplot(1,3,2);
plot(w);
subplot(1,3,3);
plot(z0);

% Compare the rocking curve offset
figure;
scatter(z0,bz);
xlabel('Fitted rocking curve');
ylabel('Optimized diffraction pattern position');

% Fit the intensity decrease
f = fit((1:n).',a,'poly3','Exclude',[0 1 0 0 1 0 0 1 0]);
figure;
plot(f,(1:n).',a);
af = f((1:n).');
figure;
plot(1:n,[a af]);

% Save the corrections
save('Ptycho_Corrections.mat','a','af','z0');


