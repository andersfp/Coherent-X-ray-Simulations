% Initialization
clear;
close all;
clc;


%% Load data
% Set the data path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_4x\';

% Load the intensity data
load([p 'Simulated_Intensity.mat']);

% Load probe and sample data
load([p 'Phantom.mat'],'obj','rx','ry','rz','fov','x','y','z');
load([p 'Exit_Fields.mat'],'E','lambda','np2d','npx','nx','ny','nz','ol','th','wp');
load([p 'Detector_Field.mat'],'rp','ws','M','wa','np','xp','yp','ap','ab');


%% Generate intermediate functions
% Start reconstruction closest to center
[~,ii] = sort(xp.^2 + yp.^2 + rp.^2);
I = I(:,:,ii);
xp = xp(ii);
yp = yp(ii);
rp = rp(ii);

% Generate the slit function
slit = @(r0,ws,r) rect((r - r0)./ws);
S = slit(permute(rp,[2 3 1]),ws,y);
S2 = S;

% Generate the probe function
prob = @(x0,y0,wp,x,y,z) sqrt(gaussRMS(x - x0,wp).*gaussRMS(y - tan(2.*th).*z - y0,wp./cos(2.*th)));

% Get the measured amplitudes
amp = sqrt(I);


%% Reconstruction
% Generate initial guess
rho = zeros(ny,nx,nz);

% Set reconstruction parameters
cyc = 200;
beta0 = 0.9;
dec = 0.995;
betamin = 0.8;
beta = beta0;
Q = NaN(np,cyc);
gamma = 0;

% Find the slice to plot examples from
ii = squeeze(sum(sum(abs(I.*x.^2.*y.^2))));
[~,ii] = max(ii);

% Convert to single precision
xp = single(xp);
yp = single(yp);
wp = single(wp);
x = single(x);
y = single(y);
z = single(z);
rho = single(rho);
nx = single(nx);
ny = single(ny);
nz = single(nz);
ap = single(ap);
ab = single(ab);
S = single(S);
S2 = single(S2);
Q = single(Q);
amp = single(amp);
beta = single(beta);
gamma = single(gamma);

% Send to GPU
xp = gpuArray(xp);
yp = gpuArray(yp);
wp = gpuArray(wp);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
rho = gpuArray(rho);
nx = gpuArray(nx);
ny = gpuArray(ny);
nz = gpuArray(nz);
ap = gpuArray(ap);
ab = gpuArray(ab);
S = gpuArray(S);
S2 = gpuArray(S2);
Q = gpuArray(Q);
amp = gpuArray(amp);
beta = gpuArray(beta);
gamma = gpuArray(gamma);

% Pre-fftshift
rho = fftshift(fftshift(rho,1),2);
ap = fftshift(ap);
ab = fftshift(ab);
S = fftshift(S,1);
S2 = fftshift(S2,1);
amp = fftshift(fftshift(amp,1),2);
x = fftshift(x);
y = fftshift(y);

% UI control
figure(101);
ui = uicontrol('style','togglebutton','String','Stop','Position',[10 10 40 20]);

% Reconstruction loop
for k = 1:cyc
    for j = 1:np
        % Generate probe j
        Pj = prob(xp(j),yp(j),wp,x,y,z);
        % Forward
        psij = rho.*Pj; % Apply probe j
        psij = sum(psij,3)./nz; % Projection onto exit beam direction
        psij = fft2(psij)./sqrt(nx.*ny); % Propagate to objective lens
        psij = ab.*ap.*psij; % Apply aperture function
        psij = fft2(psij)./sqrt(nx.*ny); % Propagate to image plane
        psij = S(:,:,j).*psij; % Apply the slit
        psij = fft2(psij)./sqrt(nx.*ny); % Propagate to the detector
        % Error metric
        Q(j,k) = sum(sum((abs(psij) - amp(:,:,j)).^2)); % Calculate RMS error
        % Update Fourier space
        psijt = (abs(psij).*(1 - ap).*gamma + amp(:,:,j)).*exp(1i.*angle(psij)); % Apply the measured amplitudes
        % Gradient
        delj = psij - psijt; % Subtract the initial guess and updated guess
        delj = ifft2(delj).*sqrt(nx.*ny); % Backpropagate to the image plane
        delj = conj(S2(:,:,j)).*delj; % Apply the slit
        delj = ifft2(delj).*sqrt(nx.*ny); % Backpropagate to the objective
        delj = conj(ap.*ab).*delj; % Reverse the lens aberrations
        delj = ifft2(delj).*sqrt(nx.*ny); % Backpropagate to the object
        delj = repmat(delj,1,1,nz); % Expansion into 3D
        delj = conj(Pj).*delj; % Apply the probe
        % Update real space
        rho = rho - beta.*delj; % Update real space
        % Get example data
        if j == ii
            A = psij;
            B = amp(:,:,j);
        end
    end
    if k == 10
        rho = fftshift(fftshift(rho,1),2);
        rho = imgaussfilt(abs(rho),2).*exp(1i.*angle(rho));
        rho = fftshift(fftshift(rho,1),2);
    end
    if beta > betamin
        beta = beta.*dec;
    end
    figure(100);
    semilogy(sum(Q));
    figure(101);
    subplot(2,2,1);
    imagesc(log10(fftshift(abs(A))),[0 3]);
    axis equal tight;
    subplot(2,2,2);
    imagesc(log10(fftshift(B)),[0 3]);
    axis equal tight;
    subplot(2,2,3);
    imagesc(fftshift(abs(rho(:,:,nz/2+1))));
    axis equal tight;
    subplot(2,2,4);
    imagesc(fftshift(angle(rho(:,:,nz/2+1))),[-pi pi]);
    axis equal tight;
    drawnow;
    if ui.Value == 1
        cyc = k;
        break
    end
end

% Post fftshift
rho = fftshift(fftshift(rho,1),2);
ap = fftshift(ap);
ab = fftshift(ab);
S = fftshift(S,1);
S2 = fftshift(S2,1);
amp = fftshift(fftshift(amp,1),2);
x = fftshift(x);
y = fftshift(y);

% Send to CPU
xp = gather(xp);
yp = gather(yp);
wp = gather(wp);
x = gather(x);
y = gather(y);
z = gather(z);
rho = gather(rho);
nx = gather(nx);
ny = gather(ny);
nz = gather(nz);
ap = gather(ap);
ab = gather(ab);
S = gather(S);
S2 = gather(S2);
Q = gather(Q);
amp = gather(amp);
beta = gather(beta);
psij = gather(psij);
psijt = gather(psijt);
delj = gather(delj);
Pj = gather(Pj);
A = gather(A);
B = gather(B);
gamma = gather(gamma);

% Plot the result
Slicer(abs(rho));
Slicer(angle(rho));

% Save the result
save([p 'Reconstruction.mat'],'rho','cyc','dec','beta0','betamin','beta','gamma','Q');


