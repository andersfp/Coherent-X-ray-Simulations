function fom = optshft(rho,I,mask,G,bx,by,bz,a)
% Get the data size
nx = size(I,2);
ny = size(I,1);
n_omega = size(I,3);

% Generate pixel coordinates
x = (-nx/2):(nx/2 - 1);
y = ((-ny/2):(ny/2 - 1)).';
z = permute((-n_omega/2):(n_omega/2 - 1),[3 1 2]);

% Calculate the phase shifts
px = -1i.*2.*pi.*bx./nx.*x;
py = -1i.*2.*pi.*by./ny.*y;
pz = -1i.*2.*pi.*bz./n_omega.*z;

% Send phase shifts to GPU
px = gpuArray(px);
py = gpuArray(py);
pz = gpuArray(pz);

% Apply the phase factors
tmp = rho.*exp(px).*exp(py).*exp(pz);

% Fourier transform to the detector space
tmp = fftshift(fftn(fftshift(tmp)));

% Calculate intensity
tmp = abs(tmp).^2;

% Apply detector mask and pupil function
tmp = tmp.*mask.*G;

% Calculate RMS difference between simulated and measured intensity
fom = gather(sum(sum(sum(sum(abs(I - a.*tmp).^2)))))./a.^2;

