function u2 = propAS1gpu(u1,lambda,d1,z)
% Propagate - using the angular spectrum method. Special case using m = 1.

% Get the size of the input
N = size(u1,1); % assume square grid

% Get the wavenumber
q = 2*pi/lambda; % optical wavevector

% Generate spatial frequency (of source plane) coordinates
df1 = 1/(N*d1);
[fX,fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
fsq = fX.^2 + fY.^2;
clear fX fY;

% Make shift array to shift data
shft = [(N/2+1):N 1:(N/2)];

% Calculate quadratic phase factors
Q2 = exp(-1i*pi^2*2*z/q*fsq);
clear fsq;
Q2 = Q2(shft,shft);

% Make the transform
u2(shft,shft) = fft2gpu(Q2(shft,shft).*fft2gpu(u1(shft,shft),1),-1);

end

function res = fft2gpu(fc,flag)
% Make a 2D Fourier transform on the GPU on a large array

% Do the 1D transform
res = fft1gpu(fc,flag);

% Do the second 1D transform
res = res';
res = fft1gpu(res,flag);
res = res';

end

function res = fft1gpu(fc,flag)
% Make a 1D Fourier transform on the GPU on a large array

% Get the size of the square input
N = size(fc,1);

% Split the array if too large to fit in GPU memory
elmax = 1e8;
k = ceil(numel(fc)/elmax);
m = ceil(N/k);
madd = k*m - N;
fc = cat(2,fc,zeros(N,madd));
fc = reshape(fc,N,m,k);

% Pre-allocate memory for the result
res = zeros(N,m,k);

% Do the Fourier transform
if flag == 1
    for i = 1:k
        fg = gpuArray(fc(:,:,i));
        fg = fft(fg);
        res(:,:,i) = gather(fg);
    end
elseif flag == -1
    for i = 1:k
        fg = gpuArray(fc(:,:,i));
        fg = ifft(fg);
        res(:,:,i) = gather(fg);
    end
else
    error('The forward/inverse FFT flag must be either 1 or -1.');
end

% Unfold the array and remove zeros
res = reshape(res,N,m*k,1);
res(:,end-madd+1:end) = [];

end


