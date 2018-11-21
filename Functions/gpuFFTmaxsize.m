function N = gpuFFTmaxsize()
% Calculate the power of 2 of the data size the gives the maximum FFT
% performance of the GPU.
%
% Example of usage:
% N = gpuFFTmaxsize();
% m = 2^N;
%

% Set the initial power of 2 and speed
n = 0;
s0 = 0;
s = 1;

% Set number of repetitions at each step
k = 10;
t = zeros(k,1);

% Compute the FFT speed until it slows down
while s > 0.1*s0
    n = n + 1;
    m = 2^n;
    A = rand(m,1);
    B = gpuArray(A);
    C = fft(B);
    for j = 1:k
        tic;
        C = fft(B);
        t(j) = toc;
    end
    s0 = s;
    s = m.*n./mean(t);
end

% Set the output
N = n - 1;

