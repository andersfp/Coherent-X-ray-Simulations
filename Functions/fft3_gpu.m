function F = fft3_gpu(f,el_opt)
% 3D FFT of 3D array performed on GPU, maximizing computation speed. The
% optional parameter 'el_opt' specifies the optimum number of elements that
% maximizes computation speed on a specific GPU.

% Set the default value of el_opt (GTX Titan X (Pascal))
if nargin == 1
    %el_opt = 4.98e7;
    el_opt = 1.95e8;
end

% Perform the FFT along the first dimension
F = fft1_gpu(f,el_opt);

% Permute the array
F = permute(F,[2 3 1]);

% Perform the FFT along the second dimension
F = fft1_gpu(F,el_opt);

% Permute the array
F = permute(F,[2 3 1]);

% Perform the FFT along the third dimension
F = fft1_gpu(F,el_opt);

% Permute the array
F = permute(F,[2 3 1]);

end

function F = fft1_gpu(f,el_opt)
% Perform the 1D FFT on a 3D input, splitting along the 3rd dimension,
% optimizing performance.

% Get the size of the input
s1 = size(f,1);
s2 = size(f,2);
s3 = size(f,3);

% Get indices for the 
m = floor(el_opt/(s1*s2));
k = ceil(s3/m);
i1 = (0:(k - 1))*m + 1;
i2 = (1:k)*m;
i2(end) = s3;

% Perform the stepwise FFT
F = zeros(s1,s2,s3);
for i = 1:k
    % Send data to GPU
    fg = gpuArray(f(:,:,i1(i):i2(i)));
    
    % Perform the FFT
    fg = fft(fg);
    
    % Return the result to the CPU
    F(:,:,i1(i):i2(i)) = gather(fg);
end

end

