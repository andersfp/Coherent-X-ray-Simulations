function u2 = propAS1(u1,lambda,d1,z)
% Propagate - using the angular spectrum method. Special case using m = 1.

% Get the size of the input
N = size(u1,1); % assume square grid

% Decide which implementation to use
if N^2 <= 79210000
    u2 = propAS1gpus(u1,lambda,d1,z);
elseif N^2 <= 110250000
    u2 = propAS1gpuss(u1,lambda,d1,z);
elseif N > 40000
    u2 = propAS1gpu(u1,lambda,d1,z);
else
    u2 = propAS1cpu(u1,lambda,d1,z);
end

