function B = bin3(A,b)
% Bin the N-D image A by b = [b1 b2 ...] (b must be positive integers). The
% size of A must be an integer time b on all sides.

% Check the number of pixels in A
if any(mod(size(A),b))
    error('The number of pixels in all dimensions of A must be an integer times b = [b1 b2 ...].');
end

% Separate the binning parameters
b1 = b(1);
b2 = b(2);
b3 = b(3);

% Bin the image
B = zeros(size(A)./b);
for i = 1:b1
    for j = 1:b2
        for k = 1:b3
            B = B + A(i:b1:end,j:b2:end,k:b3:end);
        end
    end
end

% Normalize the output
B = B./prod(b);

