function B = bin2(A,b)
% Bin the 2D image A by bxb (b must be positive integer). The size of A
% must be an integer time b on both sides.

% Check the number of pixels in A
if any(mod(size(A),b))
    error('The number of pixels in both dimensions of A must be an integer times b.');
end

% Bin the image
B = zeros(size(A)/b);
for i = 1:b
    for j = 1:b
        B = B + A(i:b:end,j:b:end);
    end
end

% Normalize the output
B = B./b.^2;

