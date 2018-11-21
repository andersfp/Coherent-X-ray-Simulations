function A = shear(A,shft)
% Shear the 3D array A in the 1st direction (vertical) by wrapping each
% slice consequtively by the amount shft.
%
% Example of usage:
% A = shear(A,shft);

% Get the number of pixels in the 3rd dimension
n = size(A,3);

% Shift the object
for i = 1:n
    A(:,:,i) = circshift(A(:,:,i),-round((i - n/2)*shft),1);
end


