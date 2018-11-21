function A = shear_interp(A,shft)
% The the 3D object A in the vertical-depth plane by an amount of shft per
% pixel in the 3rd dimension.
%
% Example of usage:
% A = shear_interp(A,shft);
%

% Get the size of the input
s1 = size(A,1);
s2 = size(A,2);
s3 = size(A,3);

% Generate axes
r1 = (-s1/2):(s1/2 - 1);
r2 = (-s2/2):(s2/2 - 1);
r3 = (-s3/2):(s3/2 - 1);

% Generate a meshgrid
[r1q,r2q,r3q] = ndgrid(r1,r2,r3);

% Correct the meshgrid
r1q = r1q + shft.*r3q;

% Perform the interpolation
A = interpn(r1,r2,r3,A,r1q,r2q,r3q,'linear',0);

