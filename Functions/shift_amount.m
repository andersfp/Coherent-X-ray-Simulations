function shft = shift_amount(d1,d3,th)
% Calculates the number of pixels to shear the 3D volume to make it
% orthogonal. The input is the pixel size in the vertical direction (1st
% dimension), the pixel size in the depth direction (3rd dimension), and
% the Bragg angle th.
%
% Example of usage:
% shft = shift_amount(d1,d3,th);
%

% Calculate amount to shift in the vertical direction
shft = d3.*sin(th)./d1;

