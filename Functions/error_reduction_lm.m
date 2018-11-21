function [field,object,chi_square_real,chi_square_reci] = error_reduction_lm(data_input,support,field,object,mask,chi_square_real,chi_square_reci)
% Use the error reduction algorithm to update the real space and reciprocal
% space.

% Generate estimate of reciprocal space
field = data_input.*exp(1i.*angle(field)).*mask + field.*(1 - mask);

% Transform to real space
object = fftshift(ifftn(fftshift(field)));

% Apply the real space support
object = object.*support;

% Get a new estimate of the reciprocal space
field = fftshift(fftn(fftshift(object)));

% Calculate error metric
chi_square_real = [chi_square_real NaN];
chi_square_reci = error_metric_data(data_input,field,chi_square_reci);

