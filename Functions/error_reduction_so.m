function [field_last,new_object_1,chi_square_real_last,chi_square_reci_last]=error_reduction_so(data_input,support,field,objet_ini,mask,chi_square_real,chi_square_reci);
% Use the error reduction algorithm to update the real space and reciprocal
% space.

% Generate estimate of reciprocal space
field_new = data_input.*exp(1i.*angle(field)).*mask + field.*(1 - mask);

% Transform to real space
new_object_0 = fftshift(ifftn(fftshift(field_new))); % fftn

% Apply the real space support
new_object_1 = new_object_0.*support;

% Get a new estimate of the reciprocal space
field_last = fftshift(fftn(fftshift(new_object_1))); % ifftn

% Calculate error metric
[chi_square_real_last] = error_metric_real(new_object_0,new_object_1,chi_square_real);
[chi_square_reci_last] = error_metric_reci(field_new,field,chi_square_reci);

