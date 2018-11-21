function [field,object,chi_square_real,chi_square_reci] = error_reduction_2(epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci)
% Use the error reduction algorithm to update the real space and reciprocal
% space.

% Extract the tolerance data
h = epsilon.h;
data_mean = epsilon.data_mean;
lim = epsilon.lim;

% Check the amplitude limits
ii = abs(field(h)) > lim;
jj = h(ii);
field(jj) = data_mean(ii).*exp(1i.*angle(field(jj)));

% Generate estimate of reciprocal space
field_new = data_input.*exp(1i.*angle(field)).*mask + field.*(1 - mask);

% Transform to real space
new_object_0 = fftshift(ifftn(fftshift(field_new)));

% Apply the real space support
object = new_object_0.*support;

% Get a new estimate of the reciprocal space
field = fftshift(fftn(fftshift(object)));

% Calculate error metric
chi_square_real = error_metric_real(new_object_0,object,chi_square_real);
chi_square_reci = error_metric_reci(field_new,field,chi_square_reci);

