function [field,object,chi_square_real,chi_square_reci] = fienupHIO_2(epsilon,beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci)

% Extract the tolerance data
h = epsilon.h;
data_mean = epsilon.data_mean;
lim = epsilon.lim;

% Check the amplitude limits
ii = abs(field(h)) > lim;
jj = h(ii);
field(jj) = data_mean(ii).*exp(1i.*angle(field(jj)));

% Apply reciprocal space contraints
field_new = data_input.*exp(1i.*angle(field)).*mask + (1 - mask).*field;

% Transform to real space
new_object_0 = fftshift(ifftn(fftshift(field_new)));

% Apply the real space constraint
if islogical(support)
    object = new_object_0.*support + (~support).*(object - beta.*new_object_0);
else
    object = new_object_0.*support + (1 - support).*(object - beta.*new_object_0);
end

% Transform to reciprocal space
field = fftshift(fftn(fftshift(object)));

% Calculate error metric
chi_square_real = error_metric_real(new_object_0,object,chi_square_real);
chi_square_reci = error_metric_reci(field_new,field,chi_square_reci);

