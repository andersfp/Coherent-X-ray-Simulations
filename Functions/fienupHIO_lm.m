function [field,object,chi_square_real,chi_square_reci] = fienupHIO_lm(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci)

% Apply reciprocal space contraints
field = data_input.*exp(1i.*angle(field)).*mask + (1 - mask).*field;

% Transform to real space
field = fftshift(ifftn(fftshift(field)));

% Apply the real space constraint
if islogical(support)
    object = field.*support + (~support).*(object - beta.*field);
else
    object = field.*support + (1 - support).*(object - beta.*field);
end

% Transform to reciprocal space
field = fftshift(fftn(fftshift(object)));

% Calculate error metric
chi_square_real = [chi_square_real NaN];
chi_square_reci = error_metric_data(data_input,field,chi_square_reci);

