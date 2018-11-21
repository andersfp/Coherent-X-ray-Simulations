function [field,object,chi_square_real_last,chi_square_reci_last]=fienupHIOso(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci);

% Save the last object
object_last = object;

% Apply reciprocal space contraints
field_new = data_input.*exp(1i.*angle(field)).*mask + (1 - mask).*field;

% Transform to real space
new_object_0 = fftshift(ifftn(fftshift(field_new))); % fftn

% Find the support condition
condition = (1 - (support == 0)); %.*(abs(objet)>=0.5*support) & (abs(objet)<=1.5*support);

% Apply the real space constraint
new_object_1 = new_object_0.*condition + (1 - condition).*(object_last - beta.*new_object_0);

% Transform to reciprocal space
field = fftshift(fftn(fftshift(new_object_1))); % ifftn

% Calculate error metric
[chi_square_real_last] = error_metric_real(new_object_0,new_object_1,chi_square_real);
[chi_square_reci_last] = error_metric_reci(field_new,field,chi_square_reci);

% Update the new object
object = new_object_1;

