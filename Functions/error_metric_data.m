function chi_square = error_metric_data(data_input,field,chi_square)

chi_square_new = sum(sum(sum((data_input - abs(field)).^2)))./numel(data_input);

chi_square = [chi_square chi_square_new];

