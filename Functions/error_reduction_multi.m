function [field,object,chi_square_real,chi_square_reci] = error_reduction_multi(data_input,support,field,object,mask,chi_square_real,chi_square_reci)
% Apply the error reduction for a number of separate data steps

% Get the number of data sets
k = size(data_input,4);

% Run the ER algorithm k times
for i = 1:k
    [field,object,chi_square_real,chi_square_reci] = error_reduction(data_input(:,:,:,i),support,field,object,mask(:,:,:,i),chi_square_real,chi_square_reci);
end

