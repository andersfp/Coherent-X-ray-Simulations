function [field,object,chi_square_real,chi_square_reci] = fienupHIO_multi(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci)
% Apply the Fienup HIO for a number of separate data steps

% Get the number of data sets
k = size(data_input,4);

% Run the Fienup algorithm k times
for i = 1:k
    [field,object,chi_square_real,chi_square_reci] = fienupHIO(beta,data_input(:,:,:,i),support,field,object,mask(:,:,:,i),chi_square_real,chi_square_reci);
end

