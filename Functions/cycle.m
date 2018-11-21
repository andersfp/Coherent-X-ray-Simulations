function [field,object,chi_square_real,chi_square_reci] = cycle(mode,iter,beta,epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci)
% Run the given number of cycles using the selected algorithm

% Select the algorithm
switch mode
    case 1
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = fienupHIO(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 2
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = error_reduction(data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 3
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = fienupHIO_2(epsilon,beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 4
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = error_reduction_2(epsilon,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 5
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = fienupHIO_lm(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 6
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = error_reduction_lm(data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 8
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = fienupHIOso(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 18
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = error_reduction_so(data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 11
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = fienupHIO_multi(beta,data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
    case 12
        for l=1:iter
            [field,object,chi_square_real,chi_square_reci] = error_reduction_multi(data_input,support,field,object,mask,chi_square_real,chi_square_reci);
        end
end

