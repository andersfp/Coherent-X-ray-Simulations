function [bx,by,bz,a] = optshft_fit(bx,by,bz,a,rho,data_input,mask,flag)
% Set the flag is not specified
if nargin == 7
    flag = 0;
end

% Get the measured intensities
I = data_input.^2;

% Get number of datasets
n = size(I,4);

% Generate a fake Gaussian mask
G = ones(size(I,1),size(I,2),'like',I);

% Optimization options
opt = optimset('Display','off','TolFun',1e10,'TolX',1e-3);

% Run the fit in different modes
switch flag
    case 0
        % Initial guess
        x0 = cat(1,bx.',by.',bz.',a.');
        
        % Optimize each dataset
        xmin = zeros(4,n);
        fval = zeros(1,n);
        for i = 1:n
            fun = @(x) optshft(rho,I(:,:,:,i),mask,G,x(1),x(2),x(3),x(4));
            [xmin(:,i),fval(i)] = fminsearch(fun,x0(:,i),opt);
        end
        
        % Set the outputs
        bx = xmin(1,:).';
        by = xmin(2,:).';
        bz = xmin(3,:).';
        a = xmin(4,:).';
        
    case 1
        % Initial guess
        x0 = cat(1,bx.',by.',a.');
        
        % Optimize each dataset
        xmin = zeros(3,n);
        fval = zeros(1,n);
        for i = 1:n
            fun = @(x) optshft(rho,I(:,:,:,i),mask,G,x(1),x(2),bz(i),x(3));
            [xmin(:,i),fval(i)] = fminsearch(fun,x0(:,i),opt);
        end
        
        % Set the outputs
        bx = xmin(1,:).';
        by = xmin(2,:).';
        a = xmin(3,:).';
end

