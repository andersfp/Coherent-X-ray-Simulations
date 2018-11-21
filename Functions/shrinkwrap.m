function [support,sig] = shrinkwrap(object,x,y,z,sig,sigmin,tol,support_old)
% Generate a new support based on the current state of the object.

% Check if support_old has been specified
if nargin == 7
    support_old = object;
end

% Generate a Gaussian profile
G = exp(-(x.^2 + y.^2 + z.^2)./(2*sig.^2));

% Convolute the object and Gaussian
support = fconvn(abs(object),G);

% Cut off the support at given tolerance level
switch class(support_old)
    case 'double'
        support = double(support > tol*max(support(:)));
    case 'single'
        support = single(support > tol*max(support(:)));
    case 'logical'
        support = support > tol*max(support(:));
    case 'gpuArray'
        switch classUnderlying(support_old)
            case 'double'
                support = double(support > tol*max(support(:)));
            case 'single'
                support = single(support > tol*max(support(:)));
            case 'logical'
                support = support > tol*max(support(:));
        end
end

% Fill voids
support = imfill(support,'holes');

% Update the Gaussian width
sig = 0.99*sig;
if sig < sigmin
    sig = sigmin;
end

