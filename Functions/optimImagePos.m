function d2c = optimImagePos(E0,x0,D,F,lambda,R0,s0,sigma_v,sigma_p,shwplt)
% Find the optimum d2 value for the true wavefront image plane.

% Set default value of shwplt if not specified
if nargin == 9
    shwplt = 0;
end

% Make an optimization function
fun = @(x) curvature(E0,x0,D,F,lambda,R0,s0,sigma_v,sigma_p,shwplt,x);

% Find the zeropoint
d2c = fzero(fun,0);

end

function C = curvature(E0,x0,D,F,lambda,R0,s0,sigma_v,sigma_p,shwplt,d2c)

% Update the distance
D(end) = D(end) + d2c;

% Calculate the propagation parameters
[a,~,Rp,sm,sp,~,~] = FrFT_parameters(D,F,lambda,R0,s0);

% Vignetting
E = E0.*sqrt(exp(-x0.^2./(2*sigma_v.^2)));

% Propagation to CRL exit plane
[E,x] = propFrFT1(E,x0,R0,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,0);

% Apply effective pupil
E = E.*sqrt(exp(-x.^2./(2*sigma_p.^2)));

% Propagate to image plane
[E,x] = propFrFT1(E,x,Inf,Rp(end),sm(end),sp(end),a(end),lambda,0);

% Fit the wavefront
R = wavefrontCurvature(x,E,lambda,shwplt);

C = 1./R;

end
