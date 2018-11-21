function r = Single_Prop_Fit1(E0,x0,sm,sp,a,lambda,EN,sigpsf,sigvig)

% Apply vignetting
attvig = sqrt(exp(-x0.^2./(2.*sigvig.^2)));
E = E0.*attvig;

% Propagate to CRL exit
[E,x] = propFrFT1(E,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,0);

% Apply aperture
attpsf = sqrt(exp(-x.^2./(2.*sigpsf.^2)));
E = E.*attpsf;

% Calculate residual
r = sum(abs(abs(EN).^2 - abs(E).^2).^2);

