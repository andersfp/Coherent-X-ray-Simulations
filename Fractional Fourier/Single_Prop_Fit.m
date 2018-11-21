function r = Single_Prop_Fit(E0,x0,sm,sp,a,lambda,m,LM,sigpsf,sigvig)

% Propagate to CRL entrance
[E,~,x] = propFrFT2(E0,x0,x0,Inf,Inf,sm(1),sp(2),a(2),lambda,0,'gpu');

% Apply vignetting
attvig = sqrt(exp(-(x.^2 + x.'.^2)./(2*(sigvig)^2)));
E = E.*attvig;

% Propagate to CRL exit
[E,~,x] = propFrFT2(E,x,x,Inf,Inf,sm(2),sp(end-1),sum(a(2:end-1)),lambda,0,'gpu');

% Apply aperture
attpsf = sqrt(exp(-(x.^2 + x.'.^2)./(2*(sigpsf)^2)));
E = E.*attpsf;

% Propagate to detector
[E,~,~] = propFrFT2(E,x,x,Inf,Inf,sm(end),sp(end),a(end),lambda,0,'gpu');

% Calculate residual
L = abs(E(:,m/2+1)).^2;
r = sum(abs(LM - L).^2);

