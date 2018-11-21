function v = PSF_Fitting(r,f,w,ws,sig0)
% Fit the PSF of a reconstruction in the Fourier domain by use of
% MTF_Fitting.m.

% Set fitting options
opt = optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');

% Perform the fitting
v = fminsearch(@(v) 10*sum(((f - MTF_Fitting(r,[1000;ws;v])).*w).^2),sig0,opt);

