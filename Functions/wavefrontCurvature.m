function R = wavefrontCurvature(x,E,lambda,shwplt)
% Given the 1D electric field, fit the radius of curvature of the phase
% within the central part (part without phase wrapping).

% Set default value of shwplt if not specified
if nargin == 3
    shwplt = 0;
end

% Get the number of points
m = length(E);

% Get the amplitude
a = abs(E);
a = a./max(a);

% Get the phase
p = angle(E.*exp(-1i.*angle(E(m/2 + 1))));

% Find the exclusion
s = [true;abs(diff(p)) > pi] | a < 0.01;

% Find the phase limit
xl = min(abs(x(s)));

% Make the fitting model
fun = @(R,x) pi.*x.^2./(lambda.*R);

% Make a starting guess for the curvature
R0 = pi.*xl.^2./(lambda.*p(x == xl));

% Fit the phase
f = fit(x,p,fun,'Exclude',abs(x) > xl,'StartPoint',R0);

% Plot the fit (debugging)
if shwplt
    figure;
    plot(f,x,p);
end

% Get the radius
R = f.R;

