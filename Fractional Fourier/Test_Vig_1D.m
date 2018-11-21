% Initialization
clear;
close all;
clc;


%% Set parameters
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42/E;

% Set CRL parameters
R_CRL = 50e-6;
delta = 1.17843774e-6;
f = R_CRL/(2*delta);
T = 1.6e-3;
N = 170;
mu = 1/24590.6e-6;
Tweb = 2e-6;

% Calculate focal length
phi = sqrt(T/f);
fN = f*phi*cot(N*phi);

% Set object and image positions (see CRL_Focal_Length_Sim)
%D1 = linspace(7,0,30);
D1 = [linspace(8.5,1,31) linspace(0.95,0.1,18) linspace(0.09,0,10)];
n = length(D1);
SA = zeros(n,1);
SV = SA;
FF = SA;
YN = SA;
YM = SA;
h = waitbar(0,'Progress');
for j = 1:n
waitbar(j/n,h);
d1 = D1(j);
d2 = 1;

% CRL parameters
S1 = N.*sinc(2.*N.*phi./pi);
S2 = sin(N.*phi).^2./phi;
d = d1./(f.*phi);

% Calculate sigma_a
%sigA = sqrt(R_CRL / (mu * N * (d1^2 + (f * phi)^2)))/sqrt(1 + 1 / N - 1 / (N * phi) * sin((N + 1) * phi) * cos((N - 1) * phi + 2 * atan(d1 / (f * phi))));
sigA = sqrt(R_CRL./mu)./(f.*phi)./sqrt(N.*(1 + d.^2) - (1 - d.^2).*S1 + 2.*d.*S2);

% Theoretical PSF width
y0 = sigA*d1;
a0 = sigA;
yn = y0*(cos((1:N)*phi) + phi/2*sin((1:N)*phi)) + a0*(f*phi*sin((1:N)*phi) - T/2*cos((1:N)*phi));
yN = yn(end);

% Calculate sigma_v
sigV = sqrt(R_CRL./mu).*sqrt((N.*(1 + d.^2) - (1 - d.^2).*S1 + 2.*d.*S2)./(N.^2 - S1.^2 - S2.^2));


%% Make the object
% Set width parameter based on CRL parameters
w0 = sqrt(R_CRL/(2*N*mu));

% Set up object parameters
m = 100;
dx = 15*w0;

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
w = 10*w0;
E0 = rect(x0/w);


%% Calculate the wave propagation
% Make the propagation distance and focal length arrays
D = [d1+T/2;T*ones(N-1,1);T/2+d2];
F = f*ones(N,1);

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);
[a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda,R0,s0);

% Initialize GPU
t = initGPU();

% Propagate in N steps with attenuation
att = @(x) sqrt(exp(-mu*x.^2./R_CRL));
EN = E0;
xN = x0;
for i = 1:N
    [EN,xN] = propFrFFT1(EN,xN,Inf,Inf,sm(i),sp(i),a(i),lambda,'gpu');
    EN = EN.*att(xN);
end

% Fit the CRL plane and image plane
fun = @(sig) Single_Prop_Fit1(E0,x0,sm,sp,a,lambda,EN,yN,sig);
opt = optimset('TolX',1e-8);
[ff,fval,exitflag,output] = fminsearch(fun,sqrt(sigV)/200,opt);

% Propagate in a single step
[E1,x1] = propFrFFT1(E0.*sqrt(exp(-x0.^2./(2.*ff.^2))),x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,'gpu');
E1 = E1.*sqrt(exp(-x1.^2./(2.*yN.^2)));

% Plot the fitting result
figure;
plot(1e6*xN,abs(EN).^2,'-',1e6*x1,abs(E1).^2,'--');
xlabel('Position [\mum]');
ylabel('Intensity [a.u.]');
legend('Lens-by-lens','Fitting');

% Save parameters
SA(j) = sigA;
SV(j) = sigV;
FF(j) = ff;
YN(j) = yN;
YM(j) = max(yn);
end
close(h);


% %% Make plots
% % Object
% figure;
% plot(1e6*x0,abs(E0).^2);
% title('Object intensity');
% xlabel('Position [\mum]');
% ylabel('Intensity [a.u.]');
% 
% % Exit field intensity (lens-by-lens)
% figure;
% plot(1e6*xN,abs(EN).^2);
% title('Exit intensity (lens by lens attenuation)');
% xlabel('Position [\mum]');
% ylabel('Intensity [a.u.]');





