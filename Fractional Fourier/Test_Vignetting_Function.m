% Initialization
clear;
close all;
clc;


%% Set parameters
% Set CRL physical dimensions
R_CRL = 50e-6;
T = 1.6e-3;
Tweb = 2e-6;
%N = 69;

% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42/E;

% Material parameters
[delta,mu] = Be_Prop(E);

% Initialize GPU
t = initGPU();

% Set object and image positions (see CRL_Focal_Length_Sim)
%D1 = linspace(7,0,30);
D1 = [linspace(8.5,1,31) linspace(0.95,0.1,18) linspace(0.09,0,10)];
%NN = 10:10:100;
NN = [2:10 15 20:10:170];
n1 = length(D1);
n2 = length(NN);
SA = zeros(n1,n2);
SV = SA;
FF = SA;
YN = SA;
YM = SA;
RE = SA;
h = waitbar(0,'Progress');
for j1 = 1:n1
    for j2 = 1:n2
        % Progress bar
        waitbar((j1 - 1)/n1 + j2/(n1*n2),h);
        
        % Variable parameters
        d1 = D1(j1);
        N = NN(j2);
        
        % Calculate CRL optical parameters
        f = R_CRL/(2*delta);
        phi = sqrt(T/f);
        fN = f*phi*cot(N*phi);
        
        % CRL parameters
        S1 = N.*sinc(2.*N.*phi./pi);
        S2 = sin(N.*phi).^2./phi;
        d = d1./(f.*phi);
        sigA = sqrt(R_CRL./mu)./(f.*phi)./sqrt(N.*(1 + d.^2) - (1 - d.^2).*S1 + 2.*d.*S2);
        sigV = sqrt(R_CRL./mu).*sqrt((N.*(1 + d.^2) - (1 - d.^2).*S1 + 2.*d.*S2)./(N.^2 - S1.^2 - S2.^2));
        
        % Theoretical PSF width
        y0 = sigA*d1;
        a0 = sigA;
        yn = y0*(cos((1:N)*phi) + phi/2*sin((1:N)*phi)) + a0*(f*phi*sin((1:N)*phi) - T/2*cos((1:N)*phi));
        yN = yn(end);
        
        % Calculate the vignetting
        tic;
        [ff,r] = Vignetting(R_CRL,N,mu,d1,T,f,lambda,yN);
        toc;
        
        % Save parameters
        SA(j1,j2) = sigA;
        SV(j1,j2) = sigV;
        FF(j1,j2) = ff;
        YN(j1,j2) = yN;
        YM(j1,j2) = max(yn);
        RE(j1,j2) = r;
    end
end
close(h);

% Save the simulations
save('Vignetting_Tests.mat','-v7.3');

% Parameters
C = FF./YN.*sqrt(SA);
%C = (C - 2*phi).*NN;
c = mean(C);
cs = std(C);
w = 1./c;


