% Initialization
clear;
close all;
clc;


%% Sample
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;


%% CRL calculations
% Get material properties
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6;
T = 1.6e-3;
Tweb = 2e-6;
N = 64;

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set object and image positions
d1 = 0.32202;
d2 = fN.*(d1 + f.*phi.*tan(N.*phi))./(d1 - fN);
M = d2.*sin(N.*phi)./(f.*phi) - cos(N.*phi);

% Calculate apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate vignetting
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);


%% Parameters
% Calculate sigma_CRL
sigma_CRL = sigma_D/fN;

% Calculate the resolution in BFP
Delta_y_BFP = 0.316*lambda/sigma_CRL;

% Calculate the angular sensitivity in BFP
Delta_alpha_BFP = 0.316*lambda*cos(N*phi)/sigma_D;

% Calculate the BFP FOV
FOV_BFP = sigma_a*fN/cos(N*phi);

% Set the beam line parameters
sigma_e = 6e-5;
theta_B = 21.88/2*pi/180;
Delta_zeta_v = 0.029e-3;
Delta_zeta_h = 0.043e-3;
Q0 = 4*pi/lambda*sin(theta_B);

% Calculate reciprocal space resolutions
Delta_Q_rock = Q0/2*sqrt(Delta_zeta_v^2 + sigma_a^2);
Delta_Q_roll = Q0/(2*sin(theta_B))*sqrt(Delta_zeta_h^2 + sigma_a^2);
Delta_Q_par = Q0/2*sqrt((2*sigma_e)^2 + cot(theta_B)^2*(Delta_zeta_v^2 + sigma_a^2));
Delta_Q = [Delta_Q_rock Delta_Q_roll Delta_Q_par];

% Calculate relative reciprocal space resolutions
Delta_Q_rock_rel = Delta_Q_rock/Q0;
Delta_Q_roll_rel = Delta_Q_roll/Q0;
Delta_Q_par_rel = Delta_Q_par/Q0;
Delta_Q_rel = [Delta_Q_rock_rel Delta_Q_roll_rel Delta_Q_par_rel];

% Calculate the FOV shift
r_shift = fN*gamma/cos(N*phi);


%% Triple crystal diffractometer
% Set silicon lattice parameter
a = 543.102e-12;

% Calculate momentum transfers
b = 2*pi/a;
q111 = sqrt(3)*b;
q220 = sqrt(8)*b;
q400 = 4*b;

% Calculate Bragg angles
E = 100e3;
lambda = 1e-10*12398.42./E;
theta_M = asind(q220*lambda/(4*pi));
theta_S = asind(q111*lambda/(4*pi));
theta_A = asind(q220*lambda/(4*pi));

% Calculate streak angles
zeta = atand((tand(theta_A) - tand(theta_M))/((tand(theta_M) + tand(theta_A) - 2*tand(theta_S))*tand(theta_S)));
zeta = zeta/2;
gam = -atand((tand(theta_A) - 2*tand(theta_M))*tand(theta_S)/(2*tand(theta_M) + tand(theta_A) - 2*tand(theta_S)));

% Plot the streaks
lw = 3;
fnt = 14;
figure;
plot([0 0],[0 1],'LineWidth',lw);
hold on;
plot([-sind(-zeta) sind(-zeta)],[cosd(-zeta) -cosd(-zeta)],'LineWidth',lw);
plot([-sind(gam) sind(gam)],[cosd(gam) -cosd(gam)],'LineWidth',lw);
plot([-1 1],[0 0],'LineWidth',lw);
plot([-sind(-theta_S) sind(-theta_S)],[cosd(-theta_S) -cosd(-theta_S)],'LineWidth',lw);
legend('\bf G','\lambda','M','S','A');
axis equal tight;
set(gca,'XLim',[-1 1],'YLim',[-1 1],'FontSize',fnt);
xlabel('q_y');
ylabel('q_x');


