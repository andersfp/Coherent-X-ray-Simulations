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
theta_M = asind(q111*lambda/(4*pi));
theta_S = asind(q111*lambda/(4*pi));
theta_A = asind(q220*lambda/(4*pi));

% Calculate streak angles
zeta = atand((tand(theta_A) - tand(theta_M))/((tand(theta_M) + tand(theta_A) - 2*tand(theta_S))*tand(theta_S)));
zeta = zeta/2;
gam = -atand((tand(theta_A) - 2*tand(theta_M))*tand(theta_S)/(2*tand(theta_M) + tand(theta_A) - 2*tand(theta_S)));

% Plot the streaks
lw = 3;
fnt = 14;
ofs = -abs(theta_A);
figure;
plot([-sind(-zeta + ofs) sind(-zeta + ofs)],[cosd(-zeta + ofs) -cosd(-zeta + ofs)],'LineWidth',lw);
hold on;
plot([-sind(gam + ofs) sind(gam + ofs)],[cosd(gam + ofs) -cosd(gam + ofs)],'LineWidth',lw);
plot([-sind(90 + ofs) sind(90 + ofs)],[cosd(90 + ofs) -cosd(90 + ofs)],'LineWidth',lw);
plot([-sind(-theta_S + ofs) sind(-theta_S + ofs)],[cosd(-theta_S + ofs) -cosd(-theta_S + ofs)],'LineWidth',lw);
legend('\lambda','mono 1','mono 2','sample','Location','NorthWest');
set(gca,'XLim',[-1e-4 1.2e-4],'YLim',[-8e-4 7e-4],'FontSize',fnt);
xlabel('\DeltaQ_{rock}/Q');
ylabel('\DeltaQ_{2\theta}/Q');


