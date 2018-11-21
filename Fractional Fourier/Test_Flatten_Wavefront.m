% Initialization (optional)
clear;
close all;
clc;


%% Decide how to calculate the FrFT
% Define method
%method = 'gpu'; % Fastest implementation, but requires an Nvidia GPU.
method = 'vec'; % Fastest implementation on CPU (for not too large images), based on vectorized array operations.
% method = 'par'; % Parallel for-loop implementation on CPU, can be faster than 'vec' for large images, and uses less memory than 'vec' requires parallel toolbox.
% method = 'for'; % Non-parallel for-loop on CPU, the slowest implementation, but has the lowest memory consumption and does not require parallel toolbox.

% Save the figure
sav = 0;


%% Define the object
% Number of pixels on each side (must be even)
m = 4096;

% Field-of-view
dx = 10e-6; % 1e-3, 1e-6

% Make coordinates
x0 = ((-m/2):(m/2 - 1)).'/m*dx;

% Make object
w = 10e-6; % 600e-6, 1e-9
%E0 = recta(x0./w);
E0 = zeros(m,1);
E0(m/2 + 1) = 1;


%% Set the simulation details
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Get material properties (delta and mu) at this energy
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6; % Radius of curvature at apex
T = 1.6e-3; % Distance between each lens element
Tweb = 2e-6; % Minimum Be thickness at the optical axis
N = 69; % Number of lenses

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set object and image plane positions (see CRL paper)
M = 10;
d2 = f.*phi.*(M + cos(N.*phi))./sin(N.*phi);
d1 = fN.*(d2 + f.*phi.*tan(N.*phi))./(d2 - fN);

% Calculate analytical apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate the effective vignetting width
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);

% Make an array of the focal lengths
F = f*ones(N,1);

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);

% Make an array of distances in between each lens (plus object-lens and lens-detector)
D = [d1 + T/2;T*ones(N-1,1);T/2 + d2];

% Find the optimal position using optimization
tic;
d2c = optimImagePos(E0,x0,D,F,lambda,R0,s0,sigma_v,sigma_p);
toc;

% Set distance variations
dd = linspace(-2e-2,2e-2,201).';
%dd = linspace(-4.65e-4,-4.6e-4,201).';

% Vary the d2 distance
n = length(dd);
R2 = zeros(n,1);
w0 = R2;
E = zeros(m,n);
X = E;
for i = 1:n
% Make an array of distances in between each lens (plus object-lens and lens-detector)
D = [d1 + T/2;T*ones(N-1,1);T/2 + d2 + dd(i)];

% Calculate the propagation parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);

% Vignetting
E1 = E0.*sqrt(exp(-x0.^2./(2*sigma_v.^2)));

% Propagate to the CRL exit plane
[E1,x1] = propFrFT1(E1,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,sum(D(1:end-1)));

% Pupil function
E1 = E1.*sqrt(exp(-x1.^2./(2*sigma_p.^2)));

% Propagate to image plane
[E2,x2] = propFrFT1(E1,x1,Inf,Rp(end),sm(end),sp(end),a(end),lambda,D(end));

% Save the results
E(:,i) = E2;
X(:,i) = x2;

% Fit the wavefront
R2(i) = wavefrontCurvature(x2,E2,lambda);

% Fit the beam waist
fun = @(w,x) exp(-2*x.^2./((w*1e-7).^2));
y = abs(E2).^2./(max(abs(E2)).^2);
ff = fit(x2,y,fun,'StartPoint',3);
w0(i) = ff.w*1e-7;

end

% Plot the curvature versus displacement
figure;
plot(dd,1./R2);

% Fit the curvature versus displacement
funR = @(w0,x0,x) 1./((x - x0).*(1 + (pi.*w0.^2./(lambda.*(x - x0))).^2));
ffR = fit(dd,1./R2,funR,'StartPoint',[min(w0) 1e-10],'Lower',[1e-10 -d2]);

% Calculate the d2 displacement
d2cc = ffR.x0;
w0R = ffR.w0;

% Plot the fitting result
figure;
plot(ffR,dd,1./R2);

% Plot the beam waist versus displacement
figure;
plot(dd,w0);

% Fit the beam waist versus displacement
funw = @(w0,x0,x) 1e6*w0.*sqrt(1 + ((x - x0).*lambda./(pi.*w0.^2)).^2);
ffw = fit(dd,1e6*w0,funw,'StartPoint',[min(w0) 0]);

% Get the fitted parameters
w0w = ffw.w0;

% Plot the fitting result
figure;
plot(ffw,dd,1e6*w0);

% Calculate the intensity
I = abs(E).^2;

% Calculate the phase
P = angle(E);
P = unwrap(P,[],1);
P = P - P(m/2 + 1,:);
%P = bsxfun(@minus,P,P(m/2 + 1,:));

% Plot the intensity map
figure;
pcolor(dd,X,I);
shading flat;
set(gca,'YLim',[-2e-6 2e-6]);

% Plot the phase map
figure;
pcolor(dd,X,P);
shading flat;
set(gca,'YLim',[-2e-6 2e-6],'CLim',[-pi pi]);

% Plot both maps
figure;
set(gcf,'Position',[50 50 1600 800]);
subplot(1,2,1);
pcolor(100*dd,1e6*X,I);
shading flat;
set(gca,'YLim',[-2 2],'XTick',-2:2,'YTick',-2:2,'FontSize',20,'OuterPosition',[0 0 0.5 1]);
xlabel('z [cm]');
ylabel('x [\mum]');
text(-1.9,1.8,'(a)','FontSize',30,'Color','white');
subplot(1,2,2);
pcolor(100*dd,1e6*X,P);
shading flat;
set(gca,'YLim',[-2 2],'CLim',[-pi pi],'XTick',-2:2,'YTick',-2:2,'FontSize',20,'OuterPosition',[0.5 0 0.5 1]);
xlabel('z [cm]');
ylabel('x [\mum]');
text(-1.9,1.8,'(b)','FontSize',30,'Color','white');

% Save the figure
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('Intensity_Phase_Map.pdf','-dpdf','-r0');
    print('Intensity_Phase_Map.png','-dpng','-r600');
end



% figure;
% imagesc((dd + 0.46e-3)*1e6,X(:,101),P,[-0.25e-3 0.25e-3]);
% set(gca,'YLim',[-1.5e-6 1.5e-6]);

