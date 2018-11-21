% Initialization
clear;
close all;
clc;


%% Load data
% Load a reconstruction
%load('Reconstruction_Ptycho_5_2018_07_20_15_41_3839.mat');
%load('Reconstruction_Ptycho_5_2018_07_20_17_38_3857.mat');
load('Reconstruction_Ptycho_5_2018_07_20_18_49_3809.mat');

% Generate 3D animation?
sav = 0;


%% Process the data
% Get the number of pixels
mx = size(rho,2);
my = size(rho,1);
m_omega = size(rho,3);

% Shift the object close to the center
cx = sum(sum(sum(abs(rho).*(1:mx))))./sum(abs(rho(:)));
cy = sum(sum(sum(abs(rho).*(1:my).')))./sum(abs(rho(:)));
cz = sum(sum(sum(abs(rho).*permute(1:m_omega,[3 1 2]))))./sum(abs(rho(:)));
rho = circshift(rho,round([my/2-cy mx/2-cx m_omega/2-cz]));

% Generate the field
field = padarray(rho,[(ny-my)/2 (nx-mx)/2 (n_omega-m_omega)/2]);
field = ifftshift(fftn(fftshift(field)));

% Determine the center diffraction spot
[~,cx] = max(max(abs(field(:,:,n_omega/2 + 1)),[],1),[],2);
[~,cy] = max(max(abs(field(:,:,n_omega/2 + 1)),[],2),[],1);
cx = cx - nx/2;
cy = cy - ny/2;

% Correct the phase
x = ((-mx/2):(mx/2-1))./nx;
y = (((-my/2):(my/2-1))./ny).';
rho = rho.*exp(-2.*pi.*1i.*cx.*x).*exp(-2.*pi.*1i.*cy.*y);
rho = rho.*exp(-1i.*angle(rho(my/2+1,mx/2+1,m_omega/2+1)));

% Set experimental parameters
M = 1.409;
D = 1.661/M;
det_dx = 55e-6;
E = 8e3;
lambda = E2lambda(E);
k = 2.*pi./lambda;
delta_omega = 0.003*pi/180;
a = 3.9242e-10;
q0 = 2*pi/a*sqrt(3);
th = asin(q0.*lambda./(4.*pi));
delta = 5.25654177e-5;
mu = 1./2.36647e-6;

% Calculate the q-range
[delta_q1,delta_q2,delta_q3,delta_r1,delta_r2,delta_r3] = q_range2(D,det_dx,th,lambda,delta_omega,nx,ny,n_omega);

% Shear the data
shft = shift_amount(delta_r1,delta_r3,th);
%object_avg = shear(object_avg,shft);
rho = shear_interp(rho,shft);

% Flip the object
%object_avg = conj(flip(flip(flip(object_avg,1),2),3));

% Generate the real space axes
delta_rx = delta_r2;
delta_ry = delta_r1;
delta_rz = delta_r3.*cos(th);
rx = ((-mx/2):(mx/2-1)).'.*delta_rx;
ry = flip(((-my/2):(my/2-1)).'.*delta_ry);
rz = ((-m_omega/2):(m_omega/2-1)).'.*delta_rz;

% Correct for attenuation by propagation through Pt
i0 = 121;
T = 2.*(ry(i0) - ry)./sin(th);
T(T < 0) = 0;
rho = rho.*exp(mu.*T);

% Plot the amplitude and phase
Slicer(abs(rho),'name','Amplitude','spacing',[delta_rx delta_ry delta_rz],'displayRange',[0 max(abs(rho(:)))]);
Slicer(angle(rho),'name','Phase','spacing',[delta_rx delta_ry delta_rz],'displayRange',[-pi pi]);

% Plot the central slice
s = 1e9;
ii = 137;
cmap = hsv(256);
rgb = complex2rgb(rho(:,:,ii)./max(max(abs(rho(:,:,ii)))),cmap);
figure;
image(s.*rx,s.*ry,rgb);
axis equal tight;
set(gca,'YDir','normal','XLim',[-500 500],'YLim',[-300 300],'XTick',-400:200:400,'YTick',-200:100:200,'FontSize',16);
xlabel('x [nm]');
ylabel('y [nm]');

% Make 3D plot of the raw reconstructed object
figure;
pp = patch(isosurface(s.*rx,s.*ry,s.*rz,abs(rho),0.1*max(abs(rho(:)))));
isonormals(s.*rx,s.*ry,s.*rz,abs(rho),pp);
pp.FaceColor = 'red';
pp.EdgeColor = 'none';
daspect([1 1 1]);
view(3);
camlight(45,45);
camlight(-135,-45);
lighting gouraud;
xlabel('x [nm]');
ylabel('y [nm]');
zlabel('z [nm]');

% Generate 3D rotation of object
if sav
    cx = s.*sum(sum(sum(rx.'.*abs(rho))))./sum(abs(rho(:)));
    cy = s.*sum(sum(sum(ry.*abs(rho))))./sum(abs(rho(:)));
    cz = s.*sum(sum(sum(permute(rz,[2 3 1]).*abs(rho))))./sum(abs(rho(:)));
    n = 360;
    mov(n) = struct('cdata',[],'colormap',[]);
    A = gca;
    A.CameraUpVector = [0 1 0];
    A.CameraTarget = [cx cy cz];
    A.CameraViewAngle = 28/10;
    axis off;
    set(gcf,'Position',[460 100 1000 800],'Color',[1 1 1]);
    for i = 1:n
        alpha = (i-1)/n*360;
        A.CameraPosition = [-sind(alpha)*1000 500 cosd(alpha)*1000]*10;
        drawnow;
        mov(i) = getframe(gcf);
    end
    vid = VideoWriter('Virtual_Ptycho_5.mp4','MPEG-4');
    vid.FrameRate = 60;
    open(vid);
    writeVideo(vid,mov);
    close(vid);
end

% Save the processed data
save('Corrected_Virtual_Ptycho_5.mat','rho','rx','ry','rz');


