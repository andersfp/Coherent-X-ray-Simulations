% Initialization
clear;
close all;
clc;


%% Load data
% Load a reconstruction
load('Reconstruction_Stitched_2_2018_05_22_15_35_5268.mat');

% Generate 3D animation?
sav = 0;


%% Process the data
% Get the number of pixels
mx = size(object,2);
my = size(object,1);
m_omega = size(object,3);

% Generate the field
field = padarray(object_avg,[(ny-my)/2 (nx-mx)/2 (n_omega-m_omega)/2]);
field = ifftshift(fftn(fftshift(field)));

% Determine the center diffraction spot
[~,cx] = max(max(abs(field(:,:,n_omega/2 + 1)),[],1),[],2);
[~,cy] = max(max(abs(field(:,:,n_omega/2 + 1)),[],2),[],1);
cx = cx - nx/2;
cy = cy - ny/2;

% Correct the phase
x = ((-mx/2):(mx/2-1))./nx;
y = (((-my/2):(my/2-1))./ny).';
object_avg = object_avg.*exp(-2.*pi.*1i.*cx.*x).*exp(-2.*pi.*1i.*cy.*y);
object_avg = object_avg.*exp(-1i.*angle(object_avg(my/2+1,mx/2+1,m_omega/2+1)));

% Set experimental parameters
M = 1.48;
D = 1.662/M;
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
object_avg = shear_interp(object_avg,shft);

% Flip the object
object_avg = conj(flip(flip(flip(object_avg,1),2),3));

% Generate the real space axes
delta_rx = delta_r2;
delta_ry = delta_r1;
delta_rz = delta_r3.*cos(th);
rx = ((-mx/2):(mx/2-1)).'.*delta_rx;
ry = flip(((-my/2):(my/2-1)).'.*delta_ry);
rz = ((-m_omega/2):(m_omega/2-1)).'.*delta_rz;

% Correct for attenuation by propagation through Pt
i0 = 21;
T = 2.*(ry(i0) - ry)./sin(th);
T(T < 0) = 0;
object_avg = object_avg.*exp(mu.*T);

% Plot the amplitude and phase
Slicer(abs(object_avg),'name','Amplitude','spacing',[delta_rx delta_ry delta_rz],'displayRange',[0 max(abs(object_avg(:)))]);
Slicer(angle(object_avg),'name','Phase','spacing',[delta_rx delta_ry delta_rz],'displayRange',[-pi pi]);

% Plot the central slice
s = 1e9;
ii = 43;
cmap = hsv(256);
rgb = complex2rgb(object_avg(:,:,ii)./max(max(abs(object_avg(:,:,ii)))),cmap);
figure;
image(s.*rx,s.*ry,rgb);
axis equal tight;
set(gca,'YDir','normal','XLim',[-500 500],'YLim',[-300 300],'XTick',-400:200:400,'YTick',-200:100:200,'FontSize',16);
xlabel('x [nm]');
ylabel('y [nm]');

% Make 3D plot of the raw reconstructed object
figure;
pp = patch(isosurface(s.*rx,s.*ry,s.*rz,abs(object_avg),0.1*max(abs(object_avg(:)))));
isonormals(s.*rx,s.*ry,s.*rz,abs(object_avg),pp);
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
    cx = s.*sum(sum(sum(rx.'.*abs(object_avg))))./sum(abs(object_avg(:)));
    cy = s.*sum(sum(sum(ry.*abs(object_avg))))./sum(abs(object_avg(:)));
    cz = s.*sum(sum(sum(permute(rz,[2 3 1]).*abs(object_avg))))./sum(abs(object_avg(:)));
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
    vid = VideoWriter('Virtual_Center.mp4','MPEG-4');
    vid.FrameRate = 60;
    open(vid);
    writeVideo(vid,mov);
    close(vid);
end



