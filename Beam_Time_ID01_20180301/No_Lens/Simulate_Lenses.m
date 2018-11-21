% Initialization
clear;
close all;
clc;


%% Load data
% Load a reconstruction
load('Reconstruction_Center_2018_03_13_23_38_300.mat');

% Set the pixel number
nx = 512;
ny = 512;
n_omega = 280;


%% Process the data
% Get the number of pixels
mx = size(object,2);
my = size(object,1);
m_omega = size(object,3);

% Generate the field
field = padarray(object_avg,[(ny-my)/2 (nx-mx)/2 (n_omega-m_omega)/2]);
field = fftshift(fftn(ifftshift(field)));

% Determine the center diffraction spot
[~,cx] = max(max(abs(field(:,:,n_omega/2 + 1)),[],1),[],2);
[~,cy] = max(max(abs(field(:,:,n_omega/2 + 1)),[],2),[],1);
cx = cx - nx/2;
cy = cy - ny/2;

% Generate the two lens pupil functions
wv = 25.2; % 24.2*0.6887
wr = 7.85; % 5.42*0.8993
qx = (-nx/2):(nx/2 - 1);
qy = ((-ny/2):(ny/2 - 1)).';
gv = gaussRMS(qx - cx,wv).*gaussRMS(qy - cy,wv);
gr = gaussRMS(qx - cx,wr).*gaussRMS(qy - cy,wr);

% Simulate the diffraction patterns
fieldv = sqrt(gv).*field;
fieldr = sqrt(gr).*field;

% Calculate the objects
objectv = fftshift(ifftn(ifftshift(fieldv)));
objectr = fftshift(ifftn(ifftshift(fieldr)));

% Remove padding from simulated objects
objectv = objectv((ny-my)/2+1:end-(ny-my)/2,(nx-mx)/2+1:end-(nx-mx)/2,(n_omega-m_omega)/2+1:end-(n_omega-m_omega)/2);
objectr = objectr((ny-my)/2+1:end-(ny-my)/2,(nx-mx)/2+1:end-(nx-mx)/2,(n_omega-m_omega)/2+1:end-(n_omega-m_omega)/2);

% Correct the phase
x = ((-mx/2):(mx/2-1))./nx;
y = (((-my/2):(my/2-1))./ny).';
object_avg = object_avg.*exp(-2.*pi.*1i.*cx.*x).*exp(-2.*pi.*1i.*cy.*y);
object_avg = object_avg.*exp(-1i.*angle(object_avg(my/2+1,mx/2+1,m_omega/2+1)));
objectv = objectv.*exp(-2.*pi.*1i.*cx.*x).*exp(-2.*pi.*1i.*cy.*y);
objectv = objectv.*exp(-1i.*angle(objectv(my/2+1,mx/2+1,m_omega/2+1)));
objectr = objectr.*exp(-2.*pi.*1i.*cx.*x).*exp(-2.*pi.*1i.*cy.*y);
objectr = objectr.*exp(-1i.*angle(objectr(my/2+1,mx/2+1,m_omega/2+1)));

% Set experimental parameters
D = 1.632;
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
object_avg = shear_interp(object_avg,shft);
objectv = shear_interp(objectv,shft);
objectr = shear_interp(objectr,shft);



% Flip the object
object_avg = conj(flip(flip(flip(object_avg,1),2),3));
objectv = conj(flip(flip(flip(objectv,1),2),3));
objectr = conj(flip(flip(flip(objectr,1),2),3));

% Generate the real space axes
delta_rx = delta_r2;
delta_ry = delta_r1;
delta_rz = delta_r3.*cos(th);
rx = ((-mx/2):(mx/2-1)).'.*delta_rx;
ry = flip(((-my/2):(my/2-1)).'.*delta_ry);
rz = ((-m_omega/2):(m_omega/2-1)).'.*delta_rz;

% Correct for attenuation by propagation through Pt
i0 = 30;
T = 2.*(ry(i0) - ry)./sin(th);
T(T < 0) = 0;
object_avg = object_avg.*exp(mu.*T);
objectv = objectv.*exp(mu.*T);
objectr = objectr.*exp(mu.*T);

% Plot the amplitude and phase
Slicer(abs(object_avg),'name','Amplitude','spacing',[delta_rx delta_ry delta_rz],'displayRange',[0 max(abs(object_avg(:)))]);
Slicer(angle(object_avg),'name','Phase','spacing',[delta_rx delta_ry delta_rz],'displayRange',[-pi pi]);
Slicer(abs(objectv),'name','Amplitude','spacing',[delta_rx delta_ry delta_rz],'displayRange',[0 max(abs(objectv(:)))]);
Slicer(angle(objectv),'name','Phase','spacing',[delta_rx delta_ry delta_rz],'displayRange',[-pi pi]);
Slicer(abs(objectr),'name','Amplitude','spacing',[delta_rx delta_ry delta_rz],'displayRange',[0 max(abs(objectr(:)))]);
Slicer(angle(objectr),'name','Phase','spacing',[delta_rx delta_ry delta_rz],'displayRange',[-pi pi]);

% Calculate central RGB slices
s = 1e9;
ii = 34;
cmap = hsv(256);
rgb = complex2rgb(object_avg(:,:,ii)./max(max(abs(object_avg(:,:,ii)))),cmap);
rgbv = complex2rgb(objectv(:,:,ii)./max(max(abs(objectv(:,:,ii)))),cmap);
rgbr = complex2rgb(objectr(:,:,ii)./max(max(abs(objectr(:,:,ii)))),cmap);

% Plot the central slices
figure;
image(s.*rx,s.*ry,rgb);
axis equal tight;
set(gca,'YDir','normal','XLim',[-500 500],'YLim',[-300 300],'XTick',-400:200:400,'YTick',-200:100:200,'FontSize',16);
xlabel('x [nm]');
ylabel('y [nm]');
figure;
image(s.*rx,s.*ry,rgbv);
axis equal tight;
set(gca,'YDir','normal','XLim',[-500 500],'YLim',[-300 300],'XTick',-400:200:400,'YTick',-200:100:200,'FontSize',16);
xlabel('x [nm]');
ylabel('y [nm]');
figure;
image(s.*rx,s.*ry,rgbr);
axis equal tight;
set(gca,'YDir','normal','XLim',[-500 500],'YLim',[-300 300],'XTick',-400:200:400,'YTick',-200:100:200,'FontSize',16);
xlabel('x [nm]');
ylabel('y [nm]');

% Make 3D plot of the raw reconstructed object
figure;
pp = patch(isosurface(s.*rx,s.*ry,s.*rz,abs(object_avg),0.3*max(abs(object_avg(:)))));
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
figure;
pp = patch(isosurface(s.*rx,s.*ry,s.*rz,abs(objectv),0.3*max(abs(objectv(:)))));
isonormals(s.*rx,s.*ry,s.*rz,abs(objectv),pp);
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
figure;
pp = patch(isosurface(s.*rx,s.*ry,s.*rz,abs(objectr),0.3*max(abs(objectr(:)))));
isonormals(s.*rx,s.*ry,s.*rz,abs(objectr),pp);
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



