% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the data sets
tic;
Ef1 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_field.bin'],[ny nx n_omega]);
toc;
tic;
Ef2 = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_2_field.bin'],[ny nx n_omega]);
toc;
tic;
El1 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_field.bin'],[ny nx n_omega]);
toc;
tic;
El2 = load_binary([p 'Reconstruction_Lens_Shrinkwrap_2_field.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures\';
sav = 0;


%% Process the data
% Shift parameters
gamma = acos(dot(q1v,q3v)./(norm(q1v).*norm(q3v)));
shftq = tan(gamma - pi/2)*delta_qy/delta_qz;
qy = sin(gamma).*qy;
qz = ((-n_omega/2):(n_omega/2 - 1)).'.*delta_q3;

% Shift the vector space to become orthogonal
tic;
for i = 1:ny
    Ef1(i,:,:) = circshift(Ef1(i,:,:),round((i - ny/2)*shftq),3);
    Ef2(i,:,:) = circshift(Ef2(i,:,:),round((i - ny/2)*shftq),3);
    El1(i,:,:) = circshift(El1(i,:,:),round((i - ny/2)*shftq),3);
    El2(i,:,:) = circshift(El2(i,:,:),round((i - ny/2)*shftq),3);
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');
toc;


%% Plot the data
% Extract slices to plot
tic;
Ef1xy = log10(abs(Ef1(:,:,n_omega/2 + 1)));
Ef1zy = log10(squeeze(abs(Ef1(:,nx/2 + 1,:))));
Ef1xz = log10(squeeze(abs(Ef1(ny/2 + 1,:,:))).');
Ef2xy = log10(abs(Ef2(:,:,n_omega/2 + 1)));
Ef2zy = log10(squeeze(abs(Ef2(:,nx/2 + 1,:))));
Ef2xz = log10(squeeze(abs(Ef2(ny/2 + 1,:,:))).');
El1xy = log10(abs(El1(:,:,n_omega/2 + 1)));
El1zy = log10(squeeze(abs(El1(:,nx/2 + 1,:))));
El1xz = log10(squeeze(abs(El1(ny/2 + 1,:,:))).');
El2xy = log10(abs(El2(:,:,n_omega/2 + 1)));
El2zy = log10(squeeze(abs(El2(:,nx/2 + 1,:))));
El2xz = log10(squeeze(abs(El2(ny/2 + 1,:,:))).');
toc;

% Scale slices
tic;
Ef1xy = (Ef1xy + 0.5)./3.5;
Ef1zy = (Ef1zy + 0.5)./3.5;
Ef1xz = (Ef1xz + 0.5)./3.5;
Ef2xy = (Ef2xy + 0.5)./3.5;
Ef2zy = (Ef2zy + 0.5)./3.5;
Ef2xz = (Ef2xz + 0.5)./3.5;
El1xy = (El1xy + 0.5)./3.5;
El1zy = (El1zy + 0.5)./3.5;
El1xz = (El1xz + 0.5)./3.5;
El2xy = (El2xy + 0.5)./3.5;
El2zy = (El2zy + 0.5)./3.5;
El2xz = (El2xz + 0.5)./3.5;
toc;

% Remove negative values
tic;
Ef1xy(Ef1xy < 0) = 0;
Ef1zy(Ef1zy < 0) = 0;
Ef1xz(Ef1xz < 0) = 0;
Ef2xy(Ef2xy < 0) = 0;
Ef2zy(Ef2zy < 0) = 0;
Ef2xz(Ef2xz < 0) = 0;
El1xy(El1xy < 0) = 0;
El1zy(El1zy < 0) = 0;
El1xz(El1xz < 0) = 0;
El2xy(El2xy < 0) = 0;
El2zy(El2zy < 0) = 0;
El2xz(El2xz < 0) = 0;
toc;

% Add phase
tic;
Ef1xy = Ef1xy.*exp(1i*angle(Ef1(:,:,n_omega/2 + 1)));
Ef1zy = Ef1zy.*exp(1i*squeeze(angle(Ef1(:,nx/2 + 1,:))));
Ef1xz = Ef1xz.*exp(1i*squeeze(angle(Ef1(ny/2 + 1,:,:))).');
Ef2xy = Ef2xy.*exp(1i*angle(Ef2(:,:,n_omega/2 + 1)));
Ef2zy = Ef2zy.*exp(1i*squeeze(angle(Ef2(:,nx/2 + 1,:))));
Ef2xz = Ef2xz.*exp(1i*squeeze(angle(Ef2(ny/2 + 1,:,:))).');
El1xy = El1xy.*exp(1i*angle(El1(:,:,n_omega/2 + 1)));
El1zy = El1zy.*exp(1i*squeeze(angle(El1(:,nx/2 + 1,:))));
El1xz = El1xz.*exp(1i*squeeze(angle(El1(ny/2 + 1,:,:))).');
El2xy = El2xy.*exp(1i*angle(El2(:,:,n_omega/2 + 1)));
El2zy = El2zy.*exp(1i*squeeze(angle(El2(:,nx/2 + 1,:))));
El2xz = El2xz.*exp(1i*squeeze(angle(El2(ny/2 + 1,:,:))).');
toc;

% Generate color image
tic;
cmap = hsv(256);
Ef1xy = complex2rgb(Ef1xy,cmap);
Ef1zy = complex2rgb(Ef1zy,cmap);
Ef1xz = complex2rgb(Ef1xz,cmap);
Ef2xy = complex2rgb(Ef2xy,cmap);
Ef2zy = complex2rgb(Ef2zy,cmap);
Ef2xz = complex2rgb(Ef2xz,cmap);
El1xy = complex2rgb(El1xy,cmap);
El1zy = complex2rgb(El1zy,cmap);
El1xz = complex2rgb(El1xz,cmap);
El2xy = complex2rgb(El2xy,cmap);
El2zy = complex2rgb(El2zy,cmap);
El2xz = complex2rgb(El2xz,cmap);
toc;

% Set plotting font size
F = 16;
s = 1e-9;
u = '[nm^{-1}]';
xt = -0.4:0.2:0.4;
yt = -0.4:0.2:0.4;
zt = -0.2:0.2:0.2;

% Plot XY planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,Ef1xy);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef1xy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef1xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,Ef2xy);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef2xy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef2xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,El1xy);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El1xy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El1xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,El2xy);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El2xy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El2xy.png'],'-dpng','-r600');
end

% Plot ZY planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qz,s*qy,Ef1zy);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef1zy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef1zy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qz,s*qy,Ef2zy);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef2zy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef2zy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qz,s*qy,El1zy);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El1zy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El1zy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qz,s*qy,El2zy);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El2zy.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El2zy.png'],'-dpng','-r600');
end

% Plot XZ planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qz,Ef1xz);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef1xz.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef1xz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qz,Ef2xz);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef2xz.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef2xz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qz,El1xz);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El1xz.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El1xz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qz,El2xz);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El2xz.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El2xz.png'],'-dpng','-r600');
end

% Make a combined free space 1,000,000 photons figure
F = 12;
l = 500;
ll = 72;
b = 15;
sy = (max(qy) - min(qy))/(max(qx) - min(qx));
sz = (max(qz) - min(qz))/(max(qx) - min(qx));
figure;
set(gcf,'Position',[400 100 2*ll+l+sz*l+b 2*ll+sz*l+sy*l+b]);
axes('Units','pixels','Position',[ll ll+sz*l+b l sy*l]);
image(s*qx,s*qy,Ef1xy);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-0.5,0.45,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+sz*l+b sz*l sy*l]);
image(s*qz,s*qy,Ef1zy);
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-0.34,0.45,'(b)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll l sz*l]);
image(s*qx,s*qz,Ef1xz);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.5,0.32,'(c)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef1.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef1.png'],'-dpng','-r600');
end

% Make a combined free space 50,000 photons figure
figure;
set(gcf,'Position',[400 100 2*ll+l+sz*l+b 2*ll+sz*l+sy*l+b]);
axes('Units','pixels','Position',[ll ll+sz*l+b l sy*l]);
image(s*qx,s*qy,Ef2xy);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-0.5,0.45,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+sz*l+b sz*l sy*l]);
image(s*qz,s*qy,Ef2zy);
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-0.34,0.45,'(b)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll l sz*l]);
image(s*qx,s*qz,Ef2xz);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.5,0.32,'(c)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_Ef2.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_Ef2.png'],'-dpng','-r600');
end

% Make a combined lens 1,000,000 photons figure
figure;
set(gcf,'Position',[400 100 2*ll+l+sz*l+b 2*ll+sz*l+sy*l+b]);
axes('Units','pixels','Position',[ll ll+sz*l+b l sy*l]);
image(s*qx,s*qy,El1xy);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-0.5,0.45,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+sz*l+b sz*l sy*l]);
image(s*qz,s*qy,El1zy);
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-0.34,0.45,'(b)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll l sz*l]);
image(s*qx,s*qz,El1xz);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.5,0.32,'(c)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El1.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El1.png'],'-dpng','-r600');
end

% Make a combined lens 50,000 photons figure
figure;
set(gcf,'Position',[400 100 2*ll+l+sz*l+b 2*ll+sz*l+sy*l+b]);
axes('Units','pixels','Position',[ll ll+sz*l+b l sy*l]);
image(s*qx,s*qy,El2xy);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-0.5,0.45,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+sz*l+b sz*l sy*l]);
image(s*qz,s*qy,El2zy);
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-0.34,0.45,'(b)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll l sz*l]);
image(s*qx,s*qz,El2xz);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.5,0.32,'(c)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Rec_Field_El2.pdf'],'-dpdf','-r0');
    print([ps 'Rec_Field_El2.png'],'-dpng','-r600');
end


