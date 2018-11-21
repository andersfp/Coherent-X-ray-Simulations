% Initialization
clear;
close all;
clc;


%% Load data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the data sets
tic;
Ef = load_binary([p 'Ef.bin'],[ny nx n_omega]);
toc;
tic;
El = load_binary([p 'El.bin'],[ny nx n_omega]);
toc;
tic;
If1 = load_binary([p 'If1.bin'],[ny nx n_omega]);
toc;
tic;
If2 = load_binary([p 'If2.bin'],[ny nx n_omega]);
toc;
tic;
Il1 = load_binary([p 'Il1.bin'],[ny nx n_omega]);
toc;
tic;
Il2 = load_binary([p 'Il2.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
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
    Ef(i,:,:) = circshift(Ef(i,:,:),round((i - ny/2)*shftq),3);
    El(i,:,:) = circshift(El(i,:,:),round((i - ny/2)*shftq),3);
    If1(i,:,:) = circshift(If1(i,:,:),round((i - ny/2)*shftq),3);
    If2(i,:,:) = circshift(If2(i,:,:),round((i - ny/2)*shftq),3);
    Il1(i,:,:) = circshift(Il1(i,:,:),round((i - ny/2)*shftq),3);
    Il2(i,:,:) = circshift(Il2(i,:,:),round((i - ny/2)*shftq),3);
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');
toc;

% Get the intensities from non-rounded electric field
If = abs(Ef).^2;
Il = abs(El).^2;

% Normalize the intensities
Il = 1e6*Il./max(If(:));
If = 1e6*If./max(If(:));


%% Plot the data
% Extract slices to plot
Ifxy = log10(If(:,:,n_omega/2 + 1));
Ifzy = log10(squeeze(If(:,nx/2 + 1,:)));
Ifxz = log10(squeeze(If(ny/2 + 1,:,:)).');
If1xy = log10(If1(:,:,n_omega/2 + 1));
If1zy = log10(squeeze(If1(:,nx/2 + 1,:)));
If1xz = log10(squeeze(If1(ny/2 + 1,:,:)).');
If2xy = log10(If2(:,:,n_omega/2 + 1));
If2zy = log10(squeeze(If2(:,nx/2 + 1,:)));
If2xz = log10(squeeze(If2(ny/2 + 1,:,:)).');
Ilxy = log10(Il(:,:,n_omega/2 + 1));
Ilzy = log10(squeeze(Il(:,nx/2 + 1,:)));
Ilxz = log10(squeeze(Il(ny/2 + 1,:,:)).');
Il1xy = log10(Il1(:,:,n_omega/2 + 1));
Il1zy = log10(squeeze(Il1(:,nx/2 + 1,:)));
Il1xz = log10(squeeze(Il1(ny/2 + 1,:,:)).');
Il2xy = log10(Il2(:,:,n_omega/2 + 1));
Il2zy = log10(squeeze(Il2(:,nx/2 + 1,:)));
Il2xz = log10(squeeze(Il2(ny/2 + 1,:,:)).');

% Set plotting font size
F = 16;
s = 1e-9;
u = '[nm^{-1}]';
xt = -0.6:0.3:0.6;
yt = -0.6:0.3:0.6;
zt = -0.3:0.3:0.3;

% Plot XY planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,Ifxy,[-1 6]);
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
    print([ps 'Ifxy.pdf'],'-dpdf','-r0');
    print([ps 'Ifxy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,If1xy,[-1 6]);
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
    print([ps 'If1xy.pdf'],'-dpdf','-r0');
    print([ps 'If1xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,If2xy,[-1 6]);
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
    print([ps 'If2xy.pdf'],'-dpdf','-r0');
    print([ps 'If2xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,Ilxy,[-1 6]);
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
    print([ps 'Ilxy.pdf'],'-dpdf','-r0');
    print([ps 'Ilxy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,Il1xy,[-1 6]);
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
    print([ps 'Il1xy.pdf'],'-dpdf','-r0');
    print([ps 'Il1xy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qy,Il2xy,[-1 6]);
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
    print([ps 'Il2xy.pdf'],'-dpdf','-r0');
    print([ps 'Il2xy.png'],'-dpng','-r600');
end

% Plot ZY planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qz,s*qy,Ifzy,[-1 6]);
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
    print([ps 'Ifzy.pdf'],'-dpdf','-r0');
    print([ps 'Ifzy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qz,s*qy,If1zy,[-1 6]);
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
    print([ps 'If1zy.pdf'],'-dpdf','-r0');
    print([ps 'If1zy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qz,s*qy,If2zy,[-1 6]);
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
    print([ps 'If2zy.pdf'],'-dpdf','-r0');
    print([ps 'If2zy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qz,s*qy,Ilzy,[-1 6]);
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
    print([ps 'Ilzy.pdf'],'-dpdf','-r0');
    print([ps 'Ilzy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qz,s*qy,Il1zy,[-1 6]);
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
    print([ps 'Il1zy.pdf'],'-dpdf','-r0');
    print([ps 'Il1zy.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qz,s*qy,Il2zy,[-1 6]);
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
    print([ps 'Il2zy.pdf'],'-dpdf','-r0');
    print([ps 'Il2zy.png'],'-dpng','-r600');
end

% Plot XZ planes
figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qz,Ifxz,[-1 6]);
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
    print([ps 'Ifxz.pdf'],'-dpdf','-r0');
    print([ps 'Ifxz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qz,If1xz,[-1 6]);
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
    print([ps 'If1xz.pdf'],'-dpdf','-r0');
    print([ps 'If1xz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qz,If2xz,[-1 6]);
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
    print([ps 'If2xz.pdf'],'-dpdf','-r0');
    print([ps 'If2xz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qz,Ilxz,[-1 6]);
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
    print([ps 'Ilxz.pdf'],'-dpdf','-r0');
    print([ps 'Ilxz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qz,Il1xz,[-1 6]);
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
    print([ps 'Il1xz.pdf'],'-dpdf','-r0');
    print([ps 'Il1xz.png'],'-dpng','-r600');
end

figure;
set(gcf,'Position',[460 100 1000 1000]);
imagesc(s*qx,s*qz,Il2xz,[-1 6]);
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
    print([ps 'Il2xz.pdf'],'-dpdf','-r0');
    print([ps 'Il2xz.png'],'-dpng','-r600');
end

% Make a combined free space figure
F = 12;
figure;
set(gcf,'Position',[460 100 1000 1000]);
axes('OuterPosition',[0 2/3 1/3 1/3]);
imagesc(s*qx,s*qy,Ifxy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(a)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 2/3 1/3 1/3]);
imagesc(s*qx,s*qy,If1xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(b)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 2/3 1/3 1/3]);
imagesc(s*qx,s*qy,If2xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(c)','FontSize',20,'Color','white');
axes('OuterPosition',[0 1/3 1/3 1/3]);
imagesc(s*qz,s*qy,Ifzy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(d)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 1/3 1/3 1/3]);
imagesc(s*qz,s*qy,If1zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(e)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 1/3 1/3 1/3]);
imagesc(s*qz,s*qy,If2zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(f)','FontSize',20,'Color','white');
axes('OuterPosition',[0 0 1/3 1/3]);
imagesc(s*qx,s*qz,Ifxz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(g)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 0 1/3 1/3]);
imagesc(s*qx,s*qz,If1xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(h)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 0 1/3 1/3]);
imagesc(s*qx,s*qz,If2xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(i)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'If.pdf'],'-dpdf','-r0');
    print([ps 'If.png'],'-dpng','-r600');
end

% Make a combined lens figure
figure;
set(gcf,'Position',[460 100 1000 1000]);
axes('OuterPosition',[0 2/3 1/3 1/3]);
imagesc(s*qx,s*qy,Ilxy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(a)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 2/3 1/3 1/3]);
imagesc(s*qx,s*qy,Il1xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(b)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 2/3 1/3 1/3]);
imagesc(s*qx,s*qy,Il2xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(c)','FontSize',20,'Color','white');
axes('OuterPosition',[0 1/3 1/3 1/3]);
imagesc(s*qz,s*qy,Ilzy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(d)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 1/3 1/3 1/3]);
imagesc(s*qz,s*qy,Il1zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(e)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 1/3 1/3 1/3]);
imagesc(s*qz,s*qy,Il2zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(f)','FontSize',20,'Color','white');
axes('OuterPosition',[0 0 1/3 1/3]);
imagesc(s*qx,s*qz,Ilxz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(g)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 0 1/3 1/3]);
imagesc(s*qx,s*qz,Il1xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(h)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 0 1/3 1/3]);
imagesc(s*qx,s*qz,Il2xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(i)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Il.pdf'],'-dpdf','-r0');
    print([ps 'Il.png'],'-dpng','-r600');
end

% Make a combined xy plane figure
figure;
set(gcf,'Position',[220 100 1530 950]);
axes('OuterPosition',[0 1/2 1/3 1/2]);
imagesc(s*qx,s*qy,Ifxy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(a)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 1/2 1/3 1/2]);
imagesc(s*qx,s*qy,If1xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(b)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 1/2 1/3 1/2]);
imagesc(s*qx,s*qy,If2xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(c)','FontSize',20,'Color','white');
axes('OuterPosition',[0 0 1/3 1/2]);
imagesc(s*qx,s*qy,Ilxy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(d)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 0 1/3 1/2]);
imagesc(s*qx,s*qy,Il1xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(e)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 0 1/3 1/2]);
imagesc(s*qx,s*qy,Il2xy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);
text(-0.81,0.78,'(f)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Ixy.pdf'],'-dpdf','-r0');
    print([ps 'Ixy.png'],'-dpng','-r600');
end

% Make a combined zy plane figure
figure;
set(gcf,'Position',[370 100 1180 1000]);
axes('OuterPosition',[0 1/2 1/3 1/2]);
imagesc(s*qz,s*qy,Ifzy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(a)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 1/2 1/3 1/2]);
imagesc(s*qz,s*qy,If1zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(b)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 1/2 1/3 1/2]);
imagesc(s*qz,s*qy,If2zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(c)','FontSize',20,'Color','white');
axes('OuterPosition',[0 0 1/3 1/2]);
imagesc(s*qz,s*qy,Ilzy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(d)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 0 1/3 1/2]);
imagesc(s*qz,s*qy,Il1zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(e)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 0 1/3 1/2]);
imagesc(s*qz,s*qy,Il2zy,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt);
text(-0.34,0.78,'(f)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Izy.pdf'],'-dpdf','-r0');
    print([ps 'Izy.png'],'-dpng','-r600');
end

% Make a combined xz plane figure
figure;
set(gcf,'Position',[110 170 1700 800]);
axes('OuterPosition',[0 1/2 1/3 1/2]);
imagesc(s*qx,s*qz,Ifxz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(a)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 1/2 1/3 1/2]);
imagesc(s*qx,s*qz,If1xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(b)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 1/2 1/3 1/2]);
imagesc(s*qx,s*qz,If2xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(c)','FontSize',20,'Color','white');
axes('OuterPosition',[0 0 1/3 1/2]);
imagesc(s*qx,s*qz,Ilxz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(d)','FontSize',20,'Color','white');
axes('OuterPosition',[1/3 0 1/3 1/2]);
imagesc(s*qx,s*qz,Il1xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(e)','FontSize',20,'Color','white');
axes('OuterPosition',[2/3 0 1/3 1/2]);
imagesc(s*qx,s*qz,Il2xz,[-1 6]);
axis equal tight;
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(f)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Ixz.pdf'],'-dpdf','-r0');
    print([ps 'Ixz.png'],'-dpng','-r600');
end

% Make a combined free space infinite photons figure
l = 500;
ll = 72;
b = 15;
sy = (max(qy) - min(qy))/(max(qx) - min(qx));
sz = (max(qz) - min(qz))/(max(qx) - min(qx));
figure;
set(gcf,'Position',[400 100 2*ll+l+sz*l+b 2*ll+sz*l+sy*l+b]);
axes('Units','pixels','Position',[ll ll l sz*l]);
imagesc(s*qx,s*qz,Ifxz,[-1 6]);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(c)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll+sz*l+b l sy*l]);
imagesc(s*qx,s*qy,Ifxy,[-1 6]);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-0.81,0.78,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+sz*l+b sz*l sy*l]);
imagesc(s*qz,s*qy,Ifzy,[-1 6]);
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-0.34,0.78,'(b)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'If_Slices.pdf'],'-dpdf','-r0');
    print([ps 'If_Slices.png'],'-dpng','-r600');
end

% Make a combined lens infinite photons figure
figure;
set(gcf,'Position',[400 100 2*ll+l+sz*l+b 2*ll+sz*l+sy*l+b]);
axes('Units','pixels','Position',[ll ll l sz*l]);
imagesc(s*qx,s*qz,Ilxz,[-1 6]);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_z ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',zt);
text(-0.81,0.30,'(c)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll ll+sz*l+b l sy*l]);
imagesc(s*qx,s*qy,Ilxy,[-1 6]);
colormap gray;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt,'XAxisLocation','top');
text(-0.81,0.78,'(a)','FontSize',20,'Color','white');
axes('Units','pixels','Position',[ll+l+b ll+sz*l+b sz*l sy*l]);
imagesc(s*qz,s*qy,Ilzy,[-1 6]);
colormap gray;
xlabel(['q_z ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',zt,'YTick',yt,'XAxisLocation','top','YAxisLocation','right');
text(-0.34,0.78,'(b)','FontSize',20,'Color','white');
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print([ps 'Il_Slices.pdf'],'-dpdf','-r0');
    print([ps 'Il_Slices.png'],'-dpng','-r600');
end



rgb = cat(3,Ifxy,Ilxy,zeros(size(Ifxy)));
rgb = rgb/6;
rgb(rgb < 0) = 0;
figure;
image(s*qx,s*qy,rgb);
axis equal tight;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);

figure;
imagesc(s*qx,s*qy,rgb(:,:,2)./rgb(:,:,1));
axis equal tight;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);

qa = 4*pi/lambda*sin(sigma_a);
rgb2 = cat(3,Ifxy.*exp(-(qx.'.^2 + qy.^2)./(2.*qa.^2)),Ilxy,zeros(size(Ifxy)));
rgb2 = rgb2/6;
rgb2(rgb2 < 0) = 0;
figure;
image(s*qx,s*qy,rgb2);
axis equal tight;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);

figure;
imagesc(s*qx,s*qy,rgb2(:,:,2)./rgb2(:,:,1));
axis equal tight;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);



figure;
imagesc(s*qx,s*qy,rgb2(:,:,1)./rgb(:,:,1));
axis equal tight;
xlabel(['q_x ' u]);
ylabel(['q_y ' u]);
set(gca,'FontSize',F,'YDir','normal','XTick',xt,'YTick',yt);

