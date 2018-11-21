% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the electric field from the free space simulation
tic;
Ef = load_binary([p 'Ef.bin'],[ny nx n_omega]);
toc;

% Load the electric field from the lens-based simulation
tic;
El = load_binary([p 'El.bin'],[ny nx n_omega]);
toc;

% Load the reconstructed field from the free space simulation
tic;
Rf = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_field.bin'],[ny nx n_omega]);
toc;

% Load the reconstructed field from the lens-based simulation
tic;
Rl = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_field.bin'],[ny nx n_omega]);
toc;

% Save the images?
ps = 'C:\Users\anfils\OneDrive\DTU\PostDoc\Papers\Lensless imaging with a lens\Figures_HE\';
sav = 0;


%% Get the central diffraction slices
% Set the phase curvatures
Rpxf = 4.0000;
Rpyf = 4.0000;
Rpxl = 6.0318;
Rpyl = 6.0318;

% Compensate the phase curvatures
tic;
xd = ((-nx/2):(nx/2 - 1)).*det_dx;
yd = ((-ny/2):(ny/2 - 1)).'.*det_dx;
Ef = Ef.*exp(-1i.*pi.*xd.^2./(lambda.*Rpxf)).*exp(-1i.*pi.*yd.^2./(lambda.*Rpyf));
El = El.*exp(-1i.*pi.*xd.^2./(lambda.*Rpxl)).*exp(-1i.*pi.*yd.^2./(lambda.*Rpyl));
toc;

% Get the slices
k = n_omega/2 + 1;
ef = Ef(:,:,k);
el = El(:,:,k);
rf = Rf(:,:,k);
rl = Rl(:,:,k);

% Center the phase
ef = ef.*exp(-1i.*angle(ef(ny/2+1,nx/2+1)));
el = el.*exp(-1i.*angle(el(ny/2+1,nx/2+1)));
rf = rf.*exp(-1i.*angle(rf(ny/2+1,nx/2+1)));
rl = rl.*exp(-1i.*angle(rl(ny/2+1,nx/2+1)));

% Normalize the amplitudes
ef = ef./max(max(abs(ef)));
el = el./max(max(abs(el)));
rf = rf./max(max(abs(rf)));
rl = rl./max(max(abs(rl)));

% Correct the phase slopes
ef = ef.*exp(1i.*0.*qx.'./delta_qx).*exp(-1i.*0.08.*qy./delta_qy);
rf = rf.*exp(-1i.*0.1115.*qx.'./delta_qx).*exp(-1i.*0.58.*qy./delta_qy);
rl = rl.*exp(-1i.*0.1115.*qx.'./delta_qx).*exp(1i.*0.035.*qy./delta_qy).*exp(1i.*0.045);


%% Make the plots
% Make plot matrices
cr = 4;
pef = log10(abs(ef))./cr + 1;
pel = log10(abs(el))./cr + 1;
prf = log10(abs(rf))./cr + 1;
prl = log10(abs(rl))./cr + 1;
pef(pef < 0) = 0;
pel(pel < 0) = 0;
prf(prf < 0) = 0;
prl(prl < 0) = 0;

% Make color-coded arrays
cmap = hsv(256);
pef = complex2rgb(pef.*exp(1i.*angle(ef)),cmap);
pel = complex2rgb(pel.*exp(1i.*angle(el)),cmap);
prf = complex2rgb(prf.*exp(1i.*angle(rf)),cmap);
prl = complex2rgb(prl.*exp(1i.*angle(rl)),cmap);

% Set plotting font size
fnt = 51;
s = 1e-9;
u = '[nm^{-1}]';
xl = [-0.2 0.2];
yl = [-0.2 0.2];
xt = -0.2:0.2:0.2;
yt = -0.2:0.2:0.2;

% Plot free space simulation
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,pef);
axis equal tight;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_ef.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_ef.png'],'-dpng','-r600');
end

% Plot lens-based simulation
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,pel);
axis equal tight;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_el.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_el.png'],'-dpng','-r600');
end

% Plot free space reconstruction
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,prf);
axis equal tight;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_rf.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_rf.png'],'-dpng','-r600');
end

% Plot lens-based reconstruction
figure;
set(gcf,'Position',[460 100 1000 1000]);
image(s*qx,s*qy,prl);
axis equal tight;
xlabel(['q_2 ' u]);
ylabel(['q_1 ' u]);
set(gca,'FontSize',fnt,'YDir','normal','XLim',xl,'XTick',xt,'XTickLabel',num2str(xt.','%.1f'),'YLim',yl,'YTick',yt,'YTickLabel',num2str(yt.','%.1f'));
if sav == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    %print([ps 'Diffraction_rl.pdf'],'-dpdf','-r0');
    print([ps 'Diffraction_rl.png'],'-dpng','-r600');
end



