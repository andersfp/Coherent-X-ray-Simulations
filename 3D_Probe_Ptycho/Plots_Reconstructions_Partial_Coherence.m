% Initialization
clear;
close all;
clc;


%% Load data
% Set the full file locations
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\';
f = 'Plotting_Data_Partial_Coherence.mat';

% Load the data
load([p f]);


%% Make plots
% Set the plotting indices
i1 = 129;
i2 = 78;

% Generate the colormap
cmap = hsv(256);

% Generate the plots
par1 = 1 - complex2rgb(par(:,:,i1),cmap);
par2 = 1 - complex2rgb(par(:,:,i2),cmap);
par3 = 1 - complex2rgb(squeeze(par(:,i1,:)),cmap);
par4 = 1 - complex2rgb(squeeze(par(:,i2,:)),cmap);

% Concatenate the plots
rgb = cat(4,par1,par2,par3,par4);

% Make a plot
lt = {{'Reconstruction','partial','coherence'}};
tt = {{'xy-plane','z = 0'},{'xy-plane','z = 1 \mum'},{'zy-plane','x = 0'},{'zy-plane','x = 1 \mum'}};
sb = 20;
fnt = 20;
oh1 = 200;
oh2 = 5;
ov1 = 5;
ov2 = 80;
h = 256;
w = 256;
figure('Position',[100 100 oh1+1024+oh2 ov1+256+ov2]);
for i = 1:4
    axes('Units','pixels','Position',[oh1+(i-1)*w ov1 w h]);
    image(rgb(:,:,:,i));
    axis equal tight off;
    rectangle('Position',[0 0 256 256],'LineWidth',3);
    if i == 1
        text(-oh1/2,129,lt{1},'FontSize',fnt,'HorizontalAlignment','center');
    end
    text(129,-ov2/2,tt{i},'FontSize',fnt,'HorizontalAlignment','center');
    if i == 4
        line([256-sb-51.2 256-sb],[256-sb 256-sb],'LineWidth',5,'Color','black');
        text(256-sb-51.2/2,215,'1 \mum','FontSize',18,'HorizontalAlignment','center');
    end
end

% Save the figure
print(gcf,'Results_Figure_Partial_Coherence.png','-dpng','-r480');


