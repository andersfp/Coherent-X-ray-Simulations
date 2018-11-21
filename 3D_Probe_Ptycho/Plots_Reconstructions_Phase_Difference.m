% Initialization
clear;
close all;
clc;


%% Load data
% Set the full file locations
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\';
f1 = 'Plotting_Data.mat';
f2 = 'Plotting_Data_Partial_Coherence.mat';

% Load the data
load([p f1]);
load([p f2]);


%% Make plots
% Set the plotting indices
i1 = 129;
i2 = 78;

% Generate the colormap
cmap = hsv(256);

% Get the difference phase maps
rec1 = angle(rec(:,:,i1).*conj(obj(:,:,i1)));
rec2 = angle(rec(:,:,i2).*conj(obj(:,:,i2)));
rec3 = squeeze(angle(rec(:,i1,:).*conj(obj(:,i1,:))));
rec4 = squeeze(angle(rec(:,i2,:).*conj(obj(:,i2,:))));
na41 = angle(na4(:,:,i1).*conj(obj(:,:,i1)));
na42 = angle(na4(:,:,i2).*conj(obj(:,:,i2)));
na43 = squeeze(angle(na4(:,i1,:).*conj(obj(:,i1,:))));
na44 = squeeze(angle(na4(:,i2,:).*conj(obj(:,i2,:))));
par1 = angle(par(:,:,i1).*conj(obj(:,:,i1)));
par2 = angle(par(:,:,i2).*conj(obj(:,:,i2)));
par3 = squeeze(angle(par(:,i1,:).*conj(obj(:,i1,:))));
par4 = squeeze(angle(par(:,i2,:).*conj(obj(:,i2,:))));

% Concatenate the plots
rgb = cat(5,cat(4,par1,par2,par3,par4),cat(4,na41,na42,na43,na44),cat(4,rec1,rec2,rec3,rec4));

% Make a plot
cs = 1;
clim = [-pi/cs pi/cs];
lt = {{'Phase diff.','partial','coherence'},{'Phase diff.','4\timesNA'},{'Phase diff.','Reconstruction'}};
tt = {{'xy-plane','z = 0'},{'xy-plane','z = 1 \mum'},{'zy-plane','x = 0'},{'zy-plane','x = 1 \mum'}};
sb = 20;
fnt = 20;
oh1 = 200;
oh2 = 5;
ov1 = 5;
ov2 = 80;
h = 256;
w = 256;
nh = size(rgb,4);
nv = size(rgb,5);
figure('Position',[100 100 oh1+nh*w+oh2 ov1+nv*h+ov2]);
for i = 1:nh
    for j = 1:nv
        axes('Units','pixels','Position',[oh1+(i-1)*w ov1+(j-1)*h w h]);
        imagesc(rgb(:,:,:,i,j),clim);
        axis equal tight off;
        rectangle('Position',[0 0 256 256],'LineWidth',3);
        if i == 1
            text(-oh1/2,129,lt{j},'FontSize',fnt,'HorizontalAlignment','center');
        end
        if j == nv
            text(129,-ov2/2,tt{i},'FontSize',fnt,'HorizontalAlignment','center');
        end
        if i == nh && j == 1
            line([256-sb-51.2 256-sb],[256-sb 256-sb],'LineWidth',5,'Color','black');
            text(256-sb-51.2/2,215,'1 \mum','FontSize',18,'HorizontalAlignment','center');
        end
%         if j == nv && i == 1
%             arrow([20 256-15],[60 256-15],'LineWidth',2);
%             arrow([15 256-20],[15 256-60],'LineWidth',2);
%             text(75,256-19,'x','FontSize',fnt,'HorizontalAlignment','center');
%             text(15,256-85,'y','FontSize',fnt,'HorizontalAlignment','center');
%         end
%         if j == nv && i == 3
%             arrow([20 256-15],[60 256-15],'LineWidth',2);
%             arrow([15 256-20],[15 256-60],'LineWidth',2);
%             text(75,256-19,'z','FontSize',fnt,'HorizontalAlignment','center');
%             text(15,256-85,'y','FontSize',fnt,'HorizontalAlignment','center');
%         end
    end
end
ch = 40;
cv = 130;
%c = colorbar('Units','pixels','Position',[oh1/2 nv*h+ov2-cv-2*ov1 40 cv],'Ticks',[clim(1) 0 clim(2)],'TickLabels',{['-\pi/' num2str(cs)],'0',['\pi/' num2str(cs)]},'FontSize',fnt);
c = colorbar('Units','pixels','Position',[oh1/2 nv*h+ov2-cv-2*ov1 40 cv],'Ticks',[clim(1) 0 clim(2)],'TickLabels',{'-\pi','0','\pi'},'FontSize',fnt);

% Save the figure
%print(gcf,'Results_Phase_Difference_Figure.png','-dpng','-r480');


