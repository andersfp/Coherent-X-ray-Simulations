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

% Generate the plots
obj1 = 1 - complex2rgb(obj(:,:,i1),cmap);
obj2 = 1 - complex2rgb(obj(:,:,i2),cmap);
obj3 = 1 - complex2rgb(squeeze(obj(:,i1,:)),cmap);
obj4 = 1 - complex2rgb(squeeze(obj(:,i2,:)),cmap);
rec1 = 1 - complex2rgb(rec(:,:,i1),cmap);
rec2 = 1 - complex2rgb(rec(:,:,i2),cmap);
rec3 = 1 - complex2rgb(squeeze(rec(:,i1,:)),cmap);
rec4 = 1 - complex2rgb(squeeze(rec(:,i2,:)),cmap);
na41 = 1 - complex2rgb(na4(:,:,i1),cmap);
na42 = 1 - complex2rgb(na4(:,:,i2),cmap);
na43 = 1 - complex2rgb(squeeze(na4(:,i1,:)),cmap);
na44 = 1 - complex2rgb(squeeze(na4(:,i2,:)),cmap);
par1 = 1 - complex2rgb(par(:,:,i1),cmap);
par2 = 1 - complex2rgb(par(:,:,i2),cmap);
par3 = 1 - complex2rgb(squeeze(par(:,i1,:)),cmap);
par4 = 1 - complex2rgb(squeeze(par(:,i2,:)),cmap);

% Concatenate the plots
rgb = cat(5,cat(4,par1,par2,par3,par4),cat(4,na41,na42,na43,na44),cat(4,rec1,rec2,rec3,rec4),cat(4,obj1,obj2,obj3,obj4));

% Make a color scale
xc = -64:63;
yc = xc.';
rc = sqrt(xc.^2 + yc.^2)./60;
ac = atan2(yc,xc);
cs = rc.*exp(1i.*ac);
cs = cs.*(rc <= 1);
c = 1 - complex2rgb(cs,cmap);

% Make a plot
lt = {{'Reconstruction','partial','coherence'},{'Reconstruction','4\timesNA'},'Reconstruction',{'Ground','truth'}};
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
figure('Position',[100 50 oh1+nh*w+oh2 ov1+nv*h+ov2]);
axes('Units','pixels','Position',[(oh1-128)/2 ov1+ov2+nv*h-128-(oh1-128)/2+15 128 128]);
image(c);
axis equal tight off;
text(128,64,'0','FontSize',fnt);
text(64,0-8,'\pi/2','FontSize',fnt,'HorizontalAlignment','center');
text(0,64,'\pi','FontSize',fnt,'HorizontalAlignment','right');
text(64,128+8,'-\pi/2','FontSize',fnt,'HorizontalAlignment','center');
for i = 1:nh
    for j = 1:nv
        axes('Units','pixels','Position',[oh1+(i-1)*w ov1+(j-1)*h w h]);
        image(rgb(:,:,:,i,j));
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
        if j == nv && i == 1
            arrow([20 256-15],[60 256-15],'LineWidth',2);
            arrow([15 256-20],[15 256-60],'LineWidth',2);
            text(75,256-19,'x','FontSize',fnt,'HorizontalAlignment','center');
            text(15,256-85,'y','FontSize',fnt,'HorizontalAlignment','center');
        end
        if j == nv && i == 3
            arrow([20 256-15],[60 256-15],'LineWidth',2);
            arrow([15 256-20],[15 256-60],'LineWidth',2);
            text(75,256-19,'z','FontSize',fnt,'HorizontalAlignment','center');
            text(15,256-85,'y','FontSize',fnt,'HorizontalAlignment','center');
        end
    end
end

% Save the figure
%print(gcf,'Results_Figure.png','-dpng','-r480');


