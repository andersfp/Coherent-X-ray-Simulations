% Initialization
clear;
close all;
clc;


%% Mesh
% Generate the basic grid
m = 10;
[x,y] = meshgrid(linspace(0,2,2*m),linspace(0,1,m));

% Add randomness to the grid
a = 0.6/(m-1);
x(2:end-1,2:end-1) = x(2:end-1,2:end-1) + a*rand(m-2,2*m-2) - a/2;
y(2:end-1,2:end-1) = y(2:end-1,2:end-1) + a*rand(m-2,2*m-2) - a/2;

% Convert coordinate grid to column vectors
x = x(:);
y = y(:);

% Perform the triangulation
DT = delaunayTriangulation(x,y);

% Extract the points and connectivity from the triangulation
p = DT.Points;
cl = DT.ConnectivityList;

% Plot all the connections
n = size(cl,1);
ls = 'k';
lw = 2;
figure;
for i = 1:n
    plot([p(cl(i,1),1) p(cl(i,2),1)],[p(cl(i,1),2) p(cl(i,2),2)],ls,'LineWidth',lw);
    hold on;
    plot([p(cl(i,2),1) p(cl(i,3),1)],[p(cl(i,2),2) p(cl(i,3),2)],ls,'LineWidth',lw);
    plot([p(cl(i,3),1) p(cl(i,1),1)],[p(cl(i,3),2) p(cl(i,1),2)],ls,'LineWidth',lw);
end
plot([0 0 2 2 0 0],[0 1 1 0 0 1],ls,'LineWidth',2*lw);
axis equal;
set(gca,'XLim',[0 2],'YLim',[0 1]);
axis off;
set(gcf,'Position',[1 41 1920 963]);

img = print('-RGBImage');
img = rgb2gray(img);

img(:,mean(img,1) == 255) = [];
img(mean(img,2) == 255,:) = [];

% Save the mesh figure
imwrite(img,'Mesh.png');


