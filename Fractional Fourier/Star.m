% Initialization
clear;
close all;
clc;

% Star
n=5; % number of star rays
square_side_pixels = 500; % square size in pixels
star_fatness = 0.5; % just try to change it a little

t = linspace(pi/2,2.5*pi,n*2+1)+0.1; % angles of all points
r=[repmat([1 star_fatness],1,n) 1]; % distance of all points from the center
[x,y] = pol2cart(t,r); % convert polar coordinates to cartesian
fill(x,y,'y') % draw the star
axis square
set(gca,'Color','b','xtick',[],'ytick',[])
set(gcf,'units','pixels');
pos = get(gcf,'position');
pos(3:4) = repmat(square_side_pixels,1,2);
set(gcf,'position',pos); % make figure square
set(gca,'Position',[0 0 1 1])
A = getframe(gcf); % save figure as image matrix
A = A.cdata;

S = double(A(:,:,1) > 200);

if square_side_pixels < 120
    m = (120 - square_side_pixels)/2;
    S = S(:,m+1:m+square_side_pixels);
end

figure;
imagesc(S);
axis equal tight;

save('Star.mat','S');

