% Initialization
clear;
close all;
clc;


%% Load the simulation results
% Load the simulation parameters
load('Sim_Parameters.mat');

% Load the .mat file
tic;
load([p 'Pink_Beam_USAF_980_10um_Processed.mat']);
toc;

% Save the result?
s = 0;


%% Process the variables
% Rename variables
x = x2;
I = Itn;
w = w.';

% Plot the animation
h = figure;
set(h,'Position',[470 110 980 980]);
for i = 1:length(w)
    imagesc(I(:,:,i));
    axis equal tight off;
    colormap gray;
    set(gca,'Position',[0 0 1 1]);
    text(200,900,['\DeltaE/E = ' num2str(round(1e4*w(i)),'%03u') '\times10^{-4}'],'FontSize',50,'Color','white');
    
    drawnow;
    fr = getframe(h);
    im = frame2im(fr);
    [imind,cm] = rgb2ind(im,256);
    
    if i == 1
        imwrite(imind,cm,'Pink_Beam_USAF.gif','gif','Loopcount',Inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,'Pink_Beam_USAF.gif','gif','DelayTime',0.2,'WriteMode','append');
    end
end



