% Initialization
clear;
close all;
clc;


%% Set a file name
% Set the data path
p = 'C:\Users\anfils\Documents\Simulation_Results\Beam_Time_Data\Pt2_g5_2\';

% File name
lst = dir([p '*.cxi']);
lst = struct2cell(lst);
lst = lst(1,:).';


%% Load the diffraction data
% Set the data path
dp = '/entry_1/instrument_1/detector_1/data';

% Read the data and sum along 3rd direction
n = length(lst);
dat = zeros(516,516,n);
for i = 1:n
    dat(:,:,i) = sum(double(h5read([p lst{i}],dp)),3);
end


%% Process the data
% Sum the array along 3rd direction and transpose
dat = sum(dat,3).';

% Make a list of strange pixels
sp = [206 80;195 50;136 203;121 248;339 73;450 204;267 280;513 398;509 481;337 492;356 272;356 273;354 273;221 428;237 483;210 468;171 384;102 295;67 288;139 363;133 375;130 380;122 374;114 364;55 406;75 419;74 410;61 402;64 401;57 394;58 393;61 392;69 394;72 397;85 393;76 388;63 379;69 377;67 375;83 377];

% Generate indices from the coordinates
ii = sub2ind(size(dat),sp(:,2),sp(:,1));

% Generate the mask
mask = dat > 0;

% Further remove the strange pixels
mask(ii) = 0;

% Plot the mask
figure;
imagesc(mask,[0 1]);
axis equal tight;
title('Mask');


%% Save the mask
% Convert the mask to uint8 to save time and space
mask = uint8(mask);

% Save the mask
save('MaxiPix_Mask.mat','mask');


