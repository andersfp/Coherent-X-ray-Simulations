% Initialization
clear;
close all;
clc;


%% Load the data
% Set the data path
dp = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\';

% Load the simulation
tic;
load([dp 'Star_Simulation_Series_1_Partial_4_1000.mat']);
toc;

% Load additional simulations
fs = {'4_1000','5_200','6_200','7_200','8_800','9_200'};
n = length(fs);
S0 = zeros([size(I) n]);
for i = 1:n
    B = load([dp 'Star_Simulation_Series_1_Partial_' fs{i} '.mat']);
    S0(:,:,:,i) = B.I;
    fprintf('.');
end
fprintf('\n');

% Clear unused data
clear B;

% Average the separate simulations
w = [1000 200 200 200 800 200];
I = sum(S0.*permute(w,[3 4 1 2]),4)./sum(w);

% Plot the simulation
Slicer(I,'displayRange',[16000 52000]);

% Save the image stack
save([dp 'Star_Simulation_Series_1_Partial_Avg.mat'],'I','x','y','d','lx','ly','-v7.3');


