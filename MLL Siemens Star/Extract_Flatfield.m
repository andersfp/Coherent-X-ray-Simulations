% Initialization
clear;
close all;
clc;


%% Define the file
% Set the data path
p1 = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_2D_17.3keV\pco_flatfield_deco_in\';
p2 = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_2D_17.3keV\pco_flatfield_deco_out\';
p3 = 'C:\Users\anfils\Documents\Simulation_Results\Siemens_Star\Data\Set4_2D_17.3keV\pco_flatfield_tf_decoh\';

% Set the base name
fp = '*0*.edf';

% Make the file lists
[~,f1] = fileNames(p1,fp);
[~,f2] = fileNames(p2,fp);
[~,f3] = fileNames(p3,fp);
f = cat(1,f1,f2,f3);


%% Load the data
% Get number of images
n = length(f);

% Pre-allocate info structure
info(n) = edf_info(f{n});

% Set the image dimensions
mx = info(n).dim_1;
my = info(n).dim_2;

% Read the images
img = zeros(my,mx,n);
fprintf([num2str(n) ': ']);
for i = 1:n
    img(:,:,i) = edf_read(f{i});
    info(i) = edf_info(f{i});
    fprintf('.');
    if rem(i,100) == 0
        fprintf('\n');
    end
end
fprintf('\n');

% Plot the flatfield images
Slicer(img,'displayRange',[100 1500]);


%% Create a mean flatfield
% Median flatfield (2 versions)
ff1 = median(img(:,:,7:11),3);
ff2 = median(img(:,:,1:6),3);

% Plot the flatfields
figure;
imagesc(ff1,[100 1500]);
axis equal tight;
colormap gray;
figure;
imagesc(ff2,[100 1500]);
axis equal tight;
colormap gray;

% Save the flatfields
save('Flatfield.mat','ff1','ff2','-v7.3');


