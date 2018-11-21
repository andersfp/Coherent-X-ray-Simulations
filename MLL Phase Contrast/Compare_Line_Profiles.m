% Initialization
clear;
close all;
clc;


%% Load the data
% Load the simulated intensity
load('Simulated_Line_Profiles.mat');

% Load the measured intensity
load('Measured_Line_Profiles.mat');

% Normalize the measured intensity
nn = 10;
k = mean(l2([1:nn end-nn:end],:));
l2 = l2./k;
img2 = img2./permute(k,[3 1 2]);

% Invert the measurement
l2 = 2 - l2;
img2 = 2 - img2;


%% Plot the data
% Plot parameters
ylim = [0.70 1.30];

% Plot the zero-position image
figure;
imagesc(img2(:,:,5),[0.6 1.4]);
colormap gray;
axis equal tight;
title('Raw image at focus position');

% Plot the simulated profiles
figure;
plot(1e3*xl,l(:,5:end));
set(gca,'YLim',ylim);
title('Positively shifted line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('0 mm','+1 mm','+2 mm','+3 mm','+4 mm');

% Plot the line profiles of 1 mm pair
figure;
plot(1e3*xl,l(:,[4 6]));
set(gca,'YLim',ylim);
title('Pair of 1 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-1 mm','+1 mm');

% Plot the line profiles of 2 mm pair
figure;
plot(1e3*xl,l(:,[3 7]));
set(gca,'YLim',ylim);
title('Pair of 2 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-2 mm','+2 mm');

% Plot the line profiles of 3 mm pair
figure;
plot(1e3*xl,l(:,[2 8]));
set(gca,'YLim',ylim);
title('Pair of 3 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-3 mm','+3 mm');

% Plot the line profiles of 4 mm pair
figure;
plot(1e3*xl,l(:,[1 9]));
set(gca,'YLim',ylim);
title('Pair of 4 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-4 mm','+4 mm');

% Plot the measured profiles
figure;
plot(1e3*xl2,l2);
set(gca,'YLim',ylim);
title('Selected line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-4 mm','-3 mm','-2 mm','-1 mm','0 mm','+1 mm','+2 mm');

% Plot the 1 mm pair
figure;
plot(1e3*xl2,l2(:,[4 5 6]));
set(gca,'YLim',ylim);
title('Pair of 1 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-1 mm','0 mm','+1 mm');

% Plot the 2 mm pair
figure;
plot(1e3*xl2,l2(:,[3 5 7]));
set(gca,'YLim',ylim);
title('Pair of 2 mm line profiles');
xlabel('Distance (mm)');
ylabel('Intensity (a.u.)');
legend('-2 mm','0 mm','+2 mm');

% Plot the corresponding measured and simulated profiles
n = size(l2,2);
for i = 1:n
    figure;
    plot(1e3*xl,l(:,i),1e3*xl2,l2(:,i));
    set(gca,'YLim',ylim);
    title(['Displacement: ' num2str(d2(i),'%+i') ' mm']);
    xlabel('Distance (mm)');
    ylabel('Normalized intensity');
    legend('Simulation','Measurement');
end

% Compare the raw images
for i = 1:n
    figure;
    subplot(1,2,1);
    imagesc(img2(:,:,i),ylim);
    axis equal tight;
    title(['Measurement: ' num2str(d2(i),'%+i') ' mm']);
    subplot(1,2,2);
    imagesc(Ic(:,:,i),ylim);
    axis equal tight;
    set(gca,'YDir','normal');
    title(['Simulation: ' num2str(1e3*d(i),'%+i') ' mm']);
end

