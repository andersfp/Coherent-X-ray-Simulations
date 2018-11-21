% Initialization
clear;
close all;
clc;


%% Cameraman example
% Load cameraman
A = double(imread('cameraman.tif'));

% Shift image
B = imshift2(A,[0.5 0]);
B = abs(B);

% Plot the images on linear scale
figure;
subplot(1,2,1);
imagesc(A,[0 255]);
axis equal tight;
subplot(1,2,2);
imagesc(B,[0 255]);
axis equal tight;

% Plot the images on log scale
figure;
subplot(1,2,1);
imagesc(log10(A),[0 2.4]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[0 2.4]);
axis equal tight;


%% Square example
% Generate square
A = recta((-128:127)./32).*recta((-128:127).'./32);

% Shift image
B = imshift2(A,[0.5 0]);
B = abs(B);

% Plot the images on linear scale
figure;
subplot(1,2,1);
imagesc(A,[0 1]);
axis equal tight;
subplot(1,2,2);
imagesc(B,[0 1]);
axis equal tight;

% Plot the images on log scale
figure;
subplot(1,2,1);
imagesc(log10(A),[-4 0]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[-4 0]);
axis equal tight;


%% Data example
% Load data
A = load('Processed_Virtual_Shifted.mat');
A = double(A.I(:,:,139,1));

% Shift image
B = imshift2(A,[0.5 0]);
B = abs(B);

% Plot the images on linear scale
figure;
subplot(1,2,1);
imagesc(A,[0 1e6]);
axis equal tight;
subplot(1,2,2);
imagesc(B,[0 1e6]);
axis equal tight;

% Plot the images on log scale
figure;
subplot(1,2,1);
imagesc(log10(A),[0 6]);
axis equal tight;
subplot(1,2,2);
imagesc(log10(B),[0 6]);
axis equal tight;


