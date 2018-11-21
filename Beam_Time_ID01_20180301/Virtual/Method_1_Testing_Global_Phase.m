% Initialization
clear;
close all;
clc;


%% Load data
% Load the amplitude
A = load('Amplitude_Ptycho_1_Individual_Reconstructions.mat');

% Load the phase
P = load('Phase_Ptycho_1_Individual_Reconstructions.mat');

% Generate the full field
fieldi = A.A.*exp(1i.*P.A);


%% Average the fields
% Do direct averaging
F1 = mean(fieldi,4);

% Generate the object
R1 = fftshift(ifftn(fftshift(F1)));

% Get the average phase
w = 3;
p = mean(mean(mean(angle(fieldi(129-w:129+w,129-w:129+w,141-w:141+w,:)))));

% Correct the global phase
fieldi2 = fieldi.*exp(-1i.*p);

% Do the averaging
F2 = mean(fieldi2,4);

% Generate the object
R2 = fftshift(ifftn(fftshift(F2)));

% Calculate the differences
DF = (abs(F2) - abs(F1))./abs(F1);
DR = (abs(R2) - abs(R1))./abs(R1);



