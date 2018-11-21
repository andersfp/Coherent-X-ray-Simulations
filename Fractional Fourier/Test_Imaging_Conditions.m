% Initialization
clear;
close all;
clc;


%% Set parameters
% Load the parameters
load('Sim_Parameters.mat');


%% Calculate parameters
% Set the side lengths
Dx = logspace(-9,0,10);
%M = logspace(3,12,10);
M = [10000 1000 500 200 100];

% Set simulation ranges
rng = [1e-3 1e-3 1e-3 1e-3 1e-2 4 4 4 4 4];
rng1 = [-4.466e-4 -4.6e-4 -5e-4 -5e-3 -8e-3 -4.0 -3.5 -3.27 -3.26662 -3.266585];
rng2 = [-4.465e-4 -4.35e-4 -4e-4 4e-3 1e-3 1.0 -3.0 -3.264 -3.26654 -3.266583];

% Make plots
k = length(M);
for j = 1:k
    % Set up object parameters
    m = M(j); % M(j), 1000
    dx = 30e-6; % Dx(j), lambda*sqrt(m)*10^(j - 1), 1000*lambda*sqrt(m)
    
    % Make the propagation distance and focal length arrays
    D = [d0 + T/2;T*ones(N-1,1);T/2 + dd];
    F = f*ones(N,1);
    
    % Calculate the initial propagation parameters
    R0 = Inf;
    s0 = dx/sqrt(m);
    
    % Optimize the d2 distance
    % D = optimize_d2(D,F,lambda,R0,s0);
    % d2c = D(end) - (T/2 + dd);
    
    % Set up the displacement array
    %d2 = linspace(rng1(j),rng2(j),1000);
    d2 = linspace(-dd,dd,1000);
    n = length(d2);
    A = zeros(n,1);
    R = A;
    
    % Calculate the FrFT parameters
    for i = 1:n
        [a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D + [zeros(N,1);d2(i)],F,lambda,R0,s0);
        A(i) = sum(a);
        R(i) = Rp(end);
    end
    
    % Plot the critical parameters
    figure;
    ax = plotyy(d2,A,d2,R);
    set(ax(1),'YLim',[-0.25 3.25],'YTick',0:0.5:3);
    set(ax(2),'YLim',[-7 7],'YTick',-6:2:6);
end


