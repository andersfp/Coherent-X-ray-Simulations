% Initialization
clear;
close all;
clc;


%% Set parameters
% Load experimental parameters
load('Sim_Parameters.mat');

% Sampling parameters
m = 1000;
dx = 10e-6;

% Get aperture
y0 = sigA*d0;
a0 = sigA;
yN = y0*(cos(N*phi) + phi/2*sin(N*phi)) + a0*(f*phi*sin(N*phi) - T/2*cos(N*phi));


%% Make the object
% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
E0 = zeros(m,m);
E0(x0==0,x0==0) = 1;

% Set the energy RMS bandwidth
w = 1e-2;
dE = w*E;

% Make the attenuation function
att = @(x) sqrt(exp(-(x.^2 + x.'.^2)/(2*yN.^2)));

% Make the propagation distance and focal length arrays
D = [d0+T/2;repmat(T,N-1,1);dd+T/2];
F = repmat(f,N,1);

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);


%% Propagation kernel
nn = 17; % delta_sig < 1e-4, w = 1e-2 -> nn = 17 -> dEs = 30.0 eV, w = 1e-3 -> nn = 5 -> dEs = 10.2, w = 1e-4 -> nn = 2 -> dEs = 2.55
h = waitbar(0,'Progress');
A = zeros(length(nn),1);
sig = A;
for j = 1:length(nn)
    waitbar(j/length(nn),h);
    % Set the half number of wavelength samples (odd)
    n = nn(j);
    
    % Generate energy and wavelength sample positions
    Es = (-n:n).'/n*3*dE + E;
    lambda = 1e-10*12398.42./Es;
    
    % Generate the distribution weight
    S = exp(-(Es - E).^2/(2*dE^2));
    
    % Calculate the total number of points
    n = 2*n + 1;
    
    % Propagate to the detector
    Ed = zeros(m,m,n);
    xd = zeros(m,n);
    %h = waitbar(0,'Progress');
    tic;
    for i = 1:n
        %    waitbar(i/n,h);
        [a,Rm,Rp,sm,sp,gm,gp] = Lens_Stack(D,F,lambda(i),R0,s0);
        [Et,~,xt] = propFrFFT(E0,x0,x0,Inf,Inf,sm(1),sp(N),sum(a(1:N)),lambda(i),'gpu');
        [Ed(:,:,i),~,xd(:,i)] = propFrFFT(Et.*att(xt),xt,xt,Inf,Inf,sm(N+1),sp(N+1),a(N+1),lambda(i),'gpu');
    end
    toc;
    %close(h);
    
    % Make the intensity
    Id = abs(Ed).^2;
    
    % Get the weighted average intensity
    I = sum(Id.*permute(S,[2 3 1]),3)/sum(S);
    
    % Get the central distance axis
    x = xd(:,(n-1)/2);
    
    % Plot the weighted intensity
%     figure;
%     imagesc(x,x,I);
%     axis equal tight;
    
    % Fit the PSF
    Ld = squeeze(Id(:,x == 0,:));
    L = Ld(:,(n-1)/2);
    fun = @(A,sig,x) A.*exp(-x.^2./(2*sig.^2));
    % A = zeros(n,1);
    % sig = A;
    % for i = 1:n
    %     ff = fit(x,Ld(:,i),fun,'StartPoint',[max(Ld(:,i)) 2.0e-7]);
    %     A(i) = ff.A;
    %     sig(i) = ff.sig;
    % %     figure;
    % %     plot(ff,x,Ld(:,i));
    % %     set(gca,'XLim',[-1e-6 1e-6]);
    % end
    fd = fit(x,L,fun,'StartPoint',[max(L) 2.0e-7]);
    A(j) = fd.A;
    sig(j) = fd.sig;
    % figure;
    % plot(fd,x,L);
    % set(gca,'XLim',[-1e-6 1e-6]);
    
    % figure;
    % plot(Es*1e-3,A,[Es(1) Es(n)]*1e-3,fd.A*[1 1]);
    % xlabel('Energy (keV)');
    % ylabel('Intensity (a.u.)');
    %
    % figure;
    % plot(Es*1e-3,1e6*sig,[Es(1) Es(n)]*1e-3,1e6*fd.sig*[1 1]);
    % xlabel('Energy (keV)');
    % ylabel('RMS width (\mum)');
    
end
close(h);

figure;
plot(nn*2+1,A);
xlabel('Energy samples');
ylabel('Intensity (a.u.)');

figure;
plot(nn*2+1,1e6*sig);
xlabel('Energy samples');
ylabel('RMS width (\mum)');

figure;
plot(diff(sig)./sig(1:end-1));
xlabel('Energy samples');
ylabel('RMS width (\mum)');

