% Initialization
clear;
close all;
clc;


%% Load the simulation results
% Load the simulation parameters
load('Sim_Parameters.mat');

% Load the .mat file
tic;
load([p 'Pink_Beam_USAF_980_10um.mat']);
toc;

% Save the result?
s = 1;


%% Process the data
% Generate new axis
% m2 = 1000;
% dx2 = 1e-4;
% x2 = ((-m2/2):(m2/2 - 1)).'/m2*dx2;
m2 = m;
x2 = x(:,(nE + 1)/2);

% Generate new dataset
tic;
I2 = zeros(m2,m2,nE);
for i = 1:nE
    I2(:,:,i) = interp2(x(:,i).',x(:,i),I(:,:,i),x2.',x2,'cubic',0);
end
toc;

% Generate the bandwidth array
%w = [1e-10 1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3 8e-3 9e-3 1e-2];
w = [1e-10 (1:100)*1e-4];
dE = w.*E;
n = length(w);

% Generate weigthing factor
fun = @(sig,x) exp(-x.^2./(2*(sig).^2));
Es = permute(Es,[3 1 2]);
S = fun(dE,Es - E);

% Generate total intensities
% It = zeros(m2,m2,n);
% h = waitbar(0,'Progress');
% tic;
% for i = 1:n
%     for j = 1:nE
%         waitbar(((i - 1)*nE + j)/(n*nE),h);
%         It(:,:,i) = It(:,:,i) + I2(:,:,j).*S(1,i,j);
%     end
% end
% toc;
% close(h);

It = zeros(m2,m2,n);
h = waitbar(0,'Progress');
tic;
for i = 1:n
    waitbar(i/n,h);
    It(:,:,i) = sum(I2.*S(1,i,:),3);
end
toc;
close(h);

% Generate normalized intensity
Itn = It./max(max(It));

% Save the result
if s == 1
    tic;
    save([p 'Pink_Beam_USAF_980_10um_Processed.mat'],'x2','It','Itn','w','-v7.3');
    toc;
end


