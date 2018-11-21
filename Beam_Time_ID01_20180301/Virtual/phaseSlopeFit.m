function [bx,by,bz] = phaseSlopeFit(A)
% Fit the slope of the phase in the x-, y-, and z-directions of the
% different sets of A along the 4th dimension, normalized to the first of
% these sets.
%
% Example of usage:
% [bx,by,bz] = phaseSlopeFit(A)
%

% Get the size of the input
nx = size(A,2);
ny = size(A,1);
nz = size(A,3);
n = size(A,4);

% Generate pixel axes
x = ((-nx/2):(nx/2 - 1)).';
y = ((-ny/2):(ny/2 - 1)).';
z = ((-nz/2):(nz/2 - 1)).';

% Extract lines
ax = squeeze(A(ny/2 + 1,:,nz/2 + 1,:));
ay = squeeze(A(:,nx/2 + 1,nz/2 + 1,:));
az = squeeze(A(ny/2 + 1,nx/2 + 1,:,:));

% Convert lines to double on the CPU
ax = double(gather(ax));
ay = double(gather(ay));
az = double(gather(az));

% Calculate the phase
px = unwrap(angle(ax./ax(:,1)));
py = unwrap(angle(ay./ay(:,1)));
pz = unwrap(angle(az./az(:,1)));

% Calculate the exclusions
e = 0.25;
ex = abs(ax) < e.*max(abs(ax));
ey = abs(ay) < e.*max(abs(ay));
ez = abs(az) < e.*max(abs(az));

% Fit the slopes
bx = zeros(n,1);
by = bx;
bz = bx;
for i = 2:n
    fx = fit(x,px(:,i),'poly1','Exclude',ex(:,i));
    fy = fit(y,py(:,i),'poly1','Exclude',ey(:,i));
    fz = fit(z,pz(:,i),'poly1','Exclude',ez(:,i));
    bx(i) = -fx.p1.*nx./(2.*pi);
    by(i) = -fy.p1.*ny./(2.*pi);
    bz(i) = -fz.p1.*nz./(2.*pi);
end

% Reshape the output
bx = permute(bx,[2 3 4 1]);
by = permute(by,[2 3 4 1]);
bz = permute(bz,[2 3 4 1]);

