% Initialization
clear;
close all;
clc;

% Set saving path
p = 'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Lens_Aberrations_Dislocations_Phase\';

% Save the phantom?
sav = 1;


%% Generate the shape
% Set the size of the phantom volume
m = 256;
A = ones(m,m,m);

% Generate coordinates
x = (-m/2):(m/2-1);
y = x.';
z = permute(x,[3 1 2]);

% Set sizes of the polyhedrons and randomness
rc = m/2*0.7; % 0.8
ro = m/2*0.7; % 0.8
rnd = 0;

% Set the cubic plane positions
c1 = ([1 0 0] + (rand(1,3) - 0.5)*rnd)*rc;
c2 = ([-1 0 0] + (rand(1,3) - 0.5)*rnd)*rc;
c3 = ([0 1 0] + (rand(1,3) - 0.5)*rnd)*rc;
c4 = ([0 -1 0] + (rand(1,3) - 0.5)*rnd)*rc;
c5 = ([0 0 1] + (rand(1,3) - 0.5)*rnd)*rc;
c6 = ([0 0 -1] + (rand(1,3) - 0.5)*rnd)*rc;

% Calculate the cubic normal vectors
c1n = c1./sqrt(dot(c1,c1));
c2n = c2./sqrt(dot(c2,c2));
c3n = c3./sqrt(dot(c3,c3));
c4n = c4./sqrt(dot(c4,c4));
c5n = c5./sqrt(dot(c5,c5));
c6n = c6./sqrt(dot(c6,c6));

% Set the octahedral plane positions
o1 = ([1 1 1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o2 = ([1 1 -1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o3 = ([1 -1 1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o4 = ([-1 1 1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o5 = ([1 -1 -1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o6 = ([-1 1 -1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o7 = ([-1 -1 1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;
o8 = ([-1 -1 -1]/sqrt(3) + (rand(1,3) - 0.5)*rnd)*ro;

% Calculate the octahedral normal vectors
o1n = o1./sqrt(dot(o1,o1));
o2n = o2./sqrt(dot(o2,o2));
o3n = o3./sqrt(dot(o3,o3));
o4n = o4./sqrt(dot(o4,o4));
o5n = o5./sqrt(dot(o5,o5));
o6n = o6./sqrt(dot(o6,o6));
o7n = o7./sqrt(dot(o7,o7));
o8n = o8./sqrt(dot(o8,o8));

% Apply the planes
A(c1n(1)*(x - c1(1)) + c1n(2)*(y - c1(2)) + c1n(3)*(z - c1(3)) > 0) = 0;
A(c2n(1)*(x - c2(1)) + c2n(2)*(y - c2(2)) + c2n(3)*(z - c2(3)) > 0) = 0;
A(c3n(1)*(x - c3(1)) + c3n(2)*(y - c3(2)) + c3n(3)*(z - c3(3)) > 0) = 0;
A(c4n(1)*(x - c4(1)) + c4n(2)*(y - c4(2)) + c4n(3)*(z - c4(3)) > 0) = 0;
A(c5n(1)*(x - c5(1)) + c5n(2)*(y - c5(2)) + c5n(3)*(z - c5(3)) > 0) = 0;
A(c6n(1)*(x - c6(1)) + c6n(2)*(y - c6(2)) + c6n(3)*(z - c6(3)) > 0) = 0;
A(o1n(1)*(x - o1(1)) + o1n(2)*(y - o1(2)) + o1n(3)*(z - o1(3)) > 0) = 0;
A(o2n(1)*(x - o2(1)) + o2n(2)*(y - o2(2)) + o2n(3)*(z - o2(3)) > 0) = 0;
A(o3n(1)*(x - o3(1)) + o3n(2)*(y - o3(2)) + o3n(3)*(z - o3(3)) > 0) = 0;
A(o4n(1)*(x - o4(1)) + o4n(2)*(y - o4(2)) + o4n(3)*(z - o4(3)) > 0) = 0;
A(o5n(1)*(x - o5(1)) + o5n(2)*(y - o5(2)) + o5n(3)*(z - o5(3)) > 0) = 0;
A(o6n(1)*(x - o6(1)) + o6n(2)*(y - o6(2)) + o6n(3)*(z - o6(3)) > 0) = 0;
A(o7n(1)*(x - o7(1)) + o7n(2)*(y - o7(2)) + o7n(3)*(z - o7(3)) > 0) = 0;
A(o8n(1)*(x - o8(1)) + o8n(2)*(y - o8(2)) + o8n(3)*(z - o8(3)) > 0) = 0;


%% Generate the phase
% Set the real space axes
fov = 5e-6;
rx = x./m.*fov;
ry = y./m.*fov;
rz = z./m.*fov;

% Set positions
d = 200e-9;
py = linspace(-fov/2,fov/2,round(fov./d)).';
d = mean(diff(py));

% Remove positions outside the object
As = sum(A,3);
py(py > max(ry.*(As(:,m/2 + 1) > 0))) = [];
py(py < min(ry.*(As(:,m/2 + 1) > 0))) = [];

% Set Al(110) material parameters
nu = 0.348;
b = 2.86e-10;
a = 4.04e-10;
G = sqrt(2).*2.*pi./a;
v = b./d.*180./pi;
s = 0.5.*pi.*b./(1 - nu)./d./d.*rx./sinh(pi.*rx./d)./sinh(pi.*rx./d);

% Calculate the displacement field (Al(110) edge dislocation wall)
u = zeros(m,m,1);
n = length(py);
for i = 1:n
    u = u + b./(8.*pi.*(1 - nu)).*((1 - 2.*nu).*log(rx.^2 + (ry - py(i)).^2) + (rx.^2 - (ry - py(i)).^2)./(rx.^2 + (ry - py(i)).^2));
end

% Remove NaNs
ii = find(isnan(u));
[iy,ix] = ind2sub([m m],ii);
u(ii) = (u(iy + 1,m/2 + 1) + u(iy - 1,m/2 + 1) + u(iy,m/2 + 2) + u(iy,m/2))./4;

% Plot the displacement field
figure;
imagesc(rx,ry,u.*(As > 0));
axis equal tight;
colorbar;

% Calculate phase
P = exp(1i.*u.*G);
figure;
imagesc(rx,ry,angle(P).*(As > 0));
axis equal tight;
colormap hsv;

% Combine into final phantom
obj = A.*repmat(P,1,1,m);

% Generate the Bragg vector
G = [3;3;3].*2.*pi./a;


%% Save the phantom
% Save phantom if specified
if sav == 1
    save([p 'Phantom.mat'],'obj','rx','ry','rz','a','b','d','G','nu','s','v','fov','x','y','z');
end


