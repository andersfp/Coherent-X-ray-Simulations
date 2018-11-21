function [delta_q1,delta_q2,delta_q3,delta_r1_prime,delta_r2_prime,delta_r3_prime] = q_range(D,det_dx,tth,lambda,omega_cen,delta_omega,nx,ny,n_omega)
% Calculate q-space and real space grid

% Incidence angles
omega = delta_omega*((-n_omega/2):(n_omega/2 - 1)) + omega_cen;
omega = permute(omega,[1 3 2]);

% Angle spanned by a single pixel (in radians)
delta_angle_ccd = det_dx/D;

% Generate a grid on the ccd
[x,y] = meshgrid((-nx/2):(nx/2-1),(-ny/2):(ny/2-1));

% Set the angles of each pixel at each omega angle
alphai = omega;
%alphaf = y.*delta_angle_ccd + tth - omega;
alphaf = bsxfun(@minus,y.*delta_angle_ccd + tth,omega);
nu = x.*delta_angle_ccd;

% Cartesian coordinates of the scattering vector
%qx = 2*pi/lambda*(cos(alphaf).*cos(nu) - cos(alphai));
qx = 2*pi/lambda*bsxfun(@minus,bsxfun(@times,cos(alphaf),cos(nu)),cos(alphai));
qy = 2*pi/lambda*sin(nu);
qy = repmat(qy,1,1,n_omega);
%qz = 2*pi/lambda*(sin(alphaf).*cos(nu) + sin(alphai));
qz = 2*pi/lambda*bsxfun(@plus,bsxfun(@times,sin(alphaf),cos(nu)),sin(alphai));

% Voxel size of q-space in measurement coordinates
delta_q1 = norm([qx(ny/2 + 1,nx/2,n_omega/2) - qx(ny/2,nx/2,n_omega/2),qy(ny/2 + 1,nx/2,n_omega/2) - qy(ny/2,nx/2,n_omega/2),qz(ny/2 + 1,nx/2,n_omega/2) - qz(ny/2,nx/2,n_omega/2)]);
delta_q2 = norm([qx(ny/2,nx/2 + 1,n_omega/2) - qx(ny/2,nx/2,n_omega/2),qy(ny/2,nx/2 + 1,n_omega/2) - qy(ny/2,nx/2,n_omega/2),qz(ny/2,nx/2 + 1,n_omega/2) - qz(ny/2,nx/2,n_omega/2)]);
delta_q3 = norm([qx(ny/2,nx/2,n_omega/2 + 1) - qx(ny/2,nx/2,n_omega/2),qy(ny/2,nx/2,n_omega/2 + 1) - qy(ny/2,nx/2,n_omega/2),qz(ny/2,nx/2,n_omega/2 + 1) - qz(ny/2,nx/2,n_omega/2)]);

% Sample space basis
uz = [1 0 0];
uy = [0 1 0];
ux = [0 0 1];

% Generate scattering vector coordinates in measurement space
alphaf_cen = tth - omega_cen ;
q1 = delta_q1*(cos(pi/2 + alphaf_cen)*ux + cos(alphaf_cen)*uz);
q2 = delta_q2*uy;
q3_prime_inter = (sin(alphaf_cen) + sin(omega_cen))*ux + (-cos(alphaf_cen) + cos(omega_cen))*uz;
q3_prime = delta_q3*q3_prime_inter/norm(q3_prime_inter);

% The real space in measurement space
V = cross(ny*q1,nx*q2)*n_omega*q3_prime';
r1_prime = ((2*pi)/V)*cross(nx*q2,n_omega*q3_prime);
r2_prime = ((2*pi)/V)*cross(n_omega*q3_prime,ny*q1);
r3_prime = ((2*pi)/V)*cross(ny*q1,nx*q2);
delta_r1_prime = norm(r1_prime);
delta_r2_prime = norm(r2_prime);
delta_r3_prime = norm(r3_prime);



