function [delta_q1,delta_q2,delta_q3,delta_r1,delta_r2,delta_r3] = q_range2(D,det_dx,th,lambda,delta_omega,nx,ny,n_omega)
% Calculate q-space and real space grid

% Angle spanned by a single pixel (in radians)
delta_angle_ccd = det_dx/D;

% Calculate the wave vector
k = 2*pi/lambda;

% Calculate the q-vectors in each measurement direction
q1 = k.*[-sin(th);0;cos(th)].*delta_angle_ccd;
q2 = k.*[0;1;0].*delta_angle_ccd;
q3 = k.*[2*sin(th);0;0].*delta_omega;

% Calculate the vector lengths
delta_q1 = norm(q1);
delta_q2 = norm(q2);
delta_q3 = norm(q3);

% The real space in measurement space
V = cross(ny*q1,nx*q2).'*(n_omega*q3);
r1 = ((2*pi)/V)*cross(nx*q2,n_omega*q3);
r2 = ((2*pi)/V)*cross(n_omega*q3,ny*q1);
r3 = ((2*pi)/V)*cross(ny*q1,nx*q2);
delta_r1 = norm(r1);
delta_r2 = norm(r2);
delta_r3 = norm(r3);

