% Initialization
clear;
close all;
clc;


%% Load data
% Load experimental parameters
load('Exp_Param.mat');


%% Generate the object
% Generate full coordinates
[X,Y,Z] = meshgrid(x,y,z);

% Make object
Oxyz = double((abs(X) <= 0.5e-6) & (abs(Y) <= 0.5e-6) & (abs(Z) <= 0.5e-6) & (Y <= 0) | (X.^2 + Y.^2 <= (0.5e-6).^2) & (abs(Z) <= 0.5e-6) & (Y >= 0));

% % Make 3D plot of the original object
% figure;
% pp = patch(isosurface(abs(Oxyz),0.1));
% isonormals(abs(Oxyz),pp);
% pp.FaceColor = 'red';
% pp.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight(45,45);
% camlight(-135,-45);
% lighting gouraud;
% xlabel('y');
% ylabel('z');
% zlabel('x');

% % Make 3D slices of the original object
% figure;
% slice(1e6*X,1e6*Y,1e6*Z,abs(Oxyz),0,0,0);
% axis equal;
% shading flat;
% xlabel('y [\mum]');
% ylabel('z [\mum]');
% zlabel('x [\mum]');


%% Convert to measurement space
% Shift the object
O123 = Oxyz;
for i = 1:n_omega
    O123(:,:,i) = circshift(Oxyz(:,:,i),round((i - n_omega/2)*shft),1);
end

% % Make 3D plot of the transformed object
% figure;
% pp = patch(isosurface(abs(O123),0.1));
% isonormals(abs(O123),pp);
% pp.FaceColor = 'red';
% pp.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight(45,45);
% camlight(-135,-45);
% lighting gouraud;
% xlabel('y');
% ylabel('z');
% zlabel('x');

% % Make 3D slices of the transformed object
% figure;
% slice(1e6*X,1e6*Y,1e6*Z,abs(O123),0,0,0);
% axis equal;
% shading flat;
% xlabel('y [\mum]');
% ylabel('z [\mum]');
% zlabel('x [\mum]');


%% Save the object
% Save both original and transformed object
%save([p 'Object.mat'],'Oxyz','O123','-v7.3');
save_binary([p 'Oxyz.bin'],Oxyz);
save_binary([p 'O123.bin'],O123);


