% Initialization
clear;
close all;
clc;


%% Load the data
% Load the experimental parameters
load('Exp_Param.mat');

% Load the actual object
Oxyz = load_binary([p 'Oxyz.bin'],[ny nx n_omega]);
O123 = load_binary([p 'O123.bin'],[ny nx n_omega]);

% Load the 3D electric fields
Ef = load_binary([p 'Ef.bin'],[ny nx n_omega]);
El = load_binary([p 'El.bin'],[ny nx n_omega]);

% Load the free space reconstruction
Efr = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_field.bin'],[ny nx n_omega]);
Of = load_binary([p 'Reconstruction_Free_Space_Shrinkwrap_1_object_xyz.bin'],[ny nx n_omega]);

% Load the lens reconstruction
Elr = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_field.bin'],[ny nx n_omega]);
Ol = load_binary([p 'Reconstruction_Lens_Shrinkwrap_1_object_xyz.bin'],[ny nx n_omega]);


%% Save the data
% Save the true object
A = Oxyz;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Original_Object.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the true object in conjugated measurement space
A = O123;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Original_Object_123.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the reconstructed amplitude
A = abs(Of);
A = A/max(A(:));
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Reconstruction_Amplitude_Free_Space_Shrinkwrap_1.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the reconstructed phase
A = angle(Of);
A(A == 0) = NaN;
A = A + pi;
A = A/(2*pi);
A(isnan(A)) = 0;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Reconstruction_Phase_Free_Space_Shrinkwrap_1.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the reconstructed amplitude
A = abs(Ol);
A = A/max(A(:));
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Reconstruction_Amplitude_Lens_Shrinkwrap_1.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the reconstructed phase
A = angle(Ol);
A(A == 0) = NaN;
A = A + pi;
A = A/(2*pi);
A(isnan(A)) = 0;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Reconstruction_Phase_Lens_Shrinkwrap_1.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the reconstructed intensity
A = abs(Efr).^2;
A = A/max(A(:));
A = log10(A);
A = A + 12;
A(A < 0) = 0;
A = A/12;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Reconstruction_Intensity_Free_Space_Shrinkwrap_1.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the reconstructed intensity
A = abs(Elr).^2;
A = A/max(A(:));
A = log10(A);
A = A + 12;
A(A < 0) = 0;
A = A/12;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Reconstruction_Intensity_Lens_Shrinkwrap_1.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the simulated intensity
A = abs(Ef).^2;
A = A/max(A(:));
A = log10(A);
A = A + 12;
A(A < 0) = 0;
A = A/12;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Simulation_Intensity_Free_Space_Shrinkwrap.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

% Save the simulated intensity
A = abs(El).^2;
A = A/max(A(:));
A = log10(A);
A = A + 12;
A(A < 0) = 0;
A = A/12;
A = A*65535;
A = round(A);
A = uint16(A);
fid = fopen([p '3D_Data\Simulation_Intensity_Lens_Shrinkwrap.raw'],'w');
fwrite(fid,A,'uint16');
fclose(fid);

disp([delta_ry delta_rx delta_rz]*1e9);


