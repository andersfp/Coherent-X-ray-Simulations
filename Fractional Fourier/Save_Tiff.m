function Save_Tiff(f,I,norm)
% Save the intenisty as an 8 bit TIFF file with absolute path name f. The
% intensity is normalized to norm.

% Get the size of the array
m = size(I,1);
n = size(I,2);

% Normalize the intensity
I = I*(255/norm);

% Convert the intensity to uint16
I = uint8(I);

% Open the tiff object
obj = Tiff(f,'w');

% Set the saving options
obj.setTag('Photometric',Tiff.Photometric.MinIsBlack);
obj.setTag('Compression',Tiff.Compression.LZW);
obj.setTag('BitsPerSample',8);
obj.setTag('SamplesPerPixel',1);
obj.setTag('SampleFormat',Tiff.SampleFormat.UInt);
obj.setTag('ImageLength',m);
obj.setTag('ImageWidth',n);
obj.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);

% Write the image file
write(obj,I);

% Close the object
obj.close();

