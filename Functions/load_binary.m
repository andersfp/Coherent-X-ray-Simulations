function A = load_binary(filename,size,precision)
% Load the file 'filename' with size 'size', which was save with precision
% 'precision'. If precision is not specified use 'double'. The result is
% returned in 'A'.

% Set default precision
if nargin == 2
    precision = 'double';
end

% Define the data sizes of precisions
switch precision
    case {'double','float64'}
        ds = 8;
    case {'single','float','float32'}
        ds = 4;
    case {'int','uint'}
        ds = 4;
    case {'int64','uint64'}
        ds = 8;
    case {'int32','uint32','long'}
        ds = 4;
    case {'int16','uint16','short'}
        ds = 2;
    case {'int8','uint8'}
        ds = 1;
end

% Get the size of the file
lst = dir(filename);
if isempty(lst)
    error('File not found.');
end
s = lst.bytes;

% Check the size of the file is compatible with the size
if prod(size)*ds == s
    c = 1;
elseif prod(size)*ds*2 == s
    c = 2;
else
    error('Input size and precision incompatible with file size.');
end

% Set the 4th dimension of size to be c
switch numel(size)
    case 1
        size(2) = 1;
        size(3) = 1;
        size(4) = c;
    case 2
        size(3) = 1;
        size(4) = c;
    case 3
        size(4) = c;
end

% Open the file for reading
fid = fopen(filename,'r');

% Read the data
A = fread(fid,Inf,['*' precision]);

% Close the file
fclose(fid);

% Reshape the data
A = reshape(A,size);

% Separate real and imaginary components if complex variable
if c == 2
    A = complex(A(:,:,:,1),A(:,:,:,2));
end


