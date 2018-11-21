function edf_write(A,filename)

if nargin<2
    error('Too few inputs!')
end

dim = size(A);

%% more or less random header...
txt = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n','{',...
                'HeaderID       = EH:000001:000000:000000 ;', ...
                'Image          = 1 ;', ...
                'ByteOrder      = LowByteFirst ;', ...
                'DataType       = UnsignedShort ;', ...
                'Dim_1          = 2048;', ...
                'Dim_2          = 2048;');
header = uint8(32*ones(1,4096));
header(1:length(txt)) = txt;
header(end-1:end) = [125 10]; % '}' and '\n'

%% handle data
if ~isa(A,'uint16')
    warning('A is not a uint16 so conversion errors may happen!');
end
A = fliplr(flipud(A)');
data = uint16(A(:));

%% check file
dowrite = 1;
if exist(filename,'file');
    overwrt = questdlg('File exists - overwrite?', 'File exist', ...
        'Yes', 'No', 'No');
    if strcmp(overwrt,'No')
        dowrite = 0;
    end
end

%% write file or not
if ~dowrite
    fprintf('%s\n','No file written.')
else
    fp = fopen(filename,'w');
    fwrite(fp,header,'uint8');
    fwrite(fp,data,'uint16');
end
    