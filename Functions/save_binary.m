function fc = save_binary(filename,A,precision)
% Save a file named 'filename' with the content of 'A' in binary format
% with the precision 'precision'. If precision is not specified, the
% precision of 'A' will be used. 'A' can have up to 3 dimensions.

% Get default precision if not specified
if nargin == 2
    precision = class(A);
end

% Reformat the array if complex
if ~isreal(A)
    A = cat(4,real(A),imag(A));
end

% Open the file for saving
fid = fopen(filename,'w');

% Save the data in a binary stream
count = fwrite(fid,A,precision);

% Check that the file was written properly
if count ~= numel(A)
    warning('Possible problem in writing the file.');
end

% Close the file
fc = fclose(fid);

