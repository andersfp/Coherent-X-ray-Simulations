function [f,F] = fileNames(p,fp)
% Find all file names matching the pattern fp in the folder p. Returns the
% file names alone in f and full paths in F.

% Search the directory to find the files
lst = dir([p fp]);

% Convert the struct to a cell
lst = struct2cell(lst);

% Extract the file names only
f = lst(1,:).';

% Combine the path and file names
F = cellfun(@(x) [p x],f,'UniformOutput',0);

