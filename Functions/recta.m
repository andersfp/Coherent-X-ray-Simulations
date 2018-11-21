function out = recta(x)
% Rectangle function. Returns 1 for x >= -1/2 & x < 1/2, 0 otherwise.
% 
% Example of usage:
% out = recta(x)
%

out = double(x >= -1/2-eps(0.5) & x < 1/2-eps(0.5));

