function[out]=tri(x)
% 
% Triangle function
% 
% Evaluates tri(x)
% 

% Create lines
t = 1 - abs(x);

% keep lines for |x|<=1, out=0 otherwise
mask = abs(x) <= 1;
out = t.*mask;

end
