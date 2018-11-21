function [h,data_mean,data_std] = maskNeighbor3(mask,data_input,m,threshold)
% Calculate the mean and standard deviation of the m nearest unmasked
% pixels. Expand the ouput to linear voxel indices.

% Get the size of the data
nx = size(data_input,2);
ny = size(data_input,1);
n_omega = size(data_input,3);

% Set the threshold (in case of non-binary mask)
if nargin == 3
    threshold = 0.5;
end

% Get the indices of the unmasked pixels on the detector
hd = find(mask >= threshold);
[id,jd] = ind2sub(size(mask),hd);
id = repmat(id,n_omega,1);
jd = repmat(jd,n_omega,1);
kd = repmat(1:n_omega,length(hd),1);
kd = kd(:);

% Get the indices of the masked pixels on the detector
hm = find(mask < threshold);
[im,jm] = ind2sub(size(mask),hm);
im = repmat(im,n_omega,1);
jm = repmat(jm,n_omega,1);
km = repmat(1:n_omega,length(hm),1);
km = km(:);

% Get the indices of the m nearest unmasked pixels for each masked pixel
hn = knnsearch([id jd kd],[im jm km],'K',m);

% Convert the indices of indices (hn) to subscripts of pixel numbers
ii = id(hn);
jj = jd(hn);
kk = kd(hn);

% Convert the pixel subscripts to indices
hh = sub2ind([ny nx n_omega],ii,jj,kk);

% Calculate the mean and standard deviation of the group of nearest pixels
dat = data_input(hh);
data_mean = mean(dat,2);
data_std = std(dat,0,2);

% Convert indices from 2D to 3D
h = sub2ind([ny nx n_omega],im,jm,km);


