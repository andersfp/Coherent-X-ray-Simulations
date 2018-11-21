function [h,data_mean,data_std] = maskNeighbor2(mask,data_input,m,threshold)
% Calculate the mean and standard deviation of the m nearest unmasked
% pixels. Expand the ouput to linear voxel indices.

% Get the size of the data
nx = size(data_input,2);
ny = size(data_input,1);
n_omega = size(data_input,3);

% Set the threshold (in case of non-binary mask) if not specified
if nargin == 3
    threshold = 0.5;
end

% Get the indices of the unmasked pixels on the detector
hd = find(mask >= threshold);
[id,jd] = ind2sub(size(mask),hd);

% Get the indices of the masked pixels on the detector
hm = find(mask < threshold);
[im,jm] = ind2sub(size(mask),hm);

% Get the indices of the m nearest unmasked pixels for each masked pixel
hn = knnsearch([id jd],[im jm],'K',m);

% Convert the indices of indices (hn) to subscripts of pixel numbers
ii = id(hn);
jj = jd(hn);

% Convert the pixel subscripts to indices
hh = sub2ind([ny nx],ii,jj);

% Calculate the mean and standard deviation of the group of nearest pixels
data_mean = zeros(length(hm),n_omega,'like',data_input);
data_std = data_mean;
if isempty(gcp('nocreate'))
    for i = 1:n_omega
        temp = data_input(:,:,i);
        data_mean(:,i) = mean(temp(hh),2);
        data_std(:,i) = std(temp(hh),0,2);
    end
else
    parfor i = 1:n_omega
        temp = data_input(:,:,i);
        data_mean(:,i) = mean(temp(hh),2);
        data_std(:,i) = std(temp(hh),0,2);
    end
end

% Convert indices from 2D to 3D
Im = repmat(im,1,n_omega);
Jm = repmat(jm,1,n_omega);
km = 1:n_omega;
Km = repmat(km,length(hm),1);
h = sub2ind([ny nx n_omega],Im,Jm,Km);

% Flatten output arrays
h = h(:);
data_mean = data_mean(:);
data_std = data_std(:);

