function Data = dxm2d_map(Data)
% Hugh Simons, Technical Univeristy of Denmark
% 3rd March 2016

ny = Data.shape(1);
nx = Data.shape(2);
nDims = length(Data.shape)-2;

%Preallocation
Data.intensityMap = nan([ny,nx]);
for iDim = 1:nDims
	Data.(sprintf('dim%d',iDim)).map = nan([ny,nx,3]);
end

%for each detector pixel (x,y), each data dimension (motors scanned/rocked)
%are evaluated: 1) intensity is summed along all other dimensions than the
%active one 2) maximum momentum method used to find peak position,
%varriance and skewness 3) these are saved in Data.dimX.map

%Looping
for ix = 1:nx %1019
        disp(['analyzing x-value ' num2str(ix) ' of ' num2str(nx)]);
    for iy = 1:ny %1264
        %reciprocalVol = squeeze(Data.intensity(iy,ix,:,:))-Data.threshold;
        reciprocalVol = squeeze(Data.intensity(iy,ix,:,:)); reciprocalVol = reciprocalVol-min(reciprocalVol(:));
        maxIntensity = max(reciprocalVol(:));      
        if maxIntensity > 0
            Data.intensityMap(iy,ix) = maxIntensity;
            for iDim = 1:nDims
                summedDims = 1:nDims; summedDims(summedDims==iDim) = [];
                intensity = double(squeeze(sumDims(reciprocalVol,summedDims))); 
                intensity = intensity(:)'/sum(intensity);
                angle = Data.(sprintf('dim%d',iDim)).angles;
                angleMean = sum(intensity.*angle);
                angleVariance = sum(intensity.*(angle-angleMean).^2);
                angleSkewness = sum(intensity.*(angle-angleMean).^3);           
                eval(sprintf('Data.dim%d.map(iy,ix,:) = [%g,%g,%g];',iDim,angleMean,angleVariance,angleSkewness));
            end %loop dimensions
        end %intensity positive
    end %loop y
end %loop x
end %function

function out = sumDims(r,dims)
out = r;
    for i = 1:length(dims)
        out = sum(out,dims(i));
    end
end