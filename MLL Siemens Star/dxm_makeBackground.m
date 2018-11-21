function background = dxm_makeBackground(fileList,type)
% Hugh Simons, Technical Univeristy of Denmark
% 3rd March 2016

nFiles = length(fileList);
for iFile = 1:nFiles
    imageStack(:,:,iFile) = edf_read(fileList{iFile});
end

switch type
    case 'mean'
        background = mean(imageStack,3);
    case 'median'
        background = median(imageStack,3);
    case 'min'
        background = min(imageStack,[],3);
    case 'max'
        background = max(imageStack,[],3);
    otherwise
        fprintf('Type not recognized. Exiting... \n');
        background = nan;
        return
end

end

