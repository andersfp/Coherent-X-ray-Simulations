function fileList = dxm_makeFileList(filePath,filePrefix)
% Hugh Simons, Technical Univeristy of Denmark
% 3rd March 2016

dirList = dir(strcat(filePath,filePrefix,'.edf')); if isempty(dirList); fprintf('Cannot locate files. Exiting...\n'); return; end

nFiles = length(dirList);

for iFile = 1:nFiles
    fileList{iFile} = strcat(filePath,dirList(iFile).name);
end

end

