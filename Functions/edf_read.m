function [im,varargout]=edf_read(filename,GT_bb)
% EDF_READ.M Reads images and 3D volumes from .edf files
% Usage:
% data=edf_read(filename[,GT_bb]);
% [data,datatype]=edf_read(filename[,GT_bb]);
% 
%  Designed to be used from imread (see edf_setup_imread.m for how this is
%  done automatically).
%
% if GT_bb is supplied (optional) - only reads the region from [originx
% originy width height] inclusive (only works for 2D EDF files
%
% if a second output variable is given, returns the datatype of the data
%
% At present, returns all data as double - used to be better, but I broke
% it :( 

% Greg Johnson
% June 2006
if nargin==0
  [fname,pname]=uigetfile('*.edf','Select an EDF file');
  filename=fullfile(pname,fname);
end
if ~exist(filename,'file')
  error(sprintf('Could not find %s',filename))
end
info=edf_info(filename);
fid=fopen(filename,'r',info.byteorder);
if fid==-1
  error(sprintf('Could not open %s\n',filename))
end
fseek(fid,info.headerlength,'bof');
if ~isfield(info,'dim_3') % this should be simplified...
  info.dim_3=1;
end


%%%%%%%%%%%%
%ROI

if exist('GT_bb','var')
  if info.dim_3>1
    error('Cannot use an ROI on 3D EDF file')
  end
  
  switch info.datatype
    case {'uint8','int8'}
      nbytes=1;
    case {'uint16','int16'}
      nbytes=2;
    case {'uint32','int32','float32'}
      nbytes=4;
    case 'float64'
      nbytes=8;
  end
	
  % ROI was specified
  pr_str=sprintf('%d*%s',GT_bb(3),info.datatype);
  skip=nbytes*(info.dim_1-GT_bb(3));
  fseek(fid,nbytes*(info.dim_1*(GT_bb(2)-1)+GT_bb(1)-1),0);
  im=fread(fid,[GT_bb(3),GT_bb(4)],pr_str,skip);
  info.dim_1=GT_bb(3);
  info.dim_2=GT_bb(4);
else
%%%%%%%%%%%%%
im=fread(fid,info.dim_1*info.dim_2*info.dim_3,info.datatype);

end
fclose(fid);
im=reshape(im,info.dim_1,info.dim_2,info.dim_3);
if info.dim_3>1 % 3D volume
  im=permute(im,[2 1 3]); % reorganise dimensions
else % 2D image
  im=transpose(im);
end

if nargout==2
  varargout{1}=info.datatype;
end
end
