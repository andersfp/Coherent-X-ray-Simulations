function [edfheader,msg]=edf_info(fname)
% EDF_INFO.M
% provides a structure containing all elements found in the header of the
% specified .edf file.  
% Usage:
% edfheader=edf_info(filename)
% or, more commonly imfinfo(filename) after edf_setup_imformats has been run

%
% GJ May 2006
% modified July 2006 to force all parameters to lower case (PyHST generates
% a number of different variations!

msg=[];
edfheader=[];
pattern='(?<parameter>.*)=(?<value>.*);';
if ~exist(fname,'file')
  error(sprintf('Could not find %s',fname))
end
fid=fopen(fname,'rt');
if fid==-1
  error(sprint('Could not open %s',fname))
end
lengthofheader=512;
endofheader=[];

while isempty(endofheader)
  htxt=fgetl(fid);
  % DIRTY FIX TO SPECIFIC PROBLEM HERE
  if strcmp(htxt,'Detector  = frelon4m (sn=29) ;'),
      htxt = 'Detector  = frelon4m (sn 29) ;';
  end
  
  
  endofheader=find(htxt=='}');
  tmp=regexp(htxt,pattern,'names');
  if ~isempty(tmp)
    field=sfCleanField(lower(strtrim(tmp.parameter)));
    edfheader=setfield(edfheader,field,strtrim(tmp.value));
  end
end
edfheader=setfield(edfheader,'headerlength',ftell(fid));
if mod(ftell(fid),512)~=0
    fprintf('Header is not multiple of 512 bytes in %s!\n',fname);
    disp('Setting correctly!')
    edfheader=setfield(edfheader,'headerlength',round(edfheader.headerlength/512)*512);
end
fclose(fid);

edfheader=sfCleanupHeaderInfo(edfheader);
end

function field=sfCleanField(str)
% if header field has weird characters (spaces, parentheses) - this removes
% them
str(find(str==' '))='_';
str(find(str=='('))=[];
str(find(str==')'))=[];
field=str;
end

function edfheader=sfCleanupHeaderInfo(edfheader)

% convert to numbers any fields that I know should not be strings
if isfield(edfheader,'dim_1')
  edfheader.dim_1=str2num(edfheader.dim_1);
end
if isfield(edfheader,'dim_2')
  edfheader.dim_2=str2num(edfheader.dim_2);
end
if isfield(edfheader,'dim_3')
  edfheader.dim_3=str2num(edfheader.dim_3);
end
if isfield(edfheader,'size')
  edfheader.size=str2num(edfheader.size);
end

if isfield(edfheader,'count_time')
  edfheader.count_time=str2num(edfheader.count_time);
end

if isfield(edfheader','timestamp_ms')
  edfheader.timestamp_ms=str2num(edfheader.timestamp_ms);
end

if isfield(edfheader,'col_beg')
  edfheader.col_beg=str2num(edfheader.col_beg);
end

if isfield(edfheader,'col_end')
  edfheader.col_end=str2num(edfheader.col_end);
end
if isfield(edfheader,'row_beg')
  edfheader.row_beg=str2num(edfheader.row_beg);
end

if isfield(edfheader,'row_end')
  edfheader.row_end=str2num(edfheader.row_end);
end

if isfield(edfheader,'col_bin')
  edfheader.col_bin=str2num(edfheader.col_bin);
end

if isfield(edfheader,'row_bin')
  edfheader.row_bin=str2num(edfheader.row_bin);
end
if isfield(edfheader,'energy')
  edfheader.energy=str2num(edfheader.energy);
end
if isfield(edfheader,'optic_used')
  edfheader.optic_used=str2num(edfheader.optic_used);
end

% convert all the motors to individual fields
if isfield(edfheader,'motor_pos')
  remainder_pos=edfheader.motor_pos;
  remainder_name=edfheader.motor_mne;
  motortmp=[];
  while ~isempty(remainder_name)
    [token_name,remainder_name]=strtok(remainder_name);
    [token_pos,remainder_pos]=strtok(remainder_pos);
    motortmp=setfield(motortmp,token_name,str2num(token_pos));
  end
  edfheader=setfield(edfheader,'motor',motortmp);
  edfheader=rmfield(edfheader,{'motor_pos','motor_mne'});
end


% convert all the counters to individual fields
if isfield(edfheader,'counter_pos')
  remainder_pos=edfheader.counter_pos;
  remainder_name=edfheader.counter_mne;
  countertmp=[];
  while ~isempty(remainder_name)
    [token_name,remainder_name]=strtok(remainder_name);
    [token_pos,remainder_pos]=strtok(remainder_pos);
    countertmp=setfield(countertmp,token_name,str2num(token_pos));
  end
  edfheader=setfield(edfheader,'counter',countertmp);
  edfheader=rmfield(edfheader,{'counter_pos','counter_mne'});
end

% convert endianness to standard representation
if strcmpi(edfheader.byteorder,'HighByteFirst')
  edfheader.byteorder='b';
elseif strcmpi(edfheader.byteorder,'LowByteFirst')
  edfheader.byteorder='l';
else
  disp('No endian-ness specified');
end
  
% convert datatype to matlab standard representation
switch lower(edfheader.datatype)
  case {'doublevalue'}
    edfheader.datatype='float64';
  case {'float','floatvalue','real'}
    edfheader.datatype='float32';
  case {'unsignedlong','unsignedinteger'}
    edfheader.datatype='uint32';
  case {'unsignedshort'}
    edfheader.datatype='uint16';
  case {'unsignedbyte'}
    edfheader.datatype='uint8';
  case {'signedbyte'}
    edfheader.datatype='int8';
  case {'signedshort'}
    edfheader.datatype='int16';
  case {'signedinteger','signedlong'}
    edfheader.datatype='int32';

  otherwise
    edfheader.datatype=[];
    disp('Type not known')
end
end
