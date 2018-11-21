function [Data] = dxm2d_load(fileList,Geometry,Settings,varargin)
% Hugh Simons, Technical Univeristy of Denmark
% 3rd March 2016

%Preliminary stuff
Data.Geometry = Geometry;
MotorList = varargin;
nMotors = length(MotorList);
if isfield(Settings,'roi')
    ny = Settings.roi(3) - Settings.roi(1) + 1; %old: ny = Settings.roi(1) - Settings.roi(3);
    nx = Settings.roi(4) - Settings.roi(2) + 1; %old: nx = Settings.roi(2) - Settings.roi(4);
else
    [ny, nx] = size(Settings.background);
end

%Load images to stack
nFiles = length(fileList);
Data.intensity = uint16(zeros(ny,nx,nFiles));
for iFile = 1:nFiles
    fprintf('Loading file #%d...',iFile);
    %Get the image and apply filters
    image = edf_read(fileList{iFile}) - Settings.background;
    Data.intensity(:,:,iFile) = filterImage(image,Settings);
    %Get the relevant motor positions for it
    info = edf_info(fileList{iFile});
    for iMotor = 1:nMotors
    	MotorList{iMotor}.positionList(iFile) = eval(sprintf('info.motor.%s',MotorList{iMotor}.name));
    end
    fprintf('Done!\n');
end

%Process scaling
if strcmp(Geometry.orientation,'vert')
    yfov = ny*Geometry.yPixelSize/Geometry.xMag/tand(Geometry.twotheta);
    xfov = nx*Geometry.xPixelSize/Geometry.xMag;
elseif strcmp(Geometry.orientation,'horz')
    yfov = ny*Geometry.yPixelSize/Geometry.xMag;
    xfov = nx*Geometry.xPixelSize/Geometry.xMag/tand(Geometry.twotheta);
else
    fprintf('Geometry orientation not specified. Exiting...\n');
    return
end
Data.threshold = Settings.threshold;
Data.twotheta = Geometry.twotheta;
Data.y = linspace(-yfov/2,yfov/2,ny);
Data.x = linspace(-xfov/2,xfov/2,nx);

%Process dimensionality
if MotorList{1}.isZap
    interval = (MotorList{1}.end - MotorList{1}.start) / MotorList{1}.nIntervals;
    Data.dim1.angles = linspace(MotorList{1}.start + interval/2, MotorList{1}.end - interval/2, MotorList{1}.nIntervals);
else
    Data.dim1.angles = unique(MotorList{1}.positionList);
end
Data.dim1.name = MotorList{1}.name;

Data.shape = [ny,nx,length(Data.dim1.angles)];
for iMotor = 2:nMotors
    iDim = sprintf('dim%d',iMotor);
    Data.(iDim).angles = unique(MotorList{iMotor}.positionList);
    Data.shape(end+1) = length(Data.(iDim).angles);
    Data.(iDim).name = MotorList{iMotor}.name;
end
Data.intensity = reshape(Data.intensity,Data.shape);

end

function fimg = filterImage(image,Settings)
if isfield(Settings,'roi')
	roi = Settings.roi;
	image = image(roi(1):roi(3),roi(2):roi(4));
end
if isfield(Settings,'median')
	image = medfilt2(image,Settings.median);
end    
if isfield(Settings,'murofi')
	image = murofi(image,Settings.murofi(1),Settings.murofi(2));
end
fimg = image;
end

function [fimg] = murofi(img,nstp,slen)
msk = ones(size(img));
stp = 180/nstp;
stack = zeros(size(img,1),size(img,2),nstp);
for n = 1:nstp
    rot = (n-1)*stp;
    imgr = imrotate(img,rot);
    mskr = imrotate(msk,rot);
    for j = 1:size(mskr,1)
        ids = find(mskr(j,:)>0);
        if ~isempty(ids)
            imgr(j,ids) = smooth(imgr(j,ids),slen);
        end
    end
    imgb = imrotate(imgr,-rot);
    idx0 = floor((size(imgb)-size(img))/2)+1;
    idx1 = floor((size(imgb)+size(img))/2);
    imgb = imgb(idx0(1):idx1(1),idx0(2):idx1(2));
    stack(:,:,n) = imgb;
end
fimg = min(stack,[],3);
end