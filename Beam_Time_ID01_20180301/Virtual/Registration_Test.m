% Initialization
clear;
close all;
clc;


%% Load data
% Load the shifted data
load('Processed_Virtual_Shifted_2.mat');

% Extract the central images
i0 = size(I,3)/2 + 1;
im = log10(squeeze(I(:,:,i0,:)));
im(isinf(im)) = 0;


%% Panorama stitching
% First image
points = detectSURFFeatures(im(:,:,1));
[features, points] = extractFeatures(im(:,:,1), points);

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = size(im,3);
tforms(numImages) = affine2d(eye(3));

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);

% Iterate over remaining image pairs
for i = 2:numImages

    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;

    % Save image size.
    imageSize(i,:) = size(im(:,:,i));

    % Detect and extract SURF features for I(n).
    points = detectSURFFeatures(im(:,:,i));
    [features, points] = extractFeatures(im(:,:,i), points);

    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);

    % Estimate the transformation between I(n) and I(n-1).
    tforms(i) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'affine', 'Confidence', 99.9, 'MaxNumTrials', 2000);

    % Compute T(n) * T(n-1) * ... * T(1)
    tforms(i).T = tforms(i).T * tforms(i-1).T;
end


xlim = zeros(numImages,2);
ylim = xlim;
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
end

maxImageSize = max(imageSize);

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width], 'like', im);


blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for i = 1:numImages

    % Transform I into the panorama.
    warpedImage = imwarp(im(:,:,i), tforms(i), 'OutputView', panoramaView);

    % Generate a binary mask.
    mask = imwarp(true(size(im(:,:,i),1),size(im(:,:,i),2)), tforms(i), 'OutputView', panoramaView);

    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, mask);
end

figure;
imshow(panorama);


