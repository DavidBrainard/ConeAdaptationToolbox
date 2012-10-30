function scene = sceneCheckerboard(varargin)
% function scene = sceneCheckerboard(varargin)
% Arguments must be ('parm', val) pairs.
% Possible parameters:
%   'scene' - scene.type 
%   'checkWavelength' - wavelength of light that will be at different
%                       intensities at each check
%   'checkIntensity' - the matrix of itensity vales of each square in the
%                      checkerboard. Location in matrix corresponds to check's location in
%                      board
%   'imageSizeDeg' - total size of the image in degrees
%   'imageSizePixels' - total size of the image in pixels
%   'boardSizePixels' - length of one side of the checkerboard in pixels
%   'checkXOffsetDeg' - changes horizontal location of checkerboard
%   'checkYOffsetDeg' - changes vertical location of checkerboard
%   'distanceToCheckerboardMeters' - 
% Replaced: function scene = sceneCheckerboard(scene,checkWavelength,checkLengthDeg,checkIntensity,imageSizeDeg,imageSizePixels,checkXOffsetDeg,checkYOffsetDeg,distanceToCheckerboardMeters)
% See sceneMonospot
% Create a scene describing a monochromatic checkerboard.
%
% The scene is assumed to be square.
% This should no longer be used for analyzing data from Radonjic 2011, as
% it does not take into account the device measurements, adjusts luminance,
% and some other issues. However, it is being kept in this folder in case
% someon wants to make a generic checkerboard some day. -ekf
%
% 7/22/11  dhb  Wrote it.
% 7/25/11  bw   A little more ISET background added
% 8/7/12   ekf  Edited for checkerboard

%% Set defaults
% I changed imageSizeDeg to a smaller default (2)
% Otherwise the image sampling is very low.
scene.type = 'scene'; 
checkWavelength = 500; 
checkIntensity = [2,1;3,4]; 
imageSizeDeg = 4.1*5;   % check size degrees * checks per side
imageSizePixels = 400; 
boardSizePixels = 200;
checkXOffsetDeg = 0; 
checkYOffsetDeg = 0; 
distanceToCheckerboardMeters = 1.2; 

%% Change values based on input
if ~isempty(varargin)
    if rem(length(varargin),2)~=0
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        parm=varargin{ii};
        val=varargin{ii+1};
        switch parm
            case 'scene'
                scene.type = val;

            case 'checkWavelength'
                checkWavelength = val;

            case 'checkIntensity'
                checkIntensity = val;

            case 'imageSizeDeg'
                imageSizeDeg = val;

            case 'imageSizePixels'
                imageSizePixels = val;
                
            case 'boardSizePixels'
                boardSizePixels = val;
                
            case 'checkXOffsetDeg'
                checkXOffsetDeg = val;

            case 'checkYOffsetDeg'
                checkYOffsetDeg = val;

            case 'distanceToCheckerboardMeters'
                distanceToCheckerboardMeters = val;

        end
    end
end

if boardSizePixels > imageSizePixels
    error('Error: Board size larger than image size')
end

%% Also, let's use sceneCreate to get the right default structure.
scene = sceneCreate;  % The data are wrong but the structure is right

% Set basic parameters
scene = sceneSet(scene,'name',sprintf('check-%d',checkWavelength));
scene = sceneSet(scene,'distance',distanceToCheckerboardMeters);
scene = sceneSet(scene,'wAngular',imageSizeDeg);

wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');

%% Get check center and size in pixels
pixelsPerDeg = imageSizePixels/imageSizeDeg;
checkXCenterPixels = round(imageSizePixels/2 + checkXOffsetDeg*pixelsPerDeg);
checkYCenterPixels = round(imageSizePixels/2 + checkYOffsetDeg*pixelsPerDeg);
topRightCorner = [checkYCenterPixels-boardSizePixels/2,checkXCenterPixels-boardSizePixels/2];


%% Set check in hyperspectal photon image
%
% Black is more than zero to prevent HDR problem with ieCompressData 
%
% Replaced: wlIndex = find(wave == checkWavelength);
wlIndex  = ieFindWaveIndex(wave,checkWavelength);
% if (isempty(wlIndex))
%     error('Oops: Desired check wavelength not one of the sampled wavelengths');
% end
% if (length(wlIndex) > 1)
%     error('Very peculiar: duplicate wavelengths in sampled wavelengths');
% end

photons = 0*ones(imageSizePixels,imageSizePixels,nWave)*1e-6;
theMonoPlane = photons(:,:,wlIndex);

% Make sure checkerboard is square
nchecks=size(checkIntensity,1)*size(checkIntensity,2);
if rem(sqrt(nchecks),1) ~= 0
    error('Checkerboard must be square')
end

checkSize=floor(boardSizePixels/sqrt(nchecks));

% throw an error if board doesn't fit
if topRightCorner(1)+boardSizePixels > imageSizePixels || topRightCorner(2)+boardSizePixels > imageSizePixels
    error('Error: Board will not fit within image at given location')
end
% insert values into monoplane
for row=1:sqrt(nchecks)
    for col=1:sqrt(nchecks)
        theMonoPlane(topRightCorner(1)+(row-1)*checkSize+1:topRightCorner(1)+row*checkSize,...
            topRightCorner(2)+(col-1)*checkSize+1:topRightCorner(2)+col*checkSize) = checkIntensity(row,col);
    end
end


photons(:,:,wlIndex) = theMonoPlane;

scene = sceneSet(scene,'cphotons',photons);

% Should set the mean luminance here
% This should probably be a parameter that is set at the input
% I have been meaning to write a set 'peak' luminance rather than 'mean'
% for cases just like this one.

meanL = 1; % cd/m2 - most is black so the check is pretty bright.
scene = sceneAdjustLuminance(scene,meanL);
    
return;
