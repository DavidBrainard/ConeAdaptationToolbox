function scene = sceneMonospot(scene,spotWavelength,spotIntensity,spotDiameterDeg,imageSizeDeg,imageSizePixels,spotXOffsetDeg,spotYOffsetDeg,distanceToSpotMeters)
% scene = sceneMonospot(scene,spotWavelength,spotDiameterDeg,spotIntensity,imageSizeDeg,imageSizePixels,spotXOffsetDeg,spotYOffsetDeg,distanceToSpotMeters)
%
% Create a scene describing a monochromatic spot.
%
% The scene is assumed to be square.
%
% 7/22/11  dhb  Wrote it.
% 7/25/11  bw       A little more ISET background added

%% Set defaults
% I changed imageSizeDeg to a smaller default (2)
% Otherwise the image sampling is very low.
if ieNotDefined('scene'), scene.type = 'scene'; end
if ieNotDefined('spotWavelength'), spotWavelength = 500; end
if ieNotDefined('spotIntensity'), spotIntensity = 1; end
if ieNotDefined('spotDiameterDeg'), spotDiameterDeg = 1; end
if ieNotDefined('imageSizeDeg'), imageSizeDeg = 2; end
if ieNotDefined('imageSizePixels'), imageSizePixels = 64; end
if ieNotDefined('spotXOffsetDeg'), spotXOffsetDeg = 0; end
if ieNotDefined('spotYOffsetDeg'), spotYOffsetDeg = 0; end
if ieNotDefined('distanceToSpotMeters'), distanceToSpotMeters = 1.2; end

%% You might consider setting the spotIntensity in units
%  I will convert them to units below.
%  Same for other parameters.

%% Also, let's use sceneCreate to get the right default structure.
scene = sceneCreate;  % The data are wrong but the structure is right

% Set basic parameters
scene = sceneSet(scene,'name',sprintf('spot-%d',spotWavelength));
scene = sceneSet(scene,'distance',distanceToSpotMeters);
scene = sceneSet(scene,'wAngular',imageSizeDeg);

wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');

%% Get spot center and size in pixels
pixelsPerDeg = imageSizePixels/imageSizeDeg;
spotXCenterPixels = round(imageSizePixels/2 + spotXOffsetDeg*pixelsPerDeg);
spotYCenterPixels = round(imageSizePixels/2 + spotYOffsetDeg*pixelsPerDeg);
spotRadiusPixels  = round(pixelsPerDeg*spotDiameterDeg/2);
if (spotRadiusPixels == 0)
    spotRadiusPixels = 1;
end

%% Set spot in hyperspectal photon image
%
% Black is more than zero to prevent HDR problem with ieCompressData 
%
% Replaced: wlIndex = find(wave == spotWavelength);
wlIndex  = ieFindWaveIndex(wave,spotWavelength);
% if (isempty(wlIndex))
%     error('Oops: Desired spot wavelength not one of the sampled wavelengths');
% end
% if (length(wlIndex) > 1)
%     error('Very peculiar: duplicate wavelengths in sampled wavelengths');
% end

[xCoords,yCoords] = meshgrid(1:imageSizePixels);
index =  ((xCoords-spotXCenterPixels).^2 + (yCoords-spotYCenterPixels).^2) <= spotRadiusPixels.^2;

photons = 0*ones(imageSizePixels,imageSizePixels,nWave)*1e-6;
theMonoPlane = photons(:,:,wlIndex);
theMonoPlane(index) = spotIntensity;
photons(:,:,wlIndex) = theMonoPlane;

scene = sceneSet(scene,'cphotons',photons);

% Should set the mean luminance here
% This should probably be a parameter that is set at the input
% I have been meaning to write a set 'peak' luminance rather than 'mean'
% for cases just like this one.

meanL = 1; % cd/m2 - most is black so the spot is pretty bright.
scene = sceneAdjustLuminance(scene,meanL);
    
return;
