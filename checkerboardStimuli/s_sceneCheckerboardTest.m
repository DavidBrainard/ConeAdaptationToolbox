function [oiPhotons, intensity]=s_sceneCheckerboardTest(varargin)
% function [oiPhotons, intensity]=s_sceneCheckerboardTest(varargin)
%
% Excercise ISET, as well as test local routines.
% experiment - 'one' or 'two'; Chooses data and arrangement from specified experiment
% ExpCol - Column number of experiment/condition data to be used 
% midSquareNum - matrix row number of the square that the test square should match 
% midSquareCol - Column of the experiment that the middle square's value
%       will be pulled from in getCheckerboardRadonjic2011Data. Usually the same
%       as the ExpCol
%
% 7/14/11  dhb, gt  Started down this road
% 7/22/11  dhb      Now working for initial monospot.
% 7/25/11  bw       A little more ISET background added
% 8/7/12   ekf      Edited for use with checkerboard


%% Initialize.  Close ISET windows effectively and open ISET again.  
% The command ieMainW will bring up the main window if you want it (probably you don't).
s_initISET;  

%% Set defaults
experiment='two';
ExpCol=2;
midSquareNum=5;
midSquareCol=ExpCol;

%% Change values based on input
if ~isempty(varargin)
    if rem(length(varargin),2)~=0
        error('Arguments must be (pair, val) pairs');
    end
    MSCDefined=0;
    for ii=1:2:(length(varargin)-1)
        parm=varargin{ii};
        val=varargin{ii+1};
        switch parm
            case 'experiment'
                experiment = val;

            case 'ExpCol'
                ExpCol = val;

            case 'midSquareNum'
                midSquareNum = val;
                
            case 'midSquareCol'
                midSquareCol = val;
                MSCDefined = 1;
        end
    end
    if ~MSCDefined
        midSquareCol=ExpCol;
    end
end
%% Bring in data from getCheckerboardRadonjic2011Data
Checkerboard=getCheckerboardRadonjic2011Data('experiment', experiment, 'ExpCol', ExpCol,'midSquareNum', midSquareNum,'midSquareCol',midSquareCol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Make stimulus of radiance photons %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load calibration data
% These live in PsychCalLocalData on our systems.  We
% can approximate the stimuli as a linear combination
% of the front screen HDR primaries.  This isn't exactly
% right, becasue there are various second order effects.
% But this is plenty close enough for analysis purposes.
cal = LoadCalFile('HDRFrontMondrianfull');
S = [400 10 31];
B = SplineSpd(cal.S_device,cal.P_device(:,1:3),S);

%% Load XYZ color matching functions
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
M_WgtsToXYZ = T_xyz*B;
M_XYZToWgts = inv(M_WgtsToXYZ);

%% Specify luminance and chromaticity of stimuli
% We'll just do two spectra here.  For the others
% just vary luminance while holding chromaticity constant.
%
% We construct a spectrum that has specified luminance and
% chromaticity.  Because of the way T_xyz was scaled, this
% spectrum is radiance in units of Watts/[sr-m2-wlband].
stimulus_xy = [0.309 0.338]';
% Luminance takes input as a row vector, so all the luminances are split
% up from the matrix into one long row vector
stimulus_Y = [Checkerboard(1,:) Checkerboard(2,:) Checkerboard(3,:) Checkerboard(4,:) Checkerboard(5,:)];
nStimuli = length(stimulus_Y);

for s = 1:nStimuli
    stimulus_xyY(:,s) = [stimulus_xy ; stimulus_Y(s)];
    stimulus_XYZ(:,s) = xyYToXYZ(stimulus_xyY(:,s));
    stimulus_Wgts(:,s) = M_XYZToWgts*stimulus_XYZ(:,s);
    stimulus_spd(:,s) = B*stimulus_Wgts(:,s);   %Watts/[sr-m2-wlband]
    stimulus_XYZCheck(:,s) = T_xyz*stimulus_spd(:,s);
    stimulus_xyYCheck(:,s) = XYZToxyY(stimulus_XYZCheck(:,s)); 
end


% Converts so that input units are what VSET wants (photons/s/m2/sr/nm) from
% Watts/[sr-m2-wlband]
% EnergyToQuanta takes input energy in nm
photons = EnergyToQuanta(S,stimulus_spd/S(2));  % returns photons/s/m2/sr/nm 

stimulus_photons(:,:)=[photons(1,1:5); photons(1,6:10); photons(1,11:15);...
    photons(1, 16:20); photons(1, 21:25)];    %initialization
% Rearranges matrix to make its geometry match VSET's desire geometry
for wavelength=2:S(3)
    temp1(1,:)=photons(wavelength,1:5);
    temp1(2,:)=photons(wavelength,6:10);
    temp1(3,:)=photons(wavelength,11:15);
    temp1(4,:)=photons(wavelength,16:20);
    temp1(5,:)=photons(wavelength,21:25);
    stimulus_photons=cat(3,stimulus_photons,temp1);
end
% Expand matrix according to how many cones recieve input from one check.
% This must be a perfect square. VSET takes each entry to be one pixel, and
% the valeton calculation assumes each pixel is one photoreceptor.
% Therefore, you must create a new entry for each photoreceptor you want,
% while maintaining the luminances of the checkerboard. This is done by
% repeating values inside the original matrix.
conesPerCheck=1;   % Value of 36 is relatively accurate for this stimulus size
sqrtCones=sqrt(conesPerCheck);

if rem(sqrtCones,1)~=0
    error('Error: number of cones per check must be an even square root')
end
for g=1:S(3)
    for k=0:sqrtCones:5*sqrtCones-1
        for j=0:sqrtCones:5*sqrtCones-1
            stimulus_photons_big(k+1:k+sqrtCones,j+1:j+sqrtCones,g)=...
                stimulus_photons(floor(k/sqrtCones)+1,floor(j/sqrtCones)+1,g);
        end
    end
end

stimulus_photons=stimulus_photons_big;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create an example scene, using an ISET preset
% We can add this to sceneCreate when you are ready.  Comments within
% sceneCheckerboard.
% Also, let's use sceneCreate to get the right default structure.

scene = sceneCreate;  % The data are wrong but the structure is right
scene.type = 'scene'; 
imageSizeDeg = 4.1*5;   % each check is 4.1 degrees * 5 checks per side 
distanceToCheckerboardMeters = 1.2; 

% Set basic parameters
scene = sceneSet(scene,'distance',distanceToCheckerboardMeters);
scene = sceneSet(scene,'wAngular',imageSizeDeg);
scene = sceneSet(scene,'cphotons',stimulus_photons);

vcAddAndSelectObject(scene); sceneWindow;


%% Define an optical system that models the human eye
oi = oiCreate('human');
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi);
oiWindow;

% Return the photon distribution matrix of the optical image
oiPhotons = oiGet(oi,'cphotons');   % Now in phot/s/m2/nm
oiPhotons11 = oiPhotons(:,:,11);    % Temporarily remove other wavelengths for plotting purposes

% plot
j=bar3(oiPhotons11);  % Input in radiance photons
% Tell handle graphics to use interpolated rather than flat shading
shading interp
colormap cool
% For each barseries, map its CData to its ZData
for i = 1:length(j)
    zdata = get(j(i),'ZData');
    set(j(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(j,'EdgeColor','k') 
end
title('Input in radiance photons')
%% Convert quantized output into trolands
S = [400 10 31]; 

intensity=zeros(size(oiPhotons11,1),size(oiPhotons11,1));   %Preallocation
quanta=zeros(1,S(3)); %Preallocation
[nRow,nCol]=size(intensity);

% Calculate intensity in trolands one pixel at a time
for k=1:nRow
    for j=1:nCol
        for l=1:S(3)
            quanta(l)= double(oiPhotons(k,j,l));
        end
        energy=QuantaToEnergy(S,quanta'); % now in watts/m2/nm
        trolands=RetIrradianceToTrolands(S(2)*energy/1e12,S,'Photopic','Human'); % converts input to watts/um2/wli
        intensity(k,j) = trolands;
    end
end

