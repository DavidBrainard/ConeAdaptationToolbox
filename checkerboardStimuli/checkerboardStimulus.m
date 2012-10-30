function [stimulus_photons,intensity]=checkerboardStimulus(varargin)
% function [stimulus_photons,intensity]=checkerboardStimulus(varargin)
%
% experiment - 'one' or 'two'; Chooses data and arrangement from specified experiment
% ExpCol - Column number of experiment/condition data to be used 
% midSquareNum - matrix row number of the square that the test square should match 
% midSquareCol - Column of the experiment that the middle square's value
%       will be pulled from in getCheckerboardRadonjic2011Data. Usually the same
%       as the ExpCol.
%
% Some computations that describe the stimuli
% used by Radonjic et al. (2011).
% Computes irradiance directly without using VSET; hence there is no
% optical blurring.
%
% 8/12/12  dhb  Wrote it.
% 8/14/12  ekf  Edited it for use with radonjic data matrix

%% Set defaults
experiment='two';
ExpCol=4;
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
%% Load matrix data 
Checkerboard=getCheckerboardRadonjic2011Data('experiment',experiment,'ExpCol',ExpCol,'midSquareNum',midSquareNum,'midSquareCol',midSquareCol);

%% Load calibration data
% These live in PsychCalLocalData on our systems.  We
% can approximate the stimuli as a linear combination
% of the front screen HDR primaries.  This isn't exactly
% right, becasue there are various second order effects.
% But this is plenty close enough for analysis purposes.
cal = LoadCalFile('HDRFrontMondrianfull');
S = [400 10 31];
B = SplineSpd(cal.S_device,cal.P_device(:,1:3),S);
figure; clf; plot(SToWls(S),B,'k');
xlabel('Wavelength (nm)');
ylabel('Power (arbitrary units)');

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
stimulus_Y = [Checkerboard(1,:) Checkerboard(2,:) Checkerboard(3,:) Checkerboard(4,:) Checkerboard(5,:)];
nStimuli = length(stimulus_Y);
spectralFig = figure; clf; hold on
for s = 1:nStimuli
    stimulus_xyY(:,s) = [stimulus_xy ; stimulus_Y(s)];
    stimulus_XYZ(:,s) = xyYToXYZ(stimulus_xyY(:,s));
    stimulus_Wgts(:,s) = M_XYZToWgts*stimulus_XYZ(:,s);
    stimulus_spd(:,s) = B*stimulus_Wgts(:,s);   %Watts/[sr-m2-wlband]
    stimulus_XYZCheck(:,s) = T_xyz*stimulus_spd(:,s);
    stimulus_xyYCheck(:,s) = XYZToxyY(stimulus_XYZCheck(:,s)); 
    fprintf('Stimulus %d, desired xyY = %0.3f,%0.3f,%0.1f; got xyY = %0.3f,%0.3f,%0.1f\n',s,...
        stimulus_xyY(1,s),stimulus_xyY(2,s),stimulus_xyY(3,s),...
        stimulus_xyYCheck(1,s),stimulus_xyYCheck(2,s),stimulus_xyYCheck(3,s));
    plot(SToWls(S),stimulus_spd(:,s),'k');
end
xlabel('Wavelength (nm)');
ylabel('Radiance (Watts/[m2-sr-wlband])');
title('Energy');

% For input comparison's sake
stimulus_quantal=EnergyToQuanta(S,stimulus_spd/S(2)); %phot/s/sr/m2/nm

%% Convert radiance to retinal irradiance
% This is my version of the calculation.  
% It would be nice if it gives the same answer
% as vset.  This is now in Watts/[um^2-wlband].
pupilDiameterMM = 3;
pupilAreaMM2 = pi*(pupilDiameterMM/2)^2;
eyeLengthMM = EyeLength('Human');
stimulus_irradiance = RadianceToRetIrradiance(stimulus_spd,S,pupilAreaMM2,eyeLengthMM);

%% Convert from Watts to photons-sec, so that we have
stimulus_photons = EnergyToQuanta(S,stimulus_irradiance);   % retinal irradiance in quanta/[um^2-wli-sec]
stimulus_photons = stimulus_photons*1e12/S(2);                %quanta/s/m2/nm

photonsFig = figure; clf; hold on 
plot(SToWls(S),stimulus_photons,'k');
xlabel('Wavelength (nm)');
ylabel('Radiance (quanta/s/m2/nm)');
title('Quanta');

%% Convert from retinal irradiance to trolands
% First we have to convert back to energy units,
% then do the trolands conversion. 
stimulus_irradiance1 = QuantaToEnergy(S,stimulus_photons);  % returns in watts/m2/nm
trolands = RetIrradianceToTrolands(S(2)*stimulus_irradiance1/1e12,S,'Photopic','Human');
trolandsCheck = stimulus_Y*pupilAreaMM2;
for s = 1:nStimuli
    fprintf('Stimulus %d, trolands from irradiance %0.3g, direct check %0.3g\n',s,trolands(s),trolandsCheck(s));
end

% Put it all back into matrix format
intensity(1,:)=trolands(1:5);
intensity(2,:)=trolands(6:10);
intensity(3,:)=trolands(11:15);
intensity(4,:)=trolands(16:20);
intensity(5,:)=trolands(21:25);

    