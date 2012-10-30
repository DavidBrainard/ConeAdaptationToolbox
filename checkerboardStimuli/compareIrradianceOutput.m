% Compare VSET and direct calculations
% Doesn't work unless conesPerCheck is set to 1 in s_sceneCheckerboardTest
clear; close all;

% outputs are being returned as a set of numbers of photons at each pixel,
% at each wavelength band value. Since this doesn't get mushed back into
% one set of photons spanning the whole visible spectrum until later in the
% two scripts, here I take only the values at 500 nm.
[oiPhotons,intensity]=s_sceneCheckerboardTest('experiment','one','ExpCol',4,'midSquareNum',5);
[stimulus_photons,intensity]=checkerboardStimulus('experiment','one','ExpCol',4,'midSquareNum',5);
% Cut output down to a single matrix and then make it into a row vector
stimulus_photons1=stimulus_photons(11,:);
irradPhots(1,:)=stimulus_photons1(1:5);
irradPhots(2,:)=stimulus_photons1(6:10);
irradPhots(3,:)=stimulus_photons1(11:15);
irradPhots(4,:)=stimulus_photons1(16:20);
irradPhots(5,:)=stimulus_photons1(21:25);

[m,n]=size(irradPhots);
[l,p]=size(oiPhotons);
% Chop off the blurred sides and make a row vector
oiPhotons1=oiPhotons((l-m)/2+1:l-((l-m)/2),(l-m)/2+1:l-((l-m)/2),11);
oiPhotons2=[oiPhotons1(1,:) oiPhotons1(2,:) oiPhotons1(3,:) oiPhotons1(4,:) oiPhotons1(5,:)];

close all;
figure;
h=bar3(abs(irradPhots-oiPhotons1)); % Tell handle graphics to use interpolated rather than flat shading
shading interp
colormap cool
% These are bars of the error between the direct calculations and VSET.
% They should all be sort of flat except for the center square, whose
% value is randomized and need not be the same between the two scripts.
%For each barseries, map its CData to its ZData
for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(h,'EdgeColor','k') 
end
title('Error between direct calculations and VSET')

% This is a comparison between the number of irradiance photons spit out by
% the two scripts.
% This should plot more or less a line if everything is working correctly.
figure;
[stimulus_photons2,idx1]=sort(stimulus_photons1);
[oiPhotons3,idx2]=sort(oiPhotons2);
plot(stimulus_photons2,oiPhotons3)
