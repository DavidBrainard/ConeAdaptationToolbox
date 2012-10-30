% Compare VSET and direct calculations
% If everything is working, the last graph should show a bunch of squares
% (error bars of height zero) with one tall bar in the center representing
% the middle square, whose value is randomized and so should not
% necessarily be the same between the two scripts.
%
% Doesn't work unless conesPerCheck is set to 1 in s_sceneCheckerboardTest
%
% 8/17/12   ekf     wrote it

clear all; close all;

[oiPhotons,intensity]=s_sceneCheckerboardTest('experiment','one','ExpCol',4,'midSquareNum',5);
vsetData=getOSRespCheckerboard(intensity);
[stimulus_photons,intensity]=checkerboardStimulus('experiment','one','ExpCol',4,'midSquareNum',5);
directData=getOSRespCheckerboard(intensity);
[m,n]=size(directData);
[l,p]=size(vsetData);
% get rid of blurred edges for comparison
vsetData1=vsetData((l-m)/2+1:l-((l-m)/2),(l-m)/2+1:l-((l-m)/2));


close all;
figure;
h=bar3(abs(directData-vsetData1));% Tell handle graphics to use interpolated rather than flat shading
shading interp
colormap cool
%For each barseries, map its CData to its ZData
for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(h,'EdgeColor','k') 
end

title('Error between direct calculations and VSET')