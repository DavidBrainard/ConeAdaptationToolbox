function [valetonOutput]=getOSRespCheckerboard(intensity)
% function [valetonOutput]=getOSRespCheckerboard(intensity)
% intensity is a matrix of troland values, where each value stimulates one
%   photoreceptor.
% Script that will call the valeton calculation on each pixel of the
% optical image created by ISET. Currently assuming each pixel==one
% photoreceptor. Lots of things to be worked out:
%   -Still lacking spatial element, and only returns output from outer
%       segment
%   -What background means in this context: Currently assumed that the
%       observer is adapted to the maximum intensity of the checkerboard



% Set up parameters
background=max(max(intensity));

timestruct.dt=.1;          
timestruct.timestart=0;
timestruct.timeon=25;
timestruct.timeend=300;

[nRow,nCol]=size(intensity);
valetonOutput=zeros(nRow,nCol); % Preallocation
% Step through each pixel and find valeton calculation output
for row=1:nRow
    for col=1:nCol
        [vector,voltage]=calcOutputValeton(background,intensity(row,col),timestruct);
        valetonOutput(row,col)=voltage;
    end
end

%% Plot outer segment response
figure;
h=bar3(valetonOutput);  % Output cone outersegment response
% Tell handle graphics to use interpolated rather than flat shading
shading interp
colormap cool
% For each barseries, map its CData to its ZData
for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(h,'EdgeColor','k') 
end

title('Cone Outer Segment Response mV')