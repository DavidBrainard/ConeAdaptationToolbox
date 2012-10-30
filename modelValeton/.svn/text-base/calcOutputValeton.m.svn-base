function [VectorVoltage,Voltage]=calcOutputValeton(background, intensity, timestruct)
% function Ios=calcOutputValeton(background, intensity, timestruct)
%
% Uses equations given on page 1540 of Valeton and van Norren 1983 to find
% the outer segment response to a given increment or decrement of intensity
% on a certain background.
%
% Assumes a stimulus duration of 150 ms
%
% Background and intensity in td
%
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus start
%   timestruct.timeend - last passed time
%
% We hope that the amplitude of the output of this function will match the
% amplitude of the output given by the van Hateren model. 
%
% I - intensity of flash
% V_im - maximum increment response
% V_bm - maximum decrement response
% V_m - maximum possible response
% sigma - half saturation intensity
% n - exponent usually between .7 and 1




I=double(intensity); 
[V_im, V_dm]=getImDm(background);
n= .74;   %between .7 and 1
V_m=V_im+abs(V_dm);
sigma=getSigma(background);
Voltage=V_m*(I^n/(I^n+sigma^n));

% A vector of zeros, with the calculated voltage value repeating from the
% time on to time off
VectorVoltage=zeros(1,length(timestruct.timestart:timestruct.dt:timestruct.timeend));
VectorVoltage((timestruct.timeon)/timestruct.dt:(timestruct.timeon+150)/timestruct.dt)=Voltage;


function [im, dm]=getImDm(background)
% function [im, dm]=getImDm(background)
% Based on figure 4 on page 1543 of Valeton and van Norren 1983
% im/dm vs background data was given as a graph in the paper, so we digitized
% it and now look through that background data to find the most relevant
% im/dm data.
% Finds the point of digitized data (1st column) that is closest to the 
% given background and returns the data in the 2nd column as im/dm

imdmData;

imlist=V_im;
while size(imlist,1)>1
    if background < 10^(imlist(round(size(imlist,1)/2),1)) 
        imlist=imlist(1:round(size(imlist,1)/2),:);
    else
        imlist=imlist(ceil(size(imlist,1)/2)+1:size(imlist,1),:);
    end
end    
im=imlist(1,2);

dmlist=V_dm;
while size(dmlist)>1
    if background < 10^(dmlist(round(size(dmlist,1)/2),1)) 
        dmlist=dmlist(1:round(size(dmlist,1)/2),:);
    else
        dmlist=dmlist(ceil(size(dmlist,1)/2)+1:size(dmlist,1),:);
    end
end    
dm=dmlist(1,2);

function sigma=getSigma(background)
% function sigma=getSigma(background)
% Finds sigma by interpolation from the first two columns of table 1 on
% page 1543 of Valeton van Norren 1983
table1valuesX=[0;10^2;10^3;10^4;10^5;10^6];
table1valuesY=[10^3.2;10^3.5;10^3.9;10^4.4;10^5.2;10^6.3];
sigma=interp1(table1valuesX,table1valuesY,background);

