% v_cAdaptValeton
%
% Compare Valeton model to van Hateren output and to data from the same
% Valeton van Norren 1983 paper.
%
% Valeton paper uses stimulus length of 150ms: this must be the pulse
% duration for the van Hateren call.

%% Set required parameters 
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus start
%   timestruct.timeend - last passed time
% 

%%%%%%%%%%%%%%%%%%
clear; close all
%%%%%%%%%%%%%%%%%%

timestruct.dt=.1;          
timestruct.timestart=0;
timestruct.timeon=25;
timestruct.timeend=300;

% Background illuminance in td, pulse duration in ms
background=800;
contrast=2;
duration=150;
intensity=(contrast+1)*background;

%% Get cAdapt
t=timestruct.timestart:timestruct.dt:timestruct.timeend;   %ms
I=background*ones(1,length(t));   %Td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = (contrast+1)*background;
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Run calculations
% van Hateren Model
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt,'arma');

% Use ohm's law to get voltage from outer segment current and
% force steady state to zero because we're only looking at amplitude
OSvoltage=(cAdapt.Ios./cAdapt.gi)-(cAdapt.Ios(1)./cAdapt.gi(1));

% Valeton model
valetonVoltage=calcOutputValeton(background, intensity, timestruct);

% Convert Valeton to Hateren units
% We don't know what the conversion between units is because of unclarity
% in the paper, so we find a conversion factor at one point and then apply
% it to all other points.
k=findConversionValetonToHateren(background,intensity);
valetonVoltageConverted=-valetonVoltage*k;  % Negative because Valeton gives amplitude, which is abs(actual signal)

%% Plot

figure
hold on
plot(cAdapt.timebase,OSvoltage)  
plot(cAdapt.timebase,valetonVoltageConverted,'r')           
hold off