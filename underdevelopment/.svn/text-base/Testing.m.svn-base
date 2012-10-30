% v_cAdaptIosHateren
%
% Part 2 of adaptation model: finding total photocurrent
% Using generic parameters, as given on page 339 of
%   vanHateren, 2009
% This should produce boxes 3, 4, 6,& 7 of Figure 6B page 336 of van
% Hateren 2005. All colors should overlap when it is working.

%% Initialize
clear; close all;

%% Set required parameters 
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus start
%   timestruct.timeend - last passed time
% 
timestruct.dt=.1;          
timestruct.timestart=0;
timestruct.timeon=25;
timestruct.timeend=300;
maxVs=[];
for contrast=1:200:10000
% Background illuminance in td, pulse duration in ms
background=100;

duration=150;

%% Get the data from the Fortran code
% See vanHaternForMat for constant values and initializations.
% ValTable=vanHatModel(background,contrast,duration,'generic',timestruct);
% 
% subplot(3,2,1); hold on
% plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{4}, 'g', 'LineWidth',3.5)
% xlim([timestruct.timestart,timestruct.timeend])
% subplot(3,2,2); hold on
% plot(timestruct.timestart:timestruct.dt:timestruct.timeend, log(ValTable{5}), 'g', 'LineWidth',3.5)
% xlim([timestruct.timestart,timestruct.timeend])
% subplot(3,2,3); hold on
% plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{8}, 'g', 'LineWidth',3.5)
% xlim([timestruct.timestart,timestruct.timeend])
% subplot(3,2,4); hold on
% plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{6}, 'g', 'LineWidth',3.5)
% xlim([timestruct.timestart,timestruct.timeend])
% subplot(3,2,5); hold on
% plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{6}, 'g', 'LineWidth',3.5)
% xlim([timestruct.timestart,timestruct.timeend])
% subplot(3,2,6); hold on
% plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{7}, 'g', 'LineWidth',3.5)
% xlim([timestruct.timestart,timestruct.timeend])
%% Get cAdapt
t=timestruct.timestart:timestruct.dt:timestruct.timeend;   %ms
I=background*ones(1,length(t));   %Td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = (contrast+1)*background;

% Create the cAdapt object.  This sets parameters and
% lets us pass to it the stimulus information.
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);


cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt, 'arma');
cAdapt=calcOutputHateren(cAdapt,'arma');
maxVs=[maxVs max(abs(cAdapt.vs-cAdapt.vs(1)))];
end

plot(1:20000:1000000,maxVs)