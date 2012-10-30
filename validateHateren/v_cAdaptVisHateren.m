% v_cAdaptVisHateren
%
% Part 3 of adaptation model: finding receptor potential
% Using constants stated in cAdaptCreate.m, as given on page 339 of
%   vanHateren, 2009
%
% Should produce boxes 9 and 10 of Figure 6b on page 336 of van Hateren
% 2005. All colors should overlap when it is working correctly. 

% Initialize
clear; close all

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

% Background illuminance in td, pulse duration in ms
background=100;
contrast=2;
duration=100;

%% Get the data from the Fortran code
% See vanHaternForMat for constant values and initializations.
ValTable=vanHatModel(background,contrast,duration,'generic',timestruct);

subplot(3,1,1); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{1}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,1,2); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{9}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,1,3); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{10}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
%% Get cAdapt
t=timestruct.timestart:timestruct.dt:timestruct.timeend;    %ms
I=background*ones(1,length(t));   %input in Td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = (contrast+1)*background;
%get cAdapt with those specified values above
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Run calculations and plot
%%%Iterative calculations and plot
%Iterative output should closely match other two methods' output
cAdapt=calcEstarHateren(cAdapt,'iterative');
cAdapt=calcXHateren(cAdapt,'iterative');
cAdapt=calcVisHateren(cAdapt, 'iterative');

subplot(3,1,1); hold on
plot(cAdapt.timebase, cAdapt.stimulus,'b', 'LineWidth', 2.5)
xlim([timestruct.timestart,timestruct.timeend])
title('Stimulus')
subplot(3,1,2); hold on
plot(cAdapt.timebase, cAdapt.vis,'b', 'LineWidth', 2.5)
xlim([timestruct.timestart,timestruct.timeend])
title('V_i_s')
subplot(3,1,3); hold on
plot(cAdapt.timebase, cAdapt.gi,'b', 'LineWidth', 2.5)
xlim([timestruct.timestart,timestruct.timeend])
title('g_i')

%%% ARMA calculations and plot overtop
%ARMA output should exactly match van Hateren output.
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt, 'arma');

subplot(3,1,1);
plot(cAdapt.timebase, cAdapt.stimulus, 'r', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,1,2); 
plot(cAdapt.timebase, cAdapt.vis, 'r', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,1,3);
plot(cAdapt.timebase, cAdapt.gi, 'r', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
hold off