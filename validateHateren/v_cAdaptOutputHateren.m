% v_cAdaptOutputHateren
%
% Final part of adaptation model: finding output to bipolar cells using
% horizontal cell feedback
% Using generic parameter values, as given on page 339 of
%   vanHateren, 2009
% All colors should overlap if everything is working. This will produce the
% boxes 11, 12, 13, and 14 of Figure 6b on page 336 of van Hateren 2005, plus
% the stimulus box.

%% Initialize
clear; close all;

%% Set required parameters 
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus starts
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

%% FORTRAN calculations and plot 
% See vanHaternForMat for constant values and initializations. 
ValTable=vanHatModel(background,contrast,duration,'generic',timestruct);

subplot(3,2,1); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{1}, 'r', 'LineWidth',3)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,2); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{11}, 'r', 'LineWidth',3)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,3); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{12}, 'r', 'LineWidth',3)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,4); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{13}, 'r', 'LineWidth',3)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,5); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{14}, 'r', 'LineWidth',3)
xlim([timestruct.timestart,timestruct.timeend])

%% Get cAdapt
% Stimulus is is specified on a timebase
% t with step dt.
t =timestruct.timestart:timestruct.dt:timestruct.timeend;   %timebase in ms   
I=background*ones(1,length(t));     %input intensity in td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = (contrast+1)*background;
%get cAdapt with those specified values above
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);


%% Run calculations and plot
%%%Iterative calculations and plot
%Iterative output should closely match other two methods' output
cAdapt=calcEstarHateren(cAdapt,'iterative');
cAdapt=calcXHateren(cAdapt,'iterative');
cAdapt=calcVisHateren(cAdapt, 'iterative');
cAdapt=calcOutputHateren(cAdapt, 'iterative');

%plot results
subplot(3,2,1); hold on
plot(cAdapt.timebase, cAdapt.stimulus,'k', 'LineWidth', 2)
xlim([timestruct.timestart,timestruct.timeend])
title('Stimulus')
subplot(3,2,2); hold on
plot(cAdapt.timebase, cAdapt.vs,'k', 'LineWidth', 2)
xlim([timestruct.timestart,timestruct.timeend])
title('V_s')
subplot(3,2,3); hold on
plot(cAdapt.timebase, cAdapt.It,'k', 'LineWidth', 2)
xlim([timestruct.timestart,timestruct.timeend])
title('I_t')
subplot(3,2,4); hold on
plot(cAdapt.timebase, cAdapt.vb,'k', 'LineWidth', 2)
xlim([timestruct.timestart,timestruct.timeend])
title('V_b')
subplot(3,2,5); hold on
plot(cAdapt.timebase, cAdapt.vh,'k', 'LineWidth', 2)
xlim([timestruct.timestart,timestruct.timeend])
title('V_h')

%%% ARMA calculations and plot overtop
%ARMA output should exactly match van Hateren output.
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt, 'arma');
cAdapt=calcOutputHateren(cAdapt, 'arma');

subplot(3,2,1);  
plot(cAdapt.timebase, cAdapt.stimulus, 'c', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,2);
plot(cAdapt.timebase, cAdapt.vs, 'c', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,3); 
plot(cAdapt.timebase, cAdapt.It, 'c', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,4);
plot(cAdapt.timebase, cAdapt.vb, 'c', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,5);
plot(cAdapt.timebase, cAdapt.vh, 'c', 'LineWidth',1)
xlim([timestruct.timestart,timestruct.timeend])


hold off