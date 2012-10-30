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

% Background illuminance in td, pulse duration in ms
background=100;
contrast=2;
duration=100;

%% Get the data from the Fortran code
% See vanHaternForMat for constant values and initializations.
ValTable=vanHatModel(background,contrast,duration,'generic',timestruct);

subplot(3,2,1); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{4}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,2); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, log(ValTable{5}), 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,3); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{8}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,4); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{6}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,5); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{6}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,6); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{7}, 'g', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
%% Get cAdapt
t=timestruct.timestart:timestruct.dt:timestruct.timeend;   %ms
I=background*ones(1,length(t));   %Td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = (contrast+1)*background;

% Create the cAdapt object.  This sets parameters and
% lets us pass to it the stimulus information.
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Run calculations and plot
%Iterative method should closely match other two methods
cAdapt=calcEstarHateren(cAdapt,'iterative');
cAdapt=calcXHateren(cAdapt,'iterative');

figure
hold on
subplot(2,1,1), plot(cAdapt.timebase,cAdapt.Estar,'b')
xlabel('time (ms)')
ylabel('Amount Activated PDE')
xlim([timestruct.timestart,timestruct.timeend])
subplot(2,1,2), plot(cAdapt.timebase,cAdapt.stimulus, 'b')
xlabel('time (ms)')
ylabel('Stimlus intensity')
xlim([timestruct.timestart,timestruct.timeend])




figure(1)
subplot(3,2,1); hold on
plot(cAdapt.timebase,cAdapt.B, 'LineWidth', 2.5)
title('B')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,2); hold on
plot(cAdapt.timebase,log(cAdapt.Tx), 'LineWidth', 2.5)
title('1/B')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,3); hold on
plot(cAdapt.timebase,cAdapt.aoverb, 'LineWidth', 2.5)
title('alpha/beta')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,4); hold on
plot(cAdapt.timebase,cAdapt.X, 'LineWidth', 2.5)
title('X')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,5); hold on
plot(cAdapt.timebase,cAdapt.Ios, 'LineWidth', 2.5)
title('Ios')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,6); hold on
plot(cAdapt.timebase,cAdapt.C, 'LineWidth', 2.5)
title('C')
xlim([timestruct.timestart,timestruct.timeend])


%Run arma method, plot overtop
%Arma should exactly match van Hateren output
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');


subplot(3,2,1), plot(cAdapt.timebase,cAdapt.B,'r', 'LineWidth',1)
title('B')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,2), plot(cAdapt.timebase,log(cAdapt.Tx),'r', 'LineWidth',1)
title('1/B')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,3), plot(cAdapt.timebase,cAdapt.aoverb,'r', 'LineWidth',1)
title('alpha/beta')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,4), plot(cAdapt.timebase,cAdapt.X,'r', 'LineWidth',1)
title('X')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,5), plot(cAdapt.timebase,cAdapt.Ios,'r', 'LineWidth',1)
title('Ios')
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,2,6), plot(cAdapt.timebase,cAdapt.C,'r', 'LineWidth',1)
title('C')
xlim([timestruct.timestart,timestruct.timeend])
hold off