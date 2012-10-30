% v_cAdaptEstarHateren
%
% Validation script for first part of van Hateren adaptation model. This part finds amount of activated PDE (E*).
%
% See also cAdaptCreate, cAdaptSet, calcEstarHateren
%
% When this code is working, it should produce the first two boxees of
% figure B on page 336 of van Hateren 2005. All of the colors should
% overlap, except green, which represents the gain.
%
% 6/xx/12  ekf  Wrote it.
% 6/24/12  dhb  Review and clean.

%% Close anything left over
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

%% Set up an impulse stimulus.
%
% Stimulus is is specified on a timebase
% t with step dt.  Here the input is an impulse function
t = timestruct.timestart:timestruct.dt:timestruct.timeend;   %timebase in ms
I = zeros(1,length(t));      %input intensity in isomerizations/ms
I(t == 0) = 1;   

%% Create the cAdapt object for impulse.  This sets parameters and
% lets us pass to it the stimulus information.
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Run Estar calculations for impulse
cAdapt=calcEstarHateren(cAdapt);

%% Plot impulse responses for Rstar and Estar.
%
% We get the impulse responses in two ways.  The first way is to grab the impulse 
% responses computed inside of calcEstarHateren. These are what are used
% in the convolution computation of Rstar and Estar.  For comparison
% we use the convolution of the two impulse responses
%
% The second way is via the computation of Rstar and Estar for the impulse
% input defined above.  If everything is working, these should match the
% impulse responses used in the convolution.  Mostly this is a check that
% we understand the indexing conventions of the convolution.
%
% If everything is working, the black and red lines in the two plots should
% overlay.

% The overall impulse response of the two lowpass filters is obtained
% as the the convolution of each of the impulse responses. We compute
% that here to have something to check against in the second plot
EstarRstarImpulseResponse = conv(cAdapt.EstarImpulseResponse,cAdapt.RstarImpulseResponse);
EstarRstarImpulseResponseTimebase = cAdapt.dt*(0:length(EstarRstarImpulseResponse)-1);

% Make the plot
figure; clf; 
subplot(2,1,1); hold on
plot(cAdapt.RstarImpulseResponseTimebase,cAdapt.RstarImpulseResponse,'k','LineWidth',2);
plotTimes = find(t >= 0 & t <= max(cAdapt.RstarImpulseResponseTimebase));
plot(cAdapt.timebase(plotTimes),cAdapt.Rstar(plotTimes),'r','LineWidth',1);
xlabel('Time (ms)'); ylabel('Activations');
title('R* impulse response');
subplot(2,1,2); hold on
plot(EstarRstarImpulseResponseTimebase,EstarRstarImpulseResponse,'k','LineWidth',2);
plotTimes = find(t >= 0 & t <= max(EstarRstarImpulseResponseTimebase));
plot(cAdapt.timebase(plotTimes),cAdapt.Estar(plotTimes),'r','LineWidth',1);
xlabel('Time (ms)'); ylabel('Activations');
title('E* impulse response');

%% Generate first Figure 6 panels
% Get the data from the Fortran code transliterated into Matlab and plot.
ValTable=vanHatModel(background,contrast,duration,'generic',timestruct);

figure;
subplot(3,1,1); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{1}, 'b', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,1,2); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{2}, 'b', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])
subplot(3,1,3); hold on
plot(timestruct.timestart:timestruct.dt:timestruct.timeend, ValTable{3}, 'b', 'LineWidth',3.5)
xlim([timestruct.timestart,timestruct.timeend])

% Now with our code.
%
% Stimulus is is specified on a timebase t with step dt as defined at
% beginning of script.
t = timestruct.timestart:timestruct.dt:timestruct.timeend;   % timebase in ms
I = background*ones(1,length(t));                            % input intensity in td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = (contrast+1)*background; 

% Create the cAdapt object.  This sets parameters and
% lets us pass to it the stimulus information.
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

% Run Estar calculations and plot.  The first and third
% panels here should look like the first and second panels
% of van Hateren (2005), Figure 6.  Since there are no
% units on panel 2 of that figure, we can't quite tell
% if what we have is right.  But the time course looks
% good. The gain in the equations for Rstar and Estar is eliminated and
% taken care of in step 2.
cAdapt=calcEstarHateren(cAdapt,'iterative');

% Plot of results
subplot(3,1,1); hold on
plot(cAdapt.timebase,cAdapt.stimulus,'k','LineWidth',1.5);
xlim([timestruct.timestart, timestruct.timeend]);
xlabel('Time (ms)');
ylabel('Intensity');
title('Stimulus');

subplot(3,1,2); hold on
plot(cAdapt.timebase,cAdapt.Rstar,'k','LineWidth',1.5);
xlim([timestruct.timestart , timestruct.timeend]);
xlabel('Time (ms)')
ylabel('Activations');
title('R*');

subplot(3,1,3); hold on
plot(cAdapt.timebase,cAdapt.Estar,'k','LineWidth',1.5);
xlim([timestruct.timestart , timestruct.timeend]);
xlabel('Time (ms)')
ylabel('Activations')
title('E*');

% NOTE: Van Hateren dumps the gains Tr*Ar and Te*Ae
% out of the Rstar and Estar equation and claims that they are compensated
% for by the kb constant in computation of outer segment current.
% We may eventually want to try to put in the real gains for each step,
% but not today.
predRstarSteadyState = cAdapt.stimulus;       %*cAdapt.Ar*cAdapt.Tr;
predEstarSteadyState = predRstarSteadyState;  %*cAdapt.Ae*cAdapt.Te;
subplot(3,1,2); hold on
plot(cAdapt.timebase,predRstarSteadyState,'g','LineWidth',1);
subplot(3,1,3); hold on
plot(cAdapt.timebase,predEstarSteadyState,'g','LineWidth',1);

% Check that the times make sense
%
% We expect that 1/e of the steady state change will happen one Tr after
% the stimulus change for R*, and one (Tr + Te) after the
% stimulus change for E*.  We can mark these times on the plots
% too.
%
% This is not currently quite making sense.  I think that is because
% convolution doesn't produce the same time course to asymptote as
% simple decay.
subplot(3,1,2); hold on
minRstarSteady = min(predRstarSteadyState);
maxRstarSteady = max(predRstarSteadyState);
plot([timestruct.timestart+cAdapt.Tr timestruct.timestart+cAdapt.Tr],[minRstarSteady minRstarSteady+exp(-1)*(maxRstarSteady-minRstarSteady)],'b','LineWidth',1);
plot([timestruct.timeend+cAdapt.Tr timestruct.timeend+cAdapt.Tr],[maxRstarSteady maxRstarSteady+exp(-1)*(minRstarSteady-maxRstarSteady)],'b','LineWidth',1);

% Compute using ARMA method and add to plot.
%
% This should agree with the convolution method,
% and does.
cAdapt=calcEstarHateren(cAdapt,'arma');
subplot(3,1,2); 
plot(cAdapt.timebase,cAdapt.Rstar,'r','LineWidth',1);
subplot(3,1,3);
plot(cAdapt.timebase,cAdapt.Estar,'r','LineWidth',1);
hold off
