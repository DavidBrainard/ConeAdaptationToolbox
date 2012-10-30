function []=plotValetonDataComparison(model,whichParams)
% function []=plotValetonDataComparison(model,whichParams)
% Plot comparison between Valeton data, van Hateren and matlab model
% Valeton data taken from figure 2 in Valeton 1983 using digitizer
% Pretty sure the data is values of the outer segment current (see bottom
% of pg 1540)
%
% All data is stored in 'ValetonData'
% Within that file, there are sets of data for different stimulus illuminances.
%   A=10^6.5 td
%   B=10^5.5 td
%   C=10^4.5 td
%   D=10^4.0 td (same as background)
%   E=10^3.5 td
%   F=10^2.5 td

% If it's working, everything should overlap.
% whichParams has options: 
%   'generic' - values used to produce figure 6 of van Hateren 2005
%   'fig7' - values used to produce figure 7 of van Hateren 2005
% model currently only has option 'vanHateren'
%
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus starts
%   timestruct.timeend - last passed time

% Close
close all;

% Set time values
timestruct.timestart=0;
timestruct.timeend=300;
timestruct.timeon=25;
timestruct.dt=.1;


% Finds contrast values algebraically based on info given in graph in
% Valeton paper
conA=(10^6.5)/(10^4)-1;
conB=(10^5.5)/(10^4)-1;
conC=(10^4.5)/(10^4)-1;
conE=(10^3.5)/(10^4)-1;
conF=(10^2.5)/(10^4)-1;

background=10^4;
pulseduration=150;

        % Gets our model's numbers.  For comparison with Valeton data, we want
        % the outer segment current.  But we want to include the feedback terms
        % on this, so we have a subfunction that runs the whole model and returns
        % the outer segment current.
        AIos=getIos(model,background,conA,pulseduration,whichParams,timestruct);
        BIos=getIos(model,background,conB,pulseduration,whichParams,timestruct);
        CIos=getIos(model,background,conC,pulseduration,whichParams,timestruct);
        DIos=getIos(model,background,0,pulseduration,whichParams,timestruct);
        EIos=getIos(model,background,conE,pulseduration,whichParams,timestruct);
        FIos=getIos(model,background,conF,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatAIos=vanHatModel(background,conA,pulseduration,whichParams,timestruct);
        vanHatBIos=vanHatModel(background,conB,pulseduration,whichParams,timestruct);
        vanHatCIos=vanHatModel(background,conC,pulseduration,whichParams,timestruct);
        vanHatDIos=vanHatModel(background,0,pulseduration,whichParams,timestruct);
        vanHatEIos=vanHatModel(background,conE,pulseduration,whichParams,timestruct);
        vanHatFIos=vanHatModel(background,conF,pulseduration,whichParams,timestruct);
        
        % Plot the stimuli for cases A-F.  These happen to live
        % in the first cell of the returned vanHatModel output.
        figure
        hold on
        plot(vanHatAIos{1})
        plot(vanHatBIos{1})
        plot(vanHatCIos{1})
        plot(vanHatDIos{1})
        plot(vanHatEIos{1})
        plot(vanHatFIos{1})
        hold off
        
        % Get Valeton Data 1983
        ValetonData;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;
        
        % Make subplots for each case.  Data are in blue, model
        % in green/red.  The model should agree with itself computed
        % in two ways, and does.  The red may look black when it overlays
        % the green.
        figure; 
        title('Relative response to background of 10^4 Td')
        subplot(3,2,1); hold on
        plot(t,AIos,'g','LineWidth',3)
        plot(t,vanHatAIos{6}-vanHatAIos{6}(length(vanHatAIos{6})),'r')
        plot(A(:,1),-A(:,2)+.05)
        title('V_h Output of stimulus=10^6.5')
        xlim([timestruct.timestart timestruct.timeend])
        
        subplot(3,2,2); hold on
        plot(t,BIos,'g','LineWidth',3)
        plot(t,vanHatBIos{6}-vanHatBIos{6}(length(vanHatBIos{6})),'r')
        plot(B(:,1),-B(:,2)+.05)
        title('V_h Output of stimulus=10^5.5')
        xlim([timestruct.timestart timestruct.timeend])
        
        subplot(3,2,3); hold on
        plot(t,CIos,'g','LineWidth',3)
        plot(t,vanHatCIos{6}-vanHatCIos{6}(length(vanHatCIos{6})),'r')
        plot(C(:,1),-C(:,2)+.05)
        title('V_h Output of stimulus=10^4.5')
        xlim([timestruct.timestart timestruct.timeend])
        
        subplot(3,2,4); hold on
        plot(t,DIos,'g','LineWidth',3)
        plot(t,vanHatDIos{6}-vanHatDIos{6}(length(vanHatDIos{6})),'r')
        plot(D(:,1),-D(:,2)+.05)
        title('V_h Output of stimulus=background=10^4')
        xlim([timestruct.timestart timestruct.timeend])
        
        subplot(3,2,5); hold on
        plot(t,EIos,'g','LineWidth',3)       
        plot(t,vanHatEIos{6}-vanHatEIos{6}(length(vanHatEIos{6})),'r')        
        plot(E(:,1),-E(:,2)+.05)
        title('V_h Output of stimulus=10^3.5')
        xlim([timestruct.timestart timestruct.timeend])

        subplot(3,2,6); hold on
        plot(t,FIos,'g','LineWidth',3)       
        plot(t,vanHatFIos{6}-vanHatFIos{6}(length(vanHatFIos{6})),'r')        
        plot(F(:,1),-F(:,2)+.05)
        title('V_h Output of stimulus=10^2.5')
        xlim([timestruct.timestart timestruct.timeend])

               
        legend('Our Model', 'Fortran Transliterated', 'Valeton Data')
        hold off
        
        
function Ios=getIos(model,background,contrast,pulseduration,whichParams,timestruct)
% function Ios=getIos(model,background,contrast,pulseduration,whichParams,timestruct)
% whichParams has options: 
%   'generic' - values used to produce figure 6 of van Hateren 2005
%   'fig7' - values used to produce figure 7 of van Hateren 2005
%
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus starts
%   timestruct.timeend - last passed time
% Runs cAdapt and outputs outer cone segment data for a given background
% illuminance and contrast level


%% Create cAdapt
t = timestruct.timestart:timestruct.dt:timestruct.timeend;   %timebase in ms   
I=background*ones(1,length(t));   %input intensity in td                       
I(t >= timestruct.timeon & t < timestruct.timeon+pulseduration) = (contrast+1)*background;
cAdapt=cAdaptCreate(model,whichParams,'timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Gets numbers
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt,'arma');
cAdapt=calcOutputHateren(cAdapt,'arma');

% Subtract off steady state to make data given in Valeton paper
Ios=cAdapt.Ios-cAdapt.X(1);