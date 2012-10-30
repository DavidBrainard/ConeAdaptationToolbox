function []=plotSmithDataComparison(whichFigToPlot,whichParams)
% function []=plotSmithDataComparison(whichFigToPlot,whichParams)
%
% Plot comparison between smith data, van Hateren and matlab model
% for figure A,B,C,D,E, or F 
%
% whichParams has options: 
%   'generic' - values used to produce figure 6 of van Hateren 2005
%   'fig7' - values used to produce Figure 7 of van Hateren 2005
% The 'fig7' choice is the one that will best reproduce the agreement
% of van Hateren Figure 7, and if things are working then the lines
% in each plot should mostly overlap for this choice.  See van Hateren
% Figure 7.  
%
% Data grabbed from Smith's graphs (graph letters aren't the same, sorry!)
%   A=100Td,10ms
%   B=10Td, 10ms
%   C=1Td, 10ms
%   D=100Td, 100ms
%   E=10Td, 100ms
%   F=1Td, 100ms
% Some of the plots were hard to digitize well and there are a few imperfections
% in the digitization.
%
% Example:
%   plotSmithDataComparison('A','fig7');
%
% All data for one graph (A,B,C, etc) are stored in 'Figure(A/B/C/etc)Data'
% Within that file, there are sets of data for different contrasts.
% contrast is given by number following underscore
%    eg, A_1 is the set of data from the graph of a 100td background 
%    stimulus with contrast 1 for 10 ms

% Close
close all;

% Set timebase
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus starts
%   timestruct.timeend - last passed time
timestruct.timestart=0;
timestruct.timeend=300;
timestruct.timeon=25;
timestruct.dt=.1;

% Set model (currently only option is vanHateren)
model='vanHateren';
switch (whichFigToPlot)
    case 'A'
    %% Figure A
        background=100;
        pulseduration=10;
        % Gets our model's numbers
        A1Vh=getVh(model,background,1,pulseduration,whichParams,timestruct);
        A2Vh=getVh(model,background,2,pulseduration,whichParams,timestruct);
        A4Vh=getVh(model,background,4,pulseduration,whichParams,timestruct);
        A8Vh=getVh(model,background,8,pulseduration,whichParams,timestruct);
        A16Vh=getVh(model,background,16,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatA1Vh=vanHatModel(background,1,pulseduration,whichParams,timestruct);
        vanHatA2Vh=vanHatModel(background,2,pulseduration,whichParams,timestruct);
        vanHatA4Vh=vanHatModel(background,4,pulseduration,whichParams,timestruct);
        vanHatA8Vh=vanHatModel(background,8,pulseduration,whichParams,timestruct);
        vanHatA16Vh=vanHatModel(background,16,pulseduration,whichParams,timestruct);
       
        
        % Gets Smith Data 2008
        FigureAData;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;
        
        figure; 
        subplot(5,1,1); hold on
        plot(t,A1Vh,'g','LineWidth',3)
        plot(t,vanHatA1Vh{14}-vanHatA1Vh{14}(length(vanHatA1Vh{14})),'r')
        plot(A_1(:,2),A_1(:,1),'k')
        title('V_h Output of contrast=1')
        xlim([0,300])
        ylim([-15,5])
        
        subplot(5,1,2); hold on
        plot(t,A2Vh,'g','LineWidth',3)
        plot(t,vanHatA2Vh{14}-vanHatA2Vh{14}(length(vanHatA2Vh{14})),'r')
        plot(A_2(:,2),A_2(:,1),'k')
        title('V_h Output of contrast=2')
        xlim([0,300])
        ylim([-15,5])
        
        subplot(5,1,3); hold on
        plot(t,A4Vh,'g','LineWidth',3)
        plot(t,vanHatA4Vh{14}-vanHatA4Vh{14}(length(vanHatA4Vh{14})),'r')
        plot(A_4(:,2),A_4(:,1),'k')
        title('V_h Output of contrast=4')
        xlim([0,300])
        ylim([-15,5])
        
        subplot(5,1,4); hold on
        plot(t,A8Vh,'g','LineWidth',3)
        plot(t,vanHatA8Vh{14}-vanHatA8Vh{14}(length(vanHatA8Vh{14})),'r')
        plot(A_8(:,2),A_8(:,1),'k')
        title('V_h Output of contrast=8')
        xlim([0,300])
        ylim([-15,5])
        
        subplot(5,1,5); hold on
        plot(t,A16Vh,'g','LineWidth',3)       
        plot(t,vanHatA16Vh{14}-vanHatA16Vh{14}(length(vanHatA16Vh{14})),'r')        
        plot(A_16(:,2),A_16(:,1),'k')
        title('V_h Output of contrast=16')
        xlim([0,300])
        ylim([-15,5])
        
        legend('Our Model', 'Fortran Transliterated', 'Smith Data')
        hold off
    case 'B'
        %% Figure B
        background=10;
        pulseduration=10;
        % Gets our model's numbers
        B1Vh=getVh(model,background,1,pulseduration,whichParams,timestruct);
        B2Vh=getVh(model,background,2,pulseduration,whichParams,timestruct);
        B4Vh=getVh(model,background,4,pulseduration,whichParams,timestruct);
        B8Vh=getVh(model,background,8,pulseduration,whichParams,timestruct);
        B16Vh=getVh(model,background,16,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatB1Vh=vanHatModel(background,1,pulseduration,whichParams,timestruct);
        vanHatB2Vh=vanHatModel(background,2,pulseduration,whichParams,timestruct);
        vanHatB4Vh=vanHatModel(background,4,pulseduration,whichParams,timestruct);
        vanHatB8Vh=vanHatModel(background,8,pulseduration,whichParams,timestruct);
        vanHatB16Vh=vanHatModel(background,16,pulseduration,whichParams,timestruct);
       
        
        % Gets Smith Data 2008
        FigureBData;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;

        
        figure; 
        subplot(5,1,1); hold on
        plot(t,B1Vh,'g','LineWidth',3)
        plot(t,vanHatB1Vh{14}-vanHatB1Vh{14}(length(vanHatB1Vh{14})),'r')
        plot(B_1(:,1),B_1(:,2),'k')
        title('V_h Output of contrast=1')
        xlim([0,300])
        
        subplot(5,1,2); hold on
        plot(t,B2Vh,'g','LineWidth',3)
        plot(t,vanHatB2Vh{14}-vanHatB2Vh{14}(length(vanHatB2Vh{14})),'r')
        plot(B_2(:,1),B_2(:,2),'k')
        title('V_h Output of contrast=2')
        xlim([0,300])
        
        subplot(5,1,3); hold on
        plot(t,B4Vh,'g','LineWidth',3)
        plot(t,vanHatB4Vh{14}-vanHatB4Vh{14}(length(vanHatB4Vh{14})),'r')
        plot(B_4(:,1),B_4(:,2),'k')
        title('V_h Output of contrast=4')
        xlim([0,300])
        
        subplot(5,1,4); hold on
        plot(t,B8Vh,'g','LineWidth',3)
        plot(t,vanHatB8Vh{14}-vanHatB8Vh{14}(length(vanHatB8Vh{14})),'r')
        plot(B_8(:,1),B_8(:,2),'k')
        title('V_h Output of contrast=8')
        xlim([0,300])
        
        subplot(5,1,5); hold on
        plot(t,B16Vh,'g','LineWidth',3)       
        plot(t,vanHatB16Vh{14}-vanHatB16Vh{14}(length(vanHatB16Vh{14})),'r')        
        plot(B_16(:,1),B_16(:,2),'k')
        title('V_h Output of contrast=16')
        xlim([0,300])
        
        legend('Our Model', 'Fortran Transliterated', 'Smith Data')
        hold off
    case 'C'
        %% Figure C
        % Graphs of contrasts=1 and 2 were indecipherable
        background=1;
        pulseduration=10;
        % Gets our model's numbers
        C4Vh=getVh(model,background,4,pulseduration,whichParams,timestruct);
        C8Vh=getVh(model,background,8,pulseduration,whichParams,timestruct);
        C16Vh=getVh(model,background,16,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatC4Vh=vanHatModel(background,4,pulseduration,whichParams,timestruct);
        vanHatC8Vh=vanHatModel(background,8,pulseduration,whichParams,timestruct);
        vanHatC16Vh=vanHatModel(background,16,pulseduration,whichParams,timestruct);
       
        
        % Gets Smith Data 2008
        FigureCData;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;

        
        subplot(3,1,1); hold on
        plot(t,C4Vh,'g','LineWidth',3)
        plot(t,vanHatC4Vh{14}-vanHatC4Vh{14}(length(vanHatC4Vh{14})),'r')
        plot(C_4(:,1),C_4(:,2),'k')
        title('V_h Output of contrast=4')
        xlim([0,300])
        
        subplot(3,1,2); hold on
        plot(t,C8Vh,'g','LineWidth',3)
        plot(t,vanHatC8Vh{14}-vanHatC8Vh{14}(length(vanHatC8Vh{14})),'r')
        plot(C_8(:,1),C_8(:,2),'k')
        title('V_h Output of contrast=8')
        xlim([0,300])
        
        subplot(3,1,3); hold on
        plot(t,C16Vh,'g','LineWidth',3)       
        plot(t,vanHatC16Vh{14}-vanHatC16Vh{14}(length(vanHatC16Vh{14})),'r')        
        plot(C_16(:,1),C_16(:,2),'k')
        title('V_h Output of contrast=16')
        xlim([0,300])
        
        legend('Our Model', 'Fortran Transliterated', 'Smith Data')
        hold off
    case 'D'
        %% Figure D
        background=100;
        pulseduration=100;
        % Gets our model's numbers
        D1Vh=getVh(model,background,1,pulseduration,whichParams,timestruct);
        D2Vh=getVh(model,background,2,pulseduration,whichParams,timestruct);
        D4Vh=getVh(model,background,4,pulseduration,whichParams,timestruct);
        D8Vh=getVh(model,background,8,pulseduration,whichParams,timestruct);
        D16Vh=getVh(model,background,16,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatD1Vh=vanHatModel(background,1,pulseduration,whichParams,timestruct);
        vanHatD2Vh=vanHatModel(background,2,pulseduration,whichParams,timestruct);
        vanHatD4Vh=vanHatModel(background,4,pulseduration,whichParams,timestruct);
        vanHatD8Vh=vanHatModel(background,8,pulseduration,whichParams,timestruct);
        vanHatD16Vh=vanHatModel(background,16,pulseduration,whichParams,timestruct);
       
        
        % Gets Smith Data 2008
        FigureDData;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;

        
        figure; 
        subplot(5,1,1); hold on
        plot(t,D1Vh,'g','LineWidth',3)
        plot(t,vanHatD1Vh{14}-vanHatD1Vh{14}(length(vanHatD1Vh{14})),'r')
        plot(D_1(:,1),D_1(:,2),'k')
        title('V_h Output of contrast=1')
        xlim([0,300])
        
        subplot(5,1,2); hold on
        plot(t,D2Vh,'g','LineWidth',3)
        plot(t,vanHatD2Vh{14}-vanHatD2Vh{14}(length(vanHatD2Vh{14})),'r')
        plot(D_2(:,1),D_2(:,2),'k')
        title('V_h Output of contrast=2')
        xlim([0,300])
        
        subplot(5,1,3); hold on
        plot(t,D4Vh,'g','LineWidth',3)
        plot(t,vanHatD4Vh{14}-vanHatD4Vh{14}(length(vanHatD4Vh{14})),'r')
        plot(D_4(:,1),D_4(:,2),'k')
        title('V_h Output of contrast=4')
        xlim([0,300])
        
        subplot(5,1,4); hold on
        plot(t,D8Vh,'g','LineWidth',3)
        plot(t,vanHatD8Vh{14}-vanHatD8Vh{14}(length(vanHatD8Vh{14})),'r')
        plot(D_8(:,1),D_8(:,2),'k')
        title('V_h Output of contrast=8')
        xlim([0,300])
        
        subplot(5,1,5); hold on
        plot(t,D16Vh,'g','LineWidth',3)       
        plot(t,vanHatD16Vh{14}-vanHatD16Vh{14}(length(vanHatD16Vh{14})),'r')        
        plot(D_16(:,1),D_16(:,2),'k')
        title('V_h Output of contrast=16')
        xlim([0,300])
        
        legend('Our Model', 'Fortran Transliterated', 'Smith Data')
        hold off
    case 'E'
        %% Figure E
        background=10;
        pulseduration=100;
        % Gets our model's numbers
        E1Vh=getVh(model,background,1,pulseduration,whichParams,timestruct);
        E2Vh=getVh(model,background,2,pulseduration,whichParams,timestruct);
        E4Vh=getVh(model,background,4,pulseduration,whichParams,timestruct);
        E8Vh=getVh(model,background,8,pulseduration,whichParams,timestruct);
        E16Vh=getVh(model,background,16,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatE1Vh=vanHatModel(background,1,pulseduration,whichParams,timestruct);
        vanHatE2Vh=vanHatModel(background,2,pulseduration,whichParams,timestruct);
        vanHatE4Vh=vanHatModel(background,4,pulseduration,whichParams,timestruct);
        vanHatE8Vh=vanHatModel(background,8,pulseduration,whichParams,timestruct);
        vanHatE16Vh=vanHatModel(background,16,pulseduration,whichParams,timestruct);
       
        
        % Gets Smith Eata 2008
        FigureEdata;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;

        
        figure; 
        subplot(5,1,1); hold on
        plot(t,E1Vh,'g','LineWidth',3)
        plot(t,vanHatE1Vh{14}-vanHatE1Vh{14}(length(vanHatE1Vh{14})),'r')
        plot(E_1(:,1),E_1(:,2),'k')
        title('V_h Output of contrast=1')
        xlim([0,300])
        
        subplot(5,1,2); hold on
        plot(t,E2Vh,'g','LineWidth',3)
        plot(t,vanHatE2Vh{14}-vanHatE2Vh{14}(length(vanHatE2Vh{14})),'r')
        plot(E_2(:,1),E_2(:,2),'k')
        title('V_h Output of contrast=2')
        xlim([0,300])
        
        subplot(5,1,3); hold on
        plot(t,E4Vh,'g','LineWidth',3)
        plot(t,vanHatE4Vh{14}-vanHatE4Vh{14}(length(vanHatE4Vh{14})),'r')
        plot(E_4(:,1),E_4(:,2),'k')
        title('V_h Output of contrast=4')
        xlim([0,300])
        
        subplot(5,1,4); hold on
        plot(t,E8Vh,'g','LineWidth',3)
        plot(t,vanHatE8Vh{14}-vanHatE8Vh{14}(length(vanHatE8Vh{14})),'r')
        plot(E_8(:,1),E_8(:,2),'k')
        title('V_h Output of contrast=8')
        xlim([0,300])
        
        subplot(5,1,5); hold on
        plot(t,E16Vh,'g','LineWidth',3)       
        plot(t,vanHatE16Vh{14}-vanHatE16Vh{14}(length(vanHatE16Vh{14})),'r')        
        plot(E_16(:,1),E_16(:,2),'k')
        title('V_h Output of contrast=16')
        xlim([0,300])
        
        legend('Our Model', 'Fortran Transliterated', 'Smith Eata')
        hold off
    case 'F'
        %% Figure F
        background=1;
        pulseduration=100;
        % Gets our model's numbers
        F1Vh=getVh(model,background,1,pulseduration,whichParams,timestruct);
        F2Vh=getVh(model,background,2,pulseduration,whichParams,timestruct);
        F4Vh=getVh(model,background,4,pulseduration,whichParams,timestruct);
        F8Vh=getVh(model,background,8,pulseduration,whichParams,timestruct);
        F16Vh=getVh(model,background,16,pulseduration,whichParams,timestruct);

        % Gets van Hateren model's numbers (should be same as ours)
        vanHatF1Vh=vanHatModel(background,1,pulseduration,whichParams,timestruct);
        vanHatF2Vh=vanHatModel(background,2,pulseduration,whichParams,timestruct);
        vanHatF4Vh=vanHatModel(background,4,pulseduration,whichParams,timestruct);
        vanHatF8Vh=vanHatModel(background,8,pulseduration,whichParams,timestruct);
        vanHatF16Vh=vanHatModel(background,16,pulseduration,whichParams,timestruct);
       
        
        % Gets Smith Fata 2008
        FigureFData;
        t=timestruct.timestart:timestruct.dt:timestruct.timeend;

        
        figure; 
        subplot(5,1,1); hold on
        plot(t,F1Vh,'g','LineWidth',3)
        plot(t,vanHatF1Vh{14}-vanHatF1Vh{14}(length(vanHatF1Vh{14})),'r')
        plot(F_1(:,1),F_1(:,2),'k')
        title('V_h Output of contrast=1')
        xlim([0,300])
        
        subplot(5,1,2); hold on
        plot(t,F2Vh,'g','LineWidth',3)
        plot(t,vanHatF2Vh{14}-vanHatF2Vh{14}(length(vanHatF2Vh{14})),'r')
        plot(F_2(:,1),F_2(:,2),'k')
        title('V_h Output of contrast=2')
        xlim([0,300])
        
        subplot(5,1,3); hold on
        plot(t,F4Vh,'g','LineWidth',3)
        plot(t,vanHatF4Vh{14}-vanHatF4Vh{14}(length(vanHatF4Vh{14})),'r')
        plot(F_4(:,1),F_4(:,2),'k')
        title('V_h Output of contrast=4')
        xlim([0,300])
        
        subplot(5,1,4); hold on
        plot(t,F8Vh,'g','LineWidth',3)
        plot(t,vanHatF8Vh{14}-vanHatF8Vh{14}(length(vanHatF8Vh{14})),'r')
        plot(F_8(:,1),F_8(:,2),'k')
        title('V_h Output of contrast=8')
        xlim([0,300])
        
        subplot(5,1,5); hold on
        plot(t,F16Vh,'g','LineWidth',3)       
        plot(t,vanHatF16Vh{14}-vanHatF16Vh{14}(length(vanHatF16Vh{14})),'r')        
        plot(F_16(:,1),F_16(:,2),'k')
        title('V_h Output of contrast=16')
        xlim([0,300])
        
        legend('Our Model', 'Fortran Transliterated', 'Smith Fata')
        hold off
    otherwise 
        error('No such plot exists. Please enter the letter of the plot to display.')
end

function [vh,steadyH]=getVh(model,background,contrast,pulseduration,whichParams,timestruct)
% function [vh,steadyH]=getVh(model,background,contrast,pulseduration,whichParams,timestruct)
% whichParams has options: 
%   'generic' - values used to produce figure 6 of van Hateren 2005
%   'fig7' - values used to produce figure 7 of van Hateren 2005
%
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus starts
%   timestruct.timeend - last passed time
% Runs cAdapt and outputs horizontal cell data for a given background
% illuminance and contrast level
%
% vh returned has steady state (as found in calcOutputHateren) already subtracted off
% and delay already incorporated (by calcOutputHateren)


%% Create cAdapt
t = timestruct.timestart:timestruct.dt:timestruct.timeend;   %timebase in ms   
I=background*ones(1,length(t));   %input intensity in td                       
I(t >= timestruct.timeon & t < timestruct.timeon+pulseduration) = (contrast+1)*background;
cAdapt=cAdaptCreate(model,whichParams,'timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Gets numbers
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt, 'arma');
cAdapt=calcOutputHateren(cAdapt, 'arma');

steadyH=cAdapt.steadyH;
vh=cAdapt.vh-cAdapt.steadyH;    %Subtract steady state