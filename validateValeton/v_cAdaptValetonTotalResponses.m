% v_cAdaptValetonNormalizedResponses
% Plot total responses from calcOutputValeton against total
% responses given in figure 6 of Valeton van Norren 1983 to check accuracy
% of calculations
% Everything should be sigmmoid curves

%%%%%%%%%%%%%%
clear; close all
%%%%%%%%%%%%%%

% set up timestruct
timestruct.dt=.1;          
timestruct.timestart=0;
timestruct.timeon=25;
timestruct.timeend=300;

% set up inputs
background=10^2;    % Digitize data is for dark and 10^2, must keep it this way
pulseduration=150;  % Unchangeable, Valeton paper used only 150ms stimuli

% initialize output data
valetonData100=[];
valetonAbcissa100=[];
valetonDataDark=[];
valetonAbcissaDark=[];

%% Valeton calculations
% Valeton model at 100 td background
for k=logspace(1,7,200)
    [Ios,V]=calcOutputValeton(background,k,timestruct);      
    valetonData100=[valetonData100 V];                               
    valetonAbcissa100=[valetonAbcissa100, k];                      
end    

% Valeton model in dark
for h=logspace(1,7,200)
    [Ios,V]=calcOutputValeton(0,h,timestruct);      
    valetonDataDark=[valetonDataDark V];                               
    valetonAbcissaDark=[valetonAbcissaDark, h];                      
end

% We now call a script that will find a conversion constant that makes the 
% Hateren and Valeton graphs match up at a single point. Because we are
% unsure of the units conversion between the two papers, by artificially
% matching the units at one point we can see how well the rest of the data
% matches up without having to compare different units.
c=findConversionValetonToHateren(background,1000);
cDark=findConversionValetonToHateren(0,1000);
valetonData100Converted=valetonData100*c;
valetonDataDarkConverted=valetonDataDark*cDark;

%% Hateren Calculations
% initialize output data
haterenData100=[];
haterenAbcissa100=[];
haterenDataDark=[];
haterenAbcissaDark=[];
duration=150;

% Hateren model at 100td background
for j=logspace(1,7,200)
    t=timestruct.timestart:timestruct.dt:timestruct.timeend;   %ms
    I=background*ones(1,length(t));   %Td
    I(t >= timestruct.timeon & t < timestruct.timeon+duration) = j;
    cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);
    cAdapt=calcEstarHateren(cAdapt,'arma');
    cAdapt=calcXHateren(cAdapt,'arma');
    cAdapt=calcVisHateren(cAdapt,'arma');
    % Finds greatest difference magnitude between stimulated and steady
    % states
    maxResponse=max(abs((cAdapt.Ios./cAdapt.gi)-(cAdapt.Ios(1)./cAdapt.gi(1)))); %Current/conductance
    % Because some stimuli are smaller than the set background, the outputs
    % for those stimuli should be decrements, not increments. Since the
    % absolute value was taken above, here we correct for that sign
    % inconsistency.
    if j < background
        maxResponse=-maxResponse;
    end
    haterenData100=[haterenData100 maxResponse];                                  
    haterenAbcissa100=[haterenAbcissa100, j];
end    


% Hateren model dark-adapted
% Background = 0 for this state
for m=logspace(1,7,200)
    t=timestruct.timestart:timestruct.dt:timestruct.timeend;   %ms
    I=zeros(1,length(t));   %Td
    I(t >= timestruct.timeon & t < timestruct.timeon+duration) = m;
    cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',0);
    cAdapt=calcEstarHateren(cAdapt,'arma');
    cAdapt=calcXHateren(cAdapt,'arma');
    cAdapt=calcVisHateren(cAdapt,'arma');
    % Finds greatest difference magnitude between stimulated and steady
    % states
    maxResponse=max(abs((cAdapt.Ios./cAdapt.gi)-(cAdapt.Ios(1)./cAdapt.gi(1)))); %Current/conductance 
    haterenDataDark=[haterenDataDark maxResponse];                                 
    haterenAbcissaDark=[haterenAbcissaDark, m];
end 

 


%% Import digitizeIt data

% Gets data
TotalResponseDataValeton;

% Splits matrix into to columns of x and y coordinates
valetonRawData100=background_e2(:,2);
valetonRawAbcissa100=background_e2(:,1);
valetonRawDataDark=background_dark(:,2);
valetonRawAbcissaDark=background_dark(:,1);

% Converts data to match units as in valeton calculations above
valetonRawData100Converted=valetonRawData100*c;
valetonRawDataDarkConverted=valetonRawDataDark*cDark;

%% Plot
subplot(2,1,1)
semilogx(valetonAbcissa100,valetonData100Converted)
hold on
semilogx(valetonRawAbcissa100, valetonRawData100Converted,'r')
semilogx(haterenAbcissa100,haterenData100,'g')
semilogx(2000*ones(1,120),-20:99, 'k')
legend('Valeton model data','Valeton raw data digitized', 'van Hateren model data','van Hateren Accuracy cutoff')
xlabel('Total intesity, td')
ylabel('Fudged outer segment voltage response')
hold off
subplot(2,1,2)
semilogx(valetonAbcissaDark,valetonDataDarkConverted)
hold on
semilogx(valetonRawAbcissaDark, valetonRawDataDarkConverted,'r')
semilogx(haterenAbcissaDark,haterenDataDark,'g')
semilogx(2000*ones(1,100),1:100, 'k')
legend('Valeton model data','Valeton raw data digitized', 'van Hateren model data','van Hateren Accuracy cutoff')
xlabel('Total intesity, td')
ylabel('Fudged outer segment voltage response')
hold off

