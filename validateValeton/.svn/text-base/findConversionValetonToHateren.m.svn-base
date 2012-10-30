function conversion=findConversionValetonToHateren(background,intensity)
% function conversion=findConversionValetonToHateren(background,intensity)
% HaterenAmplitude=(conversion)*(ValetonAmplitude)
% background illumination and intensity in td.




%% Set up input parameters
timestruct.dt=.1;          
timestruct.timestart=0;
timestruct.timeon=25;
timestruct.timeend=300;

duration=150;   % Valeton used only stimuli of 150ms, this can not be changed.


%% Get cAdapt
t=timestruct.timestart:timestruct.dt:timestruct.timeend;   %ms
I=background*ones(1,length(t));   %Td
I(t >= timestruct.timeon & t < timestruct.timeon+duration) = intensity;
cAdapt=cAdaptCreate('vanHateren','generic','timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%% Run calculations
% van Hateren model
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt,'arma');

% Use ohm's law to get voltage from outer segment current and
% force steady state to zero because we're only looking at amplitude
% This causes an error if the arguments background and intensity are
% exactly the same.
OSvoltage=(cAdapt.Ios./cAdapt.gi)-(cAdapt.Ios(1)./cAdapt.gi(1));

% Valeton model
valetonVoltage=calcOutputValeton(background, intensity, timestruct);



%% Find amplitudes and conversion factor
% HaterenAmplitude=(conversion)*(ValetonAmplitude)

HaterenAmplitude=max(abs(OSvoltage));
ValetonAmplitude=max(abs(valetonVoltage));

conversion = HaterenAmplitude/ValetonAmplitude;