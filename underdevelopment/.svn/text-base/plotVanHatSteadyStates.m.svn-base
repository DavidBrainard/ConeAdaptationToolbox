function steadystate=plotVanHatSteadyStates(x)
% function plotVanHatSteadyStates(x)
% x is a vector of values of input in trolands, can be any length
% input should be the increase in trolands over the background, so total
% trolands in-background trolands
% Single Input Steady State Plots

model= 'vanHateren';
whichParams= 'fig7';

timestruct.timestart=0;
timestruct.timeend=300;
timestruct.timeon=25;
timestruct.dt=.1;
t = timestruct.timestart:timestruct.dt:timestruct.timeend;   %timebase in ms   

background=100; %td
pulseduration=100;  %ms

%preallocate for speed
steadystate=zeros(1,length(x));
newsteady=zeros(1,length(x));

for k=1:length(x)
contrast=x(k)/background;         %contrast  
I=background*ones(1,length(t));   %input intensity in td                       
I(t >= timestruct.timeon & t < timestruct.timeon+pulseduration) = (contrast+1)*background;  %set stimulus
%get cAdapt with those specified values above
cAdapt=cAdaptCreate(model,whichParams,'timebase',t,'stimulus',I,'dt',timestruct.dt,'background',background);

%run calculations
cAdapt=calcEstarHateren(cAdapt,'arma');
cAdapt=calcXHateren(cAdapt,'arma');
cAdapt=calcVisHateren(cAdapt, 'arma');
cAdapt=calcOutputHateren(cAdapt, 'arma');
%save steady state
steadystate(k)=cAdapt.steadyH;
end

%plot
[x, ind]=sort(x);
for j=1:length(x)
    newsteady(j)=steadystate(ind(j));
end
steadystate=newsteady;
plot(x,steadystate)