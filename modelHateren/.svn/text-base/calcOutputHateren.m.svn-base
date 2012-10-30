function cAdapt=calcOutputHateren(cAdapt,method)
% cAdapt=calcoutputHateren(cAdapt,method)
% Final step of Hateren calculations containing cone-horizontal cell
% feedback loop
% Calculates output to horizontal and bipolar cells, accounting for feedback from
% horizontal cells
% Using equations found on pages 335-336 of van Hateren 2005
%
% We have a few different numerical methods implemented
%   'iterative' - Straightforward implementation of difference equations.
%   'arma' - Smarter implementation of difference equations.  Uses a 
%            little calculus to improve numerical accuracy.
% Notes: 
%   Van Hateren had a variable called background, which we think might represent
%   the stimulus before time 0.  Or else a spatial background for some
%   future implementation.  In any case, it was always set to the first stimulus value,
%   and we've set up the code on the assumption that this will always be true.
%
% Contact ekfrance@gmail.com

%% Preallocations and initial values
% Values as given in fortran code
cAdapt.visprime=cAdapt.vis;
cAdapt.aI=(cAdapt.visprime./cAdapt.vI).^cAdapt.mu;
gripmax=cAdapt.gt./cAdapt.aI;
% uses numerical method to find steady-state value of vs
% cAdapt.vis(1) and gripmax(1) are assumed to be steady-state values of
% those variables
options = optimset('fsolve');
options = optimset(options,'Diagnostics','off','Display','off');
vs=fsolve(@(x) x-(cAdapt.vis(1)-gripmax(1)/(1+exp(-(x-cAdapt.vk)/cAdapt.vn))),-13,options);
cAdapt.vs=vs*ones(1,length(cAdapt.timebase));
cAdapt.It=(gripmax./(1+exp(-(vs-cAdapt.vk)/cAdapt.vn)));
cAdapt.It1=cAdapt.It;
cAdapt.vb=cAdapt.It;
cAdapt.vh=cAdapt.It;
%% Calculations
% Iterative equations derived from equations on pg 335-336
%
% Formulas using arma coefficients are as given at the header of table 1 in
% van Hateren 2008
switch (method)
    case 'iterative'
        % Calculate output voltage to horizontal and bipolar cells as well 
        % as horizontal cell feedback iteratively
        %
        % Equations derived from those on pgs 335-336 van Hateren
        % 2005
        % Found by solving for change in variable and adding to variable's
        % previous value
        for k=2:length(cAdapt.timebase)
        cAdapt.visprime(k)=cAdapt.visprime(k-1)+(cAdapt.dt/cAdapt.Ta)*(cAdapt.vis(k)-cAdapt.visprime(k-1));
        cAdapt.aI(k)=(cAdapt.visprime(k)/cAdapt.vI)^cAdapt.mu;
        cAdapt.vs(k)=cAdapt.vis(k)-cAdapt.vh(k-1);
        cAdapt.It(k)=(cAdapt.gt/cAdapt.aI(k))/(1+exp(-(cAdapt.vs(k)-cAdapt.vk)/cAdapt.vn));
        cAdapt.It1(k)=cAdapt.It1(k-1)+(cAdapt.dt/cAdapt.T1)*(cAdapt.It(k)-cAdapt.It1(k-1));
        cAdapt.vb(k)=cAdapt.vb(k-1)+(cAdapt.dt/(cAdapt.aI(k)*cAdapt.T2))*(cAdapt.It1(k)-cAdapt.vb(k-1));
        cAdapt.vh(k)=cAdapt.vh(k-1)+(cAdapt.dt/(cAdapt.aI(k)*cAdapt.Th))*(cAdapt.vb(k)-cAdapt.vh(k-1)); 
        end
    case 'arma'
        % Calculate output voltage to horizontal and bipolar cells as well 
        % as horizontal cell feedback using ARMA coefficients
        %
        % Calculated according to formula for First-Order Hold in table 1
        % of van Hateren 2008
        %
        %%% Calculate non-time dependent arma coeffs
        temp_T1 = cAdapt.T1/cAdapt.dt;
        f1_tau_1=exp(-1/temp_T1);   %-a1
        f2_tau_1=temp_T1-(1+temp_T1)*f1_tau_1; %b1
        f3_tau_1=1-temp_T1+temp_T1*f1_tau_1;  %b0
        
        temp_Ta = cAdapt.Ta/cAdapt.dt;
        f1_tau_a=exp(-1/temp_Ta);   %-a1
        f2_tau_a=temp_Ta-(1+temp_Ta)*f1_tau_a; %b1
        f3_tau_a=1-temp_Ta+temp_Ta*f1_tau_a;  %b0
        
        for k=2:length(cAdapt.timebase)
            %calculate aI in order to calculate two arma coeffs          
            cAdapt.visprime(k)=f1_tau_a*(cAdapt.visprime(k-1))+...
                               f2_tau_a*(cAdapt.vis(k-1))+...
                               f3_tau_a*(cAdapt.vis(k));
            cAdapt.aI(k)=(cAdapt.visprime(k)/cAdapt.vI)^cAdapt.mu;
            %calculate static arma coeffs         
            temp_T2 = cAdapt.aI(k)*cAdapt.T2/cAdapt.dt;
            f1_tau_2=exp(-1/temp_T2); %-a1
            f2_tau_2=temp_T2-(1+temp_T2)*f1_tau_2; %b1
            f3_tau_2=1-temp_T2+temp_T2*f1_tau_2;  %b0

            temp_Th = cAdapt.aI(k)*cAdapt.Th/cAdapt.dt;
            f1_tau_h=exp(-1/temp_Th); %-a1
            f2_tau_h=temp_Th-(1+temp_Th)*f1_tau_h; %b1
            f3_tau_h=1-temp_Th+temp_Th*f1_tau_h;  %b0
            %calculate transmitter current, horizontal and bipolar cell
            %output
            cAdapt.vs(k)=cAdapt.vis(k)-cAdapt.vh(k-1);
            cAdapt.It(k)=(cAdapt.gt/cAdapt.aI(k))/(1+exp(-(cAdapt.vs(k)-cAdapt.vk)/cAdapt.vn));
            cAdapt.It1(k)=f1_tau_1*cAdapt.It1(k-1)+...
                          f2_tau_1*cAdapt.It(k-1)+...
                          f3_tau_1*cAdapt.It(k);
            cAdapt.vb(k)=f1_tau_2*cAdapt.vb(k-1)+...
                         f2_tau_2*cAdapt.It1(k-1)+...
                         f3_tau_2*cAdapt.It1(k);
            cAdapt.vh(k)=f1_tau_h*cAdapt.vh(k-1)+...
                         f2_tau_h*cAdapt.vb(k-1)+...
                         f3_tau_h*cAdapt.vb(k);
 
        end
    otherwise
        error('Unknown method passed');
end

%Finds steady state output voltage to horizontal cells then incorporates
%2.8 ms delay as in van Hateren 2005 by simply shifting all values over and
%adding copies of the initial value to the front, and cutting off extra
%information at the end.
cAdapt.steadyH=cAdapt.vh(length(cAdapt.timebase));
delay=round(2.82/cAdapt.dt);
cAdapt.vh=[cAdapt.vh(1)*ones(1,delay) cAdapt.vh(1:length(cAdapt.vh)-delay)];
