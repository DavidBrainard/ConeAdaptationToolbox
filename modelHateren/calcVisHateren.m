function cAdapt=calcVisHateren(cAdapt,method)
% cAdapt=calcVisHateren(cAdapt,method)
% Third step of Hateren calculations containing inner segment processes
% Using equations found on pages 334 of van Hateren 2005 
% We have a few different numerical methods implemented
%   'iterative' - Straightforward implementation of difference equations.
%   'arma' - Smarter implementation of difference equations.  Uses a 
%            little calculus to improve numerical accuracy.  
%
% The arma method is explained in van Hateren (2008).
%
% Notes: 
%   Van Hateren had a variable called background, which we think might represent
%   the stimulus before time 0.  Or else a spatial background for some
%   future implementation.  In any case, it was always set to the first stimulus value,
%   and we've set up the code on the assumption that this will always be true.

% Contact ekfrance@gmail.com


%% Preallocations and initial values
% Values as given in fortran code
cAdapt.vis=(cAdapt.Ios./cAdapt.ais).^(1/(1+cAdapt.gamma));
cAdapt.IoverGi=cAdapt.vis;
cAdapt.gis=cAdapt.vis.^cAdapt.gamma;
cAdapt.gi=cAdapt.gis.*cAdapt.ais;


%% Calculations
% Iterative equations derived from equations on pg 334
%
% Formulas using arma coefficients are as given at the header of table 1 in
% van Hateren 2008
switch (method)
    case 'iterative'
        % Calculates inner segment voltage
        % Equations found by solving for change in variable and adding to variable's
        % previous value
        for k=2:length(cAdapt.timebase)
            cAdapt.gis(k)=cAdapt.ais*(cAdapt.vis(k-1)^cAdapt.gamma);
            cAdapt.gi(k)=cAdapt.gi(k-1)+(cAdapt.dt/cAdapt.Tis)*(cAdapt.gis(k)-cAdapt.gi(k-1));
            cAdapt.vis(k)=cAdapt.vis(k-1)+(cAdapt.dt/cAdapt.Tm)*(cAdapt.Ios(k)/cAdapt.gi(k)-cAdapt.vis(k-1));
        end
    case 'arma'
        % Calculates inner segment voltage
        % Coefficients calculated according to formula for First-Order Hold
        % as given in table 1 of van Hateren 2008
        %
        %%% Calculate non-time dependent coefficients
        temp_Tm = cAdapt.Tm/cAdapt.dt;
        f1_tau_m = exp(-1/temp_Tm);             %-a1
        f2_tau_m = temp_Tm-(1+temp_Tm)*f1_tau_m;%b1
        f3_tau_m = 1-temp_Tm+temp_Tm*f1_tau_m;  %b0
        
        temp_Tis = cAdapt.Tis/cAdapt.dt;
        f1_tau_is = exp(-1/temp_Tis);             %-a1
        f2_tau_is = cAdapt.ais*(temp_Tis-(1+temp_Tis)*f1_tau_is);%b1
        f3_tau_is = cAdapt.ais*(1-temp_Tis+temp_Tis*f1_tau_is);  %b0
        
        %%% Loop through, calculating time-dependent arma coefficient each
        %%% step
        for i=2:length(cAdapt.timebase)
            cAdapt.IoverGi(i)=(cAdapt.Ios(i)/cAdapt.gi(i-1));
            cAdapt.vis(i)=f1_tau_m*(cAdapt.vis(i-1))+...
                          f2_tau_m*cAdapt.IoverGi(i-1)+...
                          f3_tau_m*cAdapt.IoverGi(i);
            cAdapt.gis(i)=(cAdapt.vis(i)^cAdapt.gamma);
            cAdapt.gi(i)=f1_tau_is*cAdapt.gi(i-1)+...
                      f2_tau_is*cAdapt.gis(i-1)+...
                      f3_tau_is*cAdapt.gis(i);
        end         
    otherwise
           error('Unknown method passed');
end

