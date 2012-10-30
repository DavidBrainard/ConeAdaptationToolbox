function cAdapt=calcXHateren(cAdapt,method)
% cAdapt=calcXHateren(cAdapt,method)
%
% Second step of Hateren calculations containing calcium feedback loop
% Using equations found on pages 333-334 of van Hateren 2005.
%
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
%
% 7/xx/12  ekf  Wrote it.  Contact ekfrance@gmail.com
% 7/24/12  dhb, ekf        Going through it.  Not done yet.

% Starting values based on Fortran code preallocates and defines initial values
% values as in Fortran code
cAdapt.B=(cAdapt.cb+(cAdapt.kb*cAdapt.background))*ones(1,length(cAdapt.timebase));
cAdapt.Tx=1./cAdapt.B;

% Use numerical method to find steady-state value of C
options = optimset('fsolve');
options = optimset(options,'Diagnostics','off','Display','off');
C=fsolve(@(x) x-(1/(cAdapt.cb+(cAdapt.kb*cAdapt.background))/(1+(cAdapt.ac*x)^cAdapt.nc))^cAdapt.nx,14,options);
cAdapt.C=C*ones(1,length(cAdapt.timebase)); 
xcoeff=(1/(cAdapt.cb+(cAdapt.kb*cAdapt.background)))*(1/(1+(cAdapt.ac*C)^cAdapt.nc));
cAdapt.X=xcoeff*ones(1,length(cAdapt.timebase));
cAdapt.Ios=(xcoeff^cAdapt.nx)*ones(1,length(cAdapt.timebase)); 
cAdapt.alpha=(1/(1+(cAdapt.ac*C)^cAdapt.nc))*ones(1,length(cAdapt.timebase));

%% Calculations
% Iterative equations derived from equations on pg 333-334
%
% Formulas using arma coefficients are as given at the header of table 1 in
% van Hateren 2008

switch (method)
    case 'iterative'
        % Compute cGmp concentration and outer segment current iteratively
        % Equations derived from eqs 12, 13, 14, and 15 on pgs 333-334 van Hateren
        % 2005
        % Found by solving for change in variable and adding to variable's
        % previous value
        for i = 2:length(cAdapt.timebase)
            cAdapt.B(i)=cAdapt.cb+(cAdapt.kb*cAdapt.Estar(i));
            cAdapt.Tx(i)=1/cAdapt.B(i);
            cAdapt.X(i)=cAdapt.X(i-1)+cAdapt.dt*(cAdapt.alpha(i-1)-cAdapt.X(i-1)/cAdapt.Tx(i));
            cAdapt.Ios(i)=cAdapt.X(i)^cAdapt.nx;
            cAdapt.C(i)=cAdapt.C(i-1)+(cAdapt.dt/cAdapt.Tc)*(cAdapt.Ios(i)-cAdapt.C(i-1));
            cAdapt.alpha(i)=1/(1+(cAdapt.ac*cAdapt.C(i))^cAdapt.nc);
        end
        cAdapt.aoverb=cAdapt.alpha./cAdapt.B;

    case 'arma'
        % Calculate concentration and outer segment current using ARMA coefficients
        % Calculated according to formula for First-Order Hold in table 1
        % of van Hateren 2008
        %
        %%% Calculate non-time-dependent ARMA coeffs
                temp_Tc = cAdapt.Tc/cAdapt.dt;
                f1_tau_c = exp(-1/temp_Tc);             %-a1
                f2_tau_c = temp_Tc-(1+temp_Tc)*f1_tau_c;%b1
                f3_tau_c = 1-temp_Tc+temp_Tc*f1_tau_c;  %b0
        %%% Loop through, calculating ARMA coeffs for Tx each step 
        for i=2:length(cAdapt.timebase)
                cAdapt.B(i)=cAdapt.cb+(cAdapt.kb*cAdapt.Estar(i));
                cAdapt.Tx(i)=1/cAdapt.B(i);
                    temp_Tx = cAdapt.Tx(i)/cAdapt.dt;
                    f1_tau_x = exp(-1/temp_Tx);             %-a1
                    f2_tau_x = temp_Tx-(1+temp_Tx)*f1_tau_x;%b1
                    f3_tau_x = 1-temp_Tx+temp_Tx*f1_tau_x;  %b0
            cAdapt.X(i)=f1_tau_x*cAdapt.X(i-1)+...
                        f2_tau_x*cAdapt.alpha(i-1)/cAdapt.B(i-1)+...
                        f3_tau_x*cAdapt.alpha(i-1)/cAdapt.B(i);
            cAdapt.Ios(i)=cAdapt.X(i)^cAdapt.nx;
            cAdapt.C(i)=f1_tau_c*cAdapt.C(i-1)+...
                        f2_tau_c*cAdapt.Ios(i-1)+...
                        f3_tau_c*cAdapt.Ios(i);
            cAdapt.alpha(i)=1/(1+(cAdapt.ac*cAdapt.C(i))^cAdapt.nc);
        end
    otherwise
    error('Unknown method passed');
end
        