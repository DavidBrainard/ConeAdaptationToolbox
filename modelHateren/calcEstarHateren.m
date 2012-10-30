function cAdapt = calcEstarHateren(cAdapt,method)
% cAdapt = calcEstarHateren(cAdapt,method)
% method has possible values
%   'convolution' - uses convultion to solve analytically
%   'iterative' - straightforward stepping through differential eqs
%   'arma' - uses auto-regressive moving average (a numerical method) based
%   on coefficients calculated according to formula for First-Order Hold
%   as given in table 1 of van Hateren 2008
%
% Calculates Estar values based on parameter values specified in the
% passed adaptation object.
%
% Equations as listed on pages 332-333 of van Hateren, 2005
% Calculates in turn amount of activated Rhodopsin and active PDE
%
% Notes: 
%   Van Hateren had a variable called background, which we think might represent
%   the stimulus before time 0.  Or else a spatial background for some
%   future implementation.  In any case, it was always set to the first stimulus value,
%   and we've set up the code on the assumption that this will always be true.
%
% Contact ekfrance@gmail.com
%
% 6/xx/12  ekf  Wrote it.
% 6/24/12  dhb  Reworked a bit, and worried about how gain is affected by dt.

% Set default method
if (nargin < 2 || isempty(method))
    method = 'convolution';
end

%% Create impluse responses for Rstar and Estar stages.
%
% Each of these is an exponentially decaying function of
% time, with different time constants for each.
%
% Since these are simple exponentials, we can figure out length
% of the needed timebase from the time constant.  We just
% go out long enough so that the impuluse response gets
% nice and small.
%
% The premultiply by cAdapt.dt keeps the overall gain correct
% given the timebase on which we are going to do our convolutions.
expImpulseTimeFactor = 10;
cAdapt.RstarImpulseResponseTimebase = 0:cAdapt.dt:(expImpulseTimeFactor/cAdapt.Kr);
cAdapt.RstarImpulseResponse = cAdapt.dt*cAdapt.Kr*exp(-cAdapt.RstarImpulseResponseTimebase*cAdapt.Kr);

cAdapt.EstarImpulseResponseTimebase = 0:cAdapt.dt:(expImpulseTimeFactor/cAdapt.Ke);
cAdapt.EstarImpulseResponse = cAdapt.dt*cAdapt.Ke*exp(-cAdapt.EstarImpulseResponseTimebase*cAdapt.Ke);

%% Compute Rstar and Estar, according to passed method
switch (method)
    case 'convolution'  
        % Rstar and Estar calculation via convolution.
        %
        % We chop off the times at the end of the convolutions
        % to keep the computed value on the same timebase as
        % the input.
        Rstar = conv(cAdapt.stimulus,cAdapt.RstarImpulseResponse);
        cAdapt.Rstar = Rstar(1:length(cAdapt.timebase));
        
        Estar = conv(cAdapt.Rstar,cAdapt.EstarImpulseResponse);
        cAdapt.Estar = Estar(1:length(cAdapt.timebase));
        
    case 'iterative'
        % Just simulate out the difference equations.
        %
        % This is a formula a derived by discretizing
        % Eq 7 of van Hateren (2005) (and similarly for
        % Estar.  
        cAdapt.Rstar(1) = cAdapt.background;
        cAdapt.Estar(1) = cAdapt.background;
        for i = 2:length(cAdapt.timebase)
            cAdapt.Rstar(i) = cAdapt.Rstar(i-1)+(cAdapt.dt/cAdapt.Tr)*(cAdapt.stimulus(i)-cAdapt.Rstar(i-1));
            cAdapt.Estar(i) = cAdapt.Estar(i-1)+(cAdapt.dt/cAdapt.Te)*(cAdapt.Rstar(i)-cAdapt.Estar(i-1));
        end
        
    case 'arma'
        % This is our transliteration of van Hateren's fortran code.
        % Coefficients calculated according to formula for First-Order Hold
        % as given in table 1 of van Hateren 2008

        temp_Tr = cAdapt.Tr/cAdapt.dt;
        f1_tau_r = exp(-1/temp_Tr);
        f2_tau_r = temp_Tr-(1+temp_Tr)*f1_tau_r;
        f3_tau_r = 1-temp_Tr+temp_Tr*f1_tau_r;
        
        temp_Te = cAdapt.Te/cAdapt.dt;
        f1_tau_e = exp(-1/temp_Te);     %a1
        f2_tau_e = temp_Te-(1+temp_Te)*f1_tau_e;    %b1
        f3_tau_e = 1-temp_Te+temp_Te*f1_tau_e;      %b0
        
        cAdapt.Rstar(1) = cAdapt.background;
        cAdapt.Estar(1) = cAdapt.background;
        for i = 2:length(cAdapt.timebase)
            cAdapt.Rstar(i) = f1_tau_r*cAdapt.Rstar(i-1) + f2_tau_r*cAdapt.stimulus(i-1) + f3_tau_r*cAdapt.stimulus(i);
            cAdapt.Estar(i) = f1_tau_e*cAdapt.Estar(i-1) + f2_tau_e*cAdapt.Rstar(i-1) + f3_tau_e*cAdapt.Rstar(i);
        end
                
    otherwise
        error('Unknown method passed');
end

end