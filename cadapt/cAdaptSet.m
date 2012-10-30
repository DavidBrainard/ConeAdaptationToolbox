function cAdapt = cAdaptSet(cAdapt,parm,val,varargin)
% cAdapt = cAdapt(cAdapt,parm,val,varargin)
%
% Set cAdapt structure parameters.
%
% See also: cAdaptCreate
%
% Parameters:
%
%  Bookkeeping
%   'name' - Name of this object
%   'type' - Type of this object, should always be 'cAdapt'
%   'model' - Which model is implemented in this particular cAdapt structure
%
%  Stimulus
%   'time'- Timebase for input (ms)
%   'stimulus'- Stimulus value over time (units tbd)
%   'dt'- Timestep for calculations
%
%  Model parameters
%

% Arg checks and parse.
%
% The switch on what we're setting is broken out into several pieces
% below to allow use of cells, and so that autoindent does something
% reasonable with our block comment style.
if ~exist('parm','var') || isempty(parm), error('Parameter must be defined'); end
if ~exist('val','var'), error('val must be defined'); end
%parm = ieParamFormat(parm);

%% Initialize flags
cAdapt.OUTPUT_STALE = true;
DIDASET = false;

%% Bookkeeping
switch parm
    % Bookkeeping
    case 'name'
        % This specific object's name
        cAdapt.name = val;
        DIDASET = true;
        
    case 'type'
        % Type should always be 'cAdapt'
        if (~strcmp(val,'cAdapt'))
            error('Can only set type of wvf structure to ''wvf''');
        end
        cAdapt.type = val;
        DIDASET = true;
        
        case 'model'
        % Type should always be 'vanHateren'
        if (~strcmp(val,'vanHateren'))
            error('Can only set type of model structure to ''vanHateren''');
        end
        cAdapt.model = val;
        DIDASET = true;
end

%% Stimulus
%
% These specify the stimulus to compute for
switch parm
    case 'timebase'
        % Timebase in ms of the stimulus sequence
        cAdapt.timebase = val;
        cAdapt.OUTPUT_STALE = true;
        DIDASET = true;
        
    case 'background'
        % Background illuminance in trolands
        cAdapt.background = val;
        cAdapt.OUTPUT_STALE = true;
        DIDASET = true;
        
    case 'stimulus'
        % Stimulus.  Must be same length as timebase.  Units to be defined.
        cAdapt.stimulus = val;
        cAdapt.OUTPUT_STALE = true;
        DIDASET = true;
    case 'dt'
        % Timestep for calculations, in ms.
        cAdapt.dt = val;
        cAdapt.OUTPUT_STALE = true;
        DIDASET = true;
    
end

%% Catch the case where we don't know about the requested parameter
switch (parm)
    otherwise
        if (~DIDASET)
            error('Unknown parameter %s\n',parm);
        end
end

return
