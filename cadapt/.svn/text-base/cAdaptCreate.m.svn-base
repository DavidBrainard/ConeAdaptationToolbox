function cAdapt=cAdaptCreate(whichModel,whichParams,varargin)
% cAdapt=cAdaptCreate(whichModel,whichParams,varargin)
%
% Create a defaut cAdapt structure.  There are
% some choices about which parameters to use,
% because these differ with in the van Hateren
% paper.
%
% whichModel: which underlying model to use.
%   'vanHateren' - The van Hateren (2005 etc.) model.
%
% whichParams: which set of parameters to set
%   'generic' - The generic parameters used for van Hateren's (2005) Figure 1.
%   'fig7' - The parameter used to create van Hateren's (2005) Figure 7.
%
% Optional arguments as string,value pairs.  These can
% be anything that cAdaptSet understands. 
%
% See also cAdaptSet

% Bookkeeping
cAdapt.name = 'default';
cAdapt.type = 'cAdapt';
cAdapt.model = whichModel;

%% Model parameters
switch (cAdapt.model)
    case 'vanHateren'
        switch(whichParams)
            case 'generic'
                %Values as specified in original fortran code. These should produce
                %the graphs in figure B of page 336 when used with the validation
                %script.
                cAdapt.Tr=3.4;  %ms
                cAdapt.Kr=1/cAdapt.Tr;  %1/ms
                cAdapt.Kr=1/cAdapt.Tr;  %1/ms
                cAdapt.Te=8.7;  %ms
                cAdapt.Ke=1/cAdapt.Te;  %1/ms
                cAdapt.Rstar=[];        %Activated Rhodopsin
                cAdapt.Estar=[];        %Activated PDE
                cAdapt.cb=2.80*10^(-3);  %1/ms
                cAdapt.kb=1.63*10^(-4);  %1/(ms*td)
                cAdapt.B=[];    %Activity of PDE
                cAdapt.X=[];    %Concentration of cGMP
                cAdapt.Tc=2.89;    %ms
                cAdapt.ac=9.08*10^(-2);    %au
                cAdapt.eta=1;   %Scaling constant
                cAdapt.nx=1;    %??
                cAdapt.nc=4;    %??
                cAdapt.Tm=4;    %ms
                cAdapt.gamma=.678;
                cAdapt.ais=7.09*10^(-2);   %au
                cAdapt.Tis=90;  %ms
                cAdapt.gt=125;  %au
                cAdapt.vk=-10;  %mv
                cAdapt.vn=3;    %mv
                cAdapt.vI=19.7;   %mv
                cAdapt.mu=.733;
                cAdapt.Ta=250;  %ms
                cAdapt.T1=4;    %ms
                cAdapt.T2=4;    %ms
                cAdapt.Th=20;   %ms
                %Spatial stuff
                cAdapt.Titd=10000; %ms
                cAdapt.Titp=10000;   %ms
                cAdapt.s=.4;
                cAdapt.ch=.25;  %au
                cAdapt.ih=-20;  %au
                cAdapt.tpmax=25;    %ms
                cAdapt.ip=-25;  %au
                cAdapt.cp=.25;  %a
            case 'fig7'
                %Values as stated in caption of figure 7 on page 337.
                %Theoretically, these should produce the graphs in figure 7, but so
                %far have not been able to.
                cAdapt.Tr=.49;  %ms
                cAdapt.Kr=1/cAdapt.Tr;  %1/ms
                cAdapt.Te=16.3;  %ms
                cAdapt.Ke=1/cAdapt.Te;  %1/ms
                cAdapt.Rstar=[];        %Activated Rhodopsin
                cAdapt.Estar=[];        %Activated PDE
                cAdapt.cb=2.80*10^(-3);  %1/ms
                cAdapt.kb=1.63*10^(-4);  %1/(ms*td)
                cAdapt.B=[];    %Activity of PDE
                cAdapt.X=[];    %Concentration of cGMP
                cAdapt.Tc=2.89;    %ms
                cAdapt.ac=9.08*10^(-2);    %au
                cAdapt.eta=1;   %Scaling constant
                cAdapt.nx=1;    %??
                cAdapt.nc=4;    %??
                cAdapt.Tm=4;    %ms
                cAdapt.gamma=.678;
                cAdapt.ais=7.09*10^(-2);   %au
                cAdapt.Tis=56.9;  %ms
                cAdapt.gt=151.1;  %au
                cAdapt.vk=-10;  %mv
                cAdapt.vn=3;    %mv
                cAdapt.vI=19.7;   %mv
                cAdapt.mu=.733;
                cAdapt.Ta=250;  %ms
                cAdapt.T1=4;    %ms
                cAdapt.T2=4;    %ms
                cAdapt.Th=20;   %ms
                %Spatial stuff
                cAdapt.Titd=10000; %ms
                cAdapt.Titp=10000;   %ms
                cAdapt.s=.4;
                cAdapt.ch=.25;  %au
                cAdapt.ih=-20;  %au
                cAdapt.tpmax=25;    %ms
                cAdapt.ip=-25;  %au
                cAdapt.cp=.25;  %au
                
            otherwise
                error('Unknown paramter set');
        end
    otherwise
        error('Unknown model specified');
end

% Handle any additional arguments via cAdaptSet
if ~isempty(varargin)
    if rem(length(varargin),2)~=0
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        cAdapt = cAdaptSet(cAdapt,varargin{ii},varargin{ii+1});
    end
end

