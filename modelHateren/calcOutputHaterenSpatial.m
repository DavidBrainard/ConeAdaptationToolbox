function cAdapt=calcOutputHaterenSpatial(cAdapt,method)
% cAdapt=calcoutputHateren(cAdapt,method)
% Final step of Hateren calculations containing cone-horizontal cell
% feedback loop
% Calculates output to horizontal and bipolar cells, accounting for feedback from
% horizontal cells
% Using equations found on pages 335-336 of van Hateren 2005
%
% 7/2/12: Output needs horizontal shift (overall delay as given in van
% Hateren 2005) EKF

%% Preallocations and initial values
% Values as given in fortran code
cAdapt.visprime=cAdapt.vis;
cAdapt.aI=(cAdapt.visprime./cAdapt.vI).^cAdapt.mu;
gripmax=cAdapt.gt./cAdapt.aI;
%%%% Constant value for vs (12.92790) ONLY VALID IF background illuminance
%%%% set to 100td OR model is starting in negative time, so it hits steady
%%%% state by time 0.
cAdapt.vs=-12.92790*ones(1,length(cAdapt.timebase));
cAdapt.It=(gripmax./(1+exp(-(-12.92790-cAdapt.vk)/cAdapt.vn)));
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
        
        temp_Titd = cAdapt.Titd/cAdapt.dt;
        f1_tau_itd=exp(-1/temp_Titd);   %-a1
        f2_tau_itd=temp_Titd-(1+temp_Titd)*f1_tau_itd; %b1
        f3_tau_itd=1-temp_Titd+temp_Titd*f1_tau_itd;  %b0
        
        temp_Titp = cAdapt.Titp/cAdapt.dt;
        f1_tau_itp=exp(-1/temp_Titp);   %-a1
        f2_tau_itp=temp_Titp-(1+temp_Titp)*f1_tau_itp; %b1
        f3_tau_itp=1-temp_Titp+temp_Titp*f1_tau_itp;  %b0
        
        %%% Calculate non-time dependent spatial filter coefficients
        lambdaSTemp = cAdapt.lambdaS/cadapt.dt;
        g1_lambda_s = exp(-1/lambdaSTemp);
        g2_lambda_s = 1-g1_lambda_s;
        
        lambdaLTemp = cAdapt.lambdaL/cadapt.dt;
        g1_lambda_L = exp(-1/lambdaLTemp);
        g2_lambda_L = 1-g1_lambda_L;
        
        for k=2:length(cAdapt.timebase)
            %calculate aI in order to calculate two arma coeffs          
            cAdapt.visprime(k)=f1_tau_a*(cAdapt.visprime(k-1))+...
                               f2_tau_a*(cAdapt.vis(k-1))+...
                               f3_tau_a*(cAdapt.vis(k));
            cAdapt.aI(k)=(cAdapt.visprime(k)/cAdapt.vI)^cAdapt.mu;
            %calculate static arma coeffs         
            temp_T2 = cAdapt.aI(k)*cAdapt.T2/cAdapt.dt;
            f1_tau_2=exp(-1/temp_T2); %-a1clc
            f2_tau_2=temp_T2-(1+temp_T2)*f1_tau_2; %b1
            f3_tau_2=1-temp_T2+temp_T2*f1_tau_2;  %b0

            temp_Th = cAdapt.aI(k)*cAdapt.Th/cAdapt.dt;
            f1_tau_h=exp(-1/temp_Th); %-a1
            f2_tau_h=temp_Th-(1+temp_Th)*f1_tau_h; %b1
            f3_tau_h=1-temp_Th+temp_Th*f1_tau_h;  %b0

            %calculate transmitter current, horizontal and bipolar cell
            %output
            cAdapt.viz(k)=cAdapt.vis(k)-cAdapt.visdark;
                                        %(dark_resp_os/a_is)**(1./(1.+gamma_is)) 
            
            temp_Tp = cAdapt.Tp(k-1)/cAdapt.dt;
            f1_tau_p=exp(-1/temp_Tp); %-a1
            f2_tau_p=temp_Tp-(1+temp_Tp)*f1_tau_p; %b1
            f3_tau_p=1-temp_Tp+temp_Tp*f1_tau_p;  %b0
            
            cAdapt.vp(k)=f1_tau_p*cAdapt.vp(k-1)+...
                         f2_tau_p*cAdapt.viz(k-1)+...
                         f3_tau_p*cAdapt.viz(k);
            cAdapt.vs(k)=cAdapt.vp(k)-cAdapt.vd(k-1);
            cAdapt.It(k)=(cAdapt.gt/cAdapt.aI(k))/(1+exp(-(cAdapt.vs(k)-cAdapt.vk)/cAdapt.vn));
            
            cAdapt.It1(k)=f1_tau_1*cAdapt.It1(k-1)+...
                          f2_tau_1*cAdapt.It(k-1)+...
                          f3_tau_1*cAdapt.It(k);
            cAdapt.ItpostTitd(k)=f1_tau_itd*cAdapt.It1(k-1)+...
                          f2_tau_itd*cAdapt.It(k-1)+...
                          f3_tau_itd*cAdapt.It(k);
            cAdapt.ItpostTitp(k)=f1_tau_itp*cAdapt.It1(k-1)+...
                          f2_tau_itp*cAdapt.It(k-1)+...
                          f3_tau_itp*cAdapt.It(k);  
            cAdapt.gh(k)=1/(1+exp(cAdapt.ch*(cAdapt.It(k)-cAdapt.ih)));
            
            cAdapt.Tp(k)=cAdapt.tpmax/(1+exp(cAdapt.cp*(cAdapt.It(k)-cAdapt.ip)));
            cAdapt.vb(k)=f1_tau_2*cAdapt.vb(k-1)+...
                         f2_tau_2*cAdapt.It1(k-1)+...
                         f3_tau_2*cAdapt.It1(k);
            cAdapt.vh1(k)=f1_tau_h*cAdapt.vh(k-1)+...
                         f2_tau_h*cAdapt.vb(k-1)+...
                         f3_tau_h*cAdapt.vb(k);
            cAdapt.vhs(k)=g1_lambda_s*cAdapt.vhs(k-1)+...
                          g2_lambda_s*cAdapt.vh1(k);
            cAdapt.vhL(k)=g1_lambda_L*cAdapt.vhL(k-1)+...
                          g2_lambda_L*cAdapt.vh1(k);
            cAdapt.vh(k)=%Some combination of vhs and vhL
            %resp_h(ic,nc3(ic))=wlambda_s*resp_hs(ic,nc3(ic))+
                %(1.-wlambda_s)*resp_hl(ic,nc3(ic))
            cAdapt.vd(k)=cAdapt.gh*cAdapt.vh+vdbias;
            %resp_d(ic,nc3(ic))=gh(ic,nc3(ic))*resp_h(ic,np3(ic))+vdbias
 
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
