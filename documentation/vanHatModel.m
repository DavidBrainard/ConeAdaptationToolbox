function ValTable=vanHatModel(background,contrast,pulseduration,whichParams,timestruct)
% function ValTable=vanHatModel(background,contrast,pulseduration,whichParams,timestruct)
% whichParams has options 'generic' 'fig7'
% timestruct requires fields:
%   timestruct.dt -  timestep in ms
%   timestruct.timestart - first passed time
%   timestruct.timeon - time stimulus start
%   timestruct.timeend - last passed time% 
%  hc_widefield
%
%       Fortran source code translated to Matlab. Developed with an Intel Fortran90 compiler on Linux
%
%       note: implicit type convention: variables starting with i..n are
%       integers, all others are reals
%
%       9 Jan 2005, author J.H. van Hateren
%       Belongs to "A cellular and molecular model of response kinetics
%       and adaptation in primate cones and horizontal cells", J.Vision (2005)
%
%       Remarks:
%        -function rtbis not yet included (see below); insert that,
%         and change the statements for resp_c and resp_vs at
%         'set pre-stimulus values'; the values there now are only
%         correct for stim_prev=100 td
%       - for higher gains in the feedback loops, nrate may need to
%         be larger
%       - check the ARMA coefficients for very high or very small values
%         of tau, switch to double precision variables if necessary
%
% Edited for Brainard lab by DHB, EKF
% contact ekfrance@gmail.com


nrate=1/timestruct.dt*1000;                %#timesteps/sec (delta_t=100 us)
nstim_max=5*nrate;          %maximum stimulus length 5 sec

%         real stim(2)                           %I in model
%         real resp_r(2)                         %signal after tau_r
%         real resp_e(2)                         %E* in model
%         real beta(2)                           %beta in model
%         real resp_q(2)                         %Q in model
%         real tau_x(2)                          %tau_x in model
%         real resp_w(2)                         %alpha/beta in model
%         real gain_x(2)                         %alpha in model
%         real resp_x(2)                         %X in model
%         real resp_os(2)                        %I_os in model
%         real resp_c(2)                         %C in model
%         real resp_vc(2)                        %signal after division by g_i
%         real resp_is(2)                        %V_is in model
%         real gain_is(2)                        %g_is in model
%         real gain_i(2)                         %g_i in model
%         real resp_vs(2)                        %V_s in model
%         real resp_ic(2)                        %I_t in model
%         real resp_tr(2)                        %signal after tau_1
%         real resp_ih(2)                        %signal after tau_2 (V_b)
%         real resp_h(2)                         %V_h in model
%         real resp_h0(2)                        %V_is_prime in model
ratefrac=nrate/1000;                             %conversion 1 ms to 100 us timebase

%
%       make stimulus
%
nlen=fix((timestruct.timeend-timestruct.timestart+timestruct.dt)*ratefrac);  % stimulus duration
stim_ar = zeros(1,nlen);                         %stimulus array
resp_ar = zeros(1,nlen);                         %response array
resp_kl = zeros(1,nlen);                         %auxiliary variable used for delay


if timestruct.timestart < 0                 %make stimulus start after t=0
    correction=(0-timestruct.timestart)*ratefrac;
else
    correction=0;
end
ip1=correction+(timestruct.timeon+timestruct.dt)*ratefrac;   %ms stimulus latency
ip2=ip1+fix(pulseduration*ratefrac)-1;           
stim_dc=background;                                     
stim_prev=background;                                   %assumed illuminance before start stimulus
stim_ar(1:nlen)=stim_dc;
stim_ar(ip1:ip2)=stim_dc+contrast*stim_dc;             

%
%       parameters, in ms where applicable, converted to delta_t by ratefrac;
%       values are those of Figs.7 and 6A of the article

switch(whichParams)
    case 'generic'
        %Values as specified in original fortran code. These should produce
        %the graphs in figure B of page 336 when used with the validation
        %script.
        tau_r=3.4*ratefrac;                 %parameter tau_r in model
        tau_e=8.7*ratefrac;                 %parameter tau_e in model
        beta_0=2.80e-3;                      %parameter c_beta in model
        g_ex0=1.63e-4;                       %parameter k_beta in model
        rnx=1;                               %parameter n_x in model
        rnc=4;                               %parameter n_c in model
        tau_c=2.89*ratefrac;                 %parameter tau_c in model
        a_c=9.08e-2;                         %parameter a_c in model
        tau_vc=4.*ratefrac;                  %parameter tau_m in model
        gamma_is=0.678;                      %parameter gamma_is in model
        tau_is=90*ratefrac;                %parameter tau_is in model
        a_is=7.09e-2;                        %parameter a_is in model
        ripmax=125;                        %parameter g_t in model
        vk=-10;                              %parameter v_k in model
        vn=3;                                %parameter v_n in model
        tau_tr=4.*ratefrac;                  %parameter tau_1 in model
        tau_ih=4.*ratefrac;                  %parameter tau_2 in model
        tau_h=20.*ratefrac;                  %parameter tau_h in model
        rdel=2.82*ratefrac;                  %overall delay
        vh0=19.7;                            %parameter V_I in model
        rho=0.733;                           %parameter mu in model
        tau_h0=250.*ratefrac;                %parameter tau_a in model
    case 'fig7'
        %Values as stated in caption of figure 7 on page 337.
        tau_r=.49*ratefrac;                 %parameter tau_r in model
        tau_e=16.8*ratefrac;                 %parameter tau_e in model
        beta_0=2.80e-3;                      %parameter c_beta in model
        g_ex0=1.63e-4;                       %parameter k_beta in model
        rnx=1;                               %parameter n_x in model
        rnc=4;                               %parameter n_c in model
        tau_c=2.89*ratefrac;                 %parameter tau_c in model
        a_c=9.08e-2;                         %parameter a_c in model
        tau_vc=4.*ratefrac;                  %parameter tau_m in model
        gamma_is=0.678;                      %parameter gamma_is in model
        tau_is=56.9*ratefrac;                %parameter tau_is in model
        a_is=7.09e-2;                        %parameter a_is in model
        ripmax=151.1;                        %parameter g_t in model
        vk=-10;                              %parameter v_k in model
        vn=3;                                %parameter v_n in model
        tau_tr=4.*ratefrac;                  %parameter tau_1 in model
        tau_ih=4.*ratefrac;                  %parameter tau_2 in model
        tau_h=20.*ratefrac;                  %parameter tau_h in model
        rdel=2.82*ratefrac;                  %overall delay
        vh0=19.7;                            %parameter V_I in model
        rho=0.733;                           %parameter mu in model
        tau_h0=250.*ratefrac;                %parameter tau_a in model
    otherwise
        error('Invalid parameter set. Please use ''generic'' or ''fig7''')
end
%
%       calculate ARMA coefficients of fixed filters
%
f1_tau_r=exp(-1./tau_r);
f2_tau_r=(tau_r-(1.+tau_r)*f1_tau_r);
f3_tau_r=(1.-tau_r+tau_r*f1_tau_r);

f1_tau_e=exp(-1./tau_e);
f2_tau_e=(tau_e-(1.+tau_e)*f1_tau_e);
f3_tau_e=(1.-tau_e+tau_e*f1_tau_e);

f1_tau_c=exp(-1./tau_c);
f2_tau_c=(tau_c-(1.+tau_c)*f1_tau_c);
f3_tau_c=(1.-tau_c+tau_c*f1_tau_c);

f1_tau_vc=exp(-1./tau_vc);
f2_tau_vc=(tau_vc-(1.+tau_vc)*f1_tau_vc);
f3_tau_vc=(1.-tau_vc+tau_vc*f1_tau_vc);

f1_tau_is=exp(-1./tau_is);
f2_tau_is=a_is*(tau_is-(1.+tau_is)*f1_tau_is);
f3_tau_is=a_is*(1.-tau_is+tau_is*f1_tau_is);

f1_tau_tr=exp(-1./tau_tr);
f2_tau_tr=tau_tr-(1.+tau_tr)*f1_tau_tr;
f3_tau_tr=1.-tau_tr+tau_tr*f1_tau_tr;

f1_tau_h0=exp(-1./tau_h0);
f2_tau_h0=tau_h0-(1.+tau_h0)*f1_tau_h0;
f3_tau_h0=1.-tau_h0+tau_h0*f1_tau_h0;

%main loop
for it=1:nlen
    
    %       ncurr and nprev determine which element of, e.g., resp_r, is
    %       the current one, resp_r(1) or resp_r(2)
    %
    if (it == 1)
        ncurr=1;                          %current
        nprev=2;                          %previous
    else
        nkl=ncurr;                        %swap values ncurr and nprev
        ncurr=nprev;
        nprev=nkl;
    end
    
    if (it == 1)                          %set pre-stimulus values
        options = optimset('fsolve');
        options = optimset(options,'Diagnostics','off','Display','off');
        resp_c(nprev)=fsolve(@(x)  x-(1./(1.+(a_c*x)^rnc)/(beta_0+g_ex0*stim_prev))^rnx, 12, options);
        %resp_c(nprev)=14.12579;           %delete if rtbis is available
        stim(nprev)=stim_prev;
        resp_r(nprev)=stim(nprev);
        resp_e(nprev)=resp_r(nprev);
        beta(nprev)=beta_0+g_ex0*resp_e(nprev);
        resp_q(nprev)=1/beta(nprev);
        tau_x(nprev)=resp_q(nprev)*ratefrac;
        gain_x(nprev)=1./(1.+(a_c*resp_c(nprev))^rnc);
        resp_w(nprev)=gain_x(nprev)*resp_q(nprev);
        resp_x(nprev)=gain_x(nprev)*resp_q(nprev);
        resp_os(nprev)=resp_x(nprev)^rnx;
        resp_vc(nprev)=(resp_os(nprev)/a_is)^(1/(1+gamma_is));
        resp_is(nprev)=resp_vc(nprev);
        gain_is(nprev)=resp_is(nprev)^gamma_is;
        gain_i(nprev)=a_is*gain_is(nprev);
        resp_h0(nprev)=resp_is(nprev);
        gtau=(resp_h0(nprev)/vh0)^rho;                   %a_I in model
        gripmax=ripmax/gtau;
        st_is=resp_is(nprev);                %transferred to steady_vs via common
        resp_vs(nprev)=fsolve(@(x) x-(st_is-gripmax/(1.+exp(-(x-vk)/vn))),-13,options);
        %resp_vs(nprev)=-12.92790;            %delete if rtbis is available
        resp_ic(nprev)=gripmax/(1.+exp(-(resp_vs(nprev)-vk)/vn));
        resp_tr(nprev)=resp_ic(nprev);
        resp_ih(nprev)=resp_tr(nprev);
        resp_h(nprev)=resp_ih(nprev);
    end
    
    stim(ncurr)=stim_ar(it);
    
    resp_r(ncurr)=f1_tau_r*resp_r(nprev)+ ...
        f2_tau_r*stim(nprev)+ ...
        f3_tau_r*stim(ncurr);
  
    resp_e(ncurr)=f1_tau_e*resp_e(nprev)+ ...
        f2_tau_e*resp_r(nprev)+ ...
        f3_tau_e*resp_r(ncurr);
    
    beta(ncurr)=beta_0+g_ex0*resp_e(ncurr);
    resp_q(ncurr)=1./beta(ncurr);
    
    tau_x(ncurr)=resp_q(ncurr)*ratefrac;
    f1_tau_x=exp(-1./tau_x(ncurr));                   %ARMA coefficients tau_x
    f2_tau_x=tau_x(ncurr)-(1.+tau_x(ncurr))*f1_tau_x;
    f3_tau_x=1.-tau_x(ncurr)+tau_x(ncurr)*f1_tau_x;
    
    resp_x(ncurr)=f1_tau_x*resp_x(nprev)+ ...
        gain_x(nprev)*f2_tau_x*resp_q(nprev)+ ...
        gain_x(nprev)*f3_tau_x*resp_q(ncurr);
    
    resp_w(ncurr)=gain_x(nprev)*resp_q(ncurr);        %not necessary, only for figure
    
    resp_os(ncurr)=resp_x(ncurr)^rnx;
    
    resp_c(ncurr)=f1_tau_c*resp_c(nprev)+ ...
        f2_tau_c*resp_os(nprev)+ ...
        f3_tau_c*resp_os(ncurr);
    
    gain_x(ncurr)=1./(1.+(a_c*resp_c(ncurr))^rnc);
    
    resp_vc(ncurr)=resp_os(ncurr)/gain_i(nprev);
    
    resp_is(ncurr)=f1_tau_vc*resp_is(nprev)+ ...
        f2_tau_vc*resp_vc(nprev)+ ...
        f3_tau_vc*resp_vc(ncurr);
    
    gain_is(ncurr)=resp_is(ncurr)^gamma_is;
    
    gain_i(ncurr)=f1_tau_is*gain_i(nprev)+ ...
        f2_tau_is*gain_is(nprev)+f3_tau_is*gain_is(ncurr);
    
    resp_h0(ncurr)=f1_tau_h0*resp_h0(nprev)+ ...
        f2_tau_h0*resp_is(nprev)+ ...
        f3_tau_h0*resp_is(ncurr);
    
    gtau=(resp_h0(ncurr)/vh0)^rho;
    
    gripmax=ripmax/gtau;
    gtau_ih=tau_ih*gtau;                          %tau_2_prime in model
    gtau_h=tau_h*gtau;                            %tau_h_prime in model
    
    resp_vs(ncurr)=resp_is(ncurr)-resp_h(nprev);
    
    resp_ic(ncurr)=gripmax/(1.+exp(-(resp_vs(ncurr)-vk)/vn));
    
    resp_tr(ncurr)=f1_tau_tr*resp_tr(nprev)+ ...
        f2_tau_tr*resp_ic(nprev)+ ...
        f3_tau_tr*resp_ic(ncurr);
    
    f1_tau_ih=exp(-1./gtau_ih);                %ARMA coefficients tau_2_prime
    f2_tau_ih=gtau_ih-(1.+gtau_ih)*f1_tau_ih;
    f3_tau_ih=1.-gtau_ih+gtau_ih*f1_tau_ih;
    
    f1_tau_h=exp(-1./gtau_h);
    f2_tau_h=gtau_h-(1.+gtau_h)*f1_tau_h;
    f3_tau_h=1.-gtau_h+gtau_h*f1_tau_h;
    
    resp_ih(ncurr)=f1_tau_ih*resp_ih(nprev)+ ...
        f2_tau_ih*resp_tr(nprev)+ ...
        f3_tau_ih*resp_tr(ncurr);
    
    resp_h(ncurr)=f1_tau_h*resp_h(nprev)+ ...
        f2_tau_h*resp_ih(nprev)+ ...
        f3_tau_h*resp_ih(ncurr);
    
    resp_ar(it)=resp_h(ncurr);                %output of resp_h (replace for obtaining other signals)

    ValTable{1}(it)=stim(ncurr);
    ValTable{2}(it)=resp_r(ncurr);
    ValTable{3}(it)=resp_e(ncurr);
    ValTable{4}(it)=beta(ncurr);
    ValTable{5}(it)=resp_q(ncurr);
    ValTable{6}(it)=resp_x(ncurr);
    ValTable{7}(it)=resp_c(ncurr);
    ValTable{8}(it)=resp_w(ncurr);
    ValTable{9}(it)=resp_is(ncurr);
    ValTable{10}(it)=gain_i(ncurr);
    ValTable{11}(it)=resp_vs(ncurr);
    ValTable{12}(it)=resp_ic(ncurr);
    ValTable{13}(it)=resp_ih(ncurr);
    ValTable{14}(it)=resp_os(ncurr);
    
    
    
    
end


      %delay response by interpolation

for it=1:nlen
    %for each timestep
    rel=it-rdel;
        %rel=current timestep-delay in ms
    if (rel >= 0)
        %as long as current timestep is bigger than delay
        iel1=fix(rel);
        %round the new current timestep
        iel2=iel1+1;
        %new variable one greater than new timestep
        d1=rel-iel1;
        %finds the amount that got rounded off
        d2=1-d1;
        %1-amount that was rounded off
    else
        iel1=fix(rel);
        iel2=iel1-1;
        d1=iel1-rel;
        d2=1.-d1;
    end
    if (iel1 < 1) 
        iel1=1; 
    end
    %if new timestep is too small, make it smallest possible val
    if (iel1 > nlen) 
        iel1=nlen; 
    end
    %if new timestep is too big, make it largest possible val
    if (iel2 < 1) 
        iel2=1; 
    end
    %
    if (iel2 > nlen) 
        iel2=nlen; 
    end
    resp_kl(it)=d2*resp_ar(iel1)+d1*resp_ar(iel2);
    %new value at that time is weighted average of previous values at new
    %timestep and next timestep
end
resp_ar(1:nlen)=resp_kl(1:nlen);
ValTable{14}=resp_ar;

%%%%%%%%%%%%FOR EMMA'S SANITY %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%NOT IN ORIGINAL CODE%%%%%%%%%%%%%%%
 clear resp* gain*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


