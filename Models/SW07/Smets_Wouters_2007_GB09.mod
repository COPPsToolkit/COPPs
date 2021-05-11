/*
 * This file provides replication files for 
 * Smets, Frank and Wouters, Rafael (2007): "Shocks and Frictions in US Business Cycles: A Bayesian
 * DSGE Approach", American Economic Review, 97(3), 586-606, that are compatible with Dynare 4.2.5 onwards
 *
 * To replicate the full results, you have to get back to the original replication files available at
 * https://www.aeaweb.org/articles.php?doi=10.1257/aer.97.3.586 and include the respective estimation commands and mode-files.
 *
 * Notes: Please see the header to the Smets_Wouters_2007_45.mod for more details and a fully documented version.
 *
 * This file was originally written by Frank Smets and Rafeal Wouters and has been updated by
 * Johannes Pfeifer. 
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright (C) 2007-2013 Frank Smets and Raf Wouters
 * Copyright (C) 2013-15 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can receive a copy of the GNU General Public License
 * at <http://www.gnu.org/licenses/>.
 */



var labobs robs pinfobs dy dc dinve dw ewma epinfma zcapf rkf kf pkf cf 
    invef yf labf wf rrf mc zcap rk k pk c inve y lab pinf w r a b g qs ms  
    spinf sw kpf kp
    ygap dr r_aux
    pinf_4q
    obs_pinf_4q
    obs_r
    obs_r_ann
    dy_4q
    ;    
 
varexo ea eb eg eqs em epinf ew;  
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
    czcap csadjcost ctou csigma chabb  cfc % ccs cinvs
    cindw cprobw cindp cprobp csigl clandaw 
     crpi crdy cry crr % crdpi
    crhoa  crhob crhog crhols crhoqs crhoms crhopinf crhow  % crhoas
    ctrend cg
    beta_ss ss_r_ann;

// fixed parameters
ctou=.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;

% Mode from SW07 AER
constepinf = 0.81;
constepinf = 2/4; % Fed price stability mandate (2009)
constelab = -0.1;
ctrend = 0.43; 

calfa=.19;
% cbeta=.9995;
csigma=1.39;
cfc=1.61;
cgy=0.51;

csadjcost= 5.48;
chabb=    0.71;    
cprobw=   0.73;
csigl=    1.92;
cprobp=   0.65;
% cprobw=   0.9;
% cprobp=   0.9;
cindw=    0.59;
cindp=    0.22;
czcap=    0.54;
crpi=     2.03;
crr=      0.81;
cry=      0.08;
crdy=     0.22;

crhoa=    0.95;
crhob=    0.18;
crhog=    0.97;
crhols=   0.71;
crhoqs=   0.7165;
% crhoas=1; 
crhoms=0.12;
crhopinf=0.9;
crhow=0.97;
cmap = 0.74;
cmaw  = 0.88;


constebeta = 0.16; 
beta_ss = 1/(1+constebeta/100);

ss_r_ann = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100*4;

@#ifdef FlatPC
    cprobw=   0.73+1*.15;
    cprobp=   0.65+1*.15;
@#endif

var qeobs; 

model(linear); 
//deal with parameter dependencies; taken from usmodel_stst.mod 
#cpie=1+constepinf/100;
#cgamma=1+ctrend/100 ;
#cbeta=1/(1+constebeta/100);

#clandap=cfc;
#cbetabar=cbeta*cgamma^(-csigma);
#cr=cpie/(cbeta*cgamma^(-csigma));
#crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
#cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
//cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*((cbeta^(-1))*(cgamma^csigma) - (1-ctou))^calfa))^(1/(1-calfa));
#cikbar=(1-(1-ctou)/cgamma);
#cik=(1-(1-ctou)/cgamma)*cgamma;
#clk=((1-calfa)/calfa)*(crk/cw);
#cky=cfc*(clk)^(calfa-1);
#ciy=cik*cky;
#ccy=1-cg-cik*cky;
#crkky=crk*cky;
#cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
#cwly=1-crk*cky;

#conster=(cr-1)*100;

        r =  crpi*(1-crr)*pinf
               +cry*(1-crr)*(y-yf)     
               +crdy*(y-yf-y(-1)+yf(-1))
               +crr*r(-1)
               +ms  ;

// flexible economy

        0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
        zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
        rkf =  (wf)+labf-kf ;
        kf =  kpf(-1)+zcapf ;
        invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
        cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
        yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
        yf = cfc*( calfa*kf+(1-calfa)*labf +a );
        wf = csigl*labf   +(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
        kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy

        mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
        zcap =  (1/(czcap/(1-czcap)))* rk ;
        rk =  w+lab-k ;
        k =  kp(-1)+zcap ;
        inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
        c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
        y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
        y = cfc*( calfa*k+(1-calfa)*lab +a );
        pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
        w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
               + 1*sw ;
        a = crhoa*a(-1)  + ea;
        b = crhob*b(-1) + eb;
        g = crhog*(g(-1)) + eg + cgy*ea;
        qs = crhoqs*qs(-1) + eqs;
        ms = crhoms*ms(-1) + em;
        spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
            epinfma=epinf;
        sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
            ewma=ew; 
        kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

// measurment equations

dy=y-y(-1)+ctrend;
dc=c-c(-1)+ctrend;
dinve=inve-inve(-1)+ctrend;
dw=w-w(-1)+ctrend;
pinfobs = 1*(pinf) + constepinf;
robs =    1*(r) + conster;
labobs = lab + constelab;

ygap = y-yf;
r_aux = r;
dr = r - r_aux(-1);

pinf_4q = pinf + pinf(-1) + pinf(-2) + pinf(-3);
obs_pinf_4q = pinf_4q + 4*constepinf;
obs_r = r + conster;
obs_r_ann = 4*obs_r;
dy_4q = dy + dy(-1) + dy(-2) + dy(-3); 

qeobs = 0;

end; 

steady_state_model;
dy=ctrend;
dc=ctrend;
dinve=ctrend;
dw=ctrend;
pinfobs = constepinf;
robs = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
labobs = constelab;

obs_pinf_4q = 4*constepinf;
obs_r = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
obs_r_ann = 4*obs_r;
dy_4q = 4*ctrend;

ewma     =0;
epinfma  =0;
zcapf    =0;
rkf      =0;
kf       =0;
pkf      =0;
cf       =0;
invef    =0;
yf       =0;
labf     =0;
wf       =0;
rrf      =0;
mc       =0;
zcap     =0;
rk       =0;
k        =0;
pk       =0;
c        =0;
inve     =0;
y        =0;
lab      =0;
pinf     =0;
w        =0;
r        =0;
a        =0;
b        =0;
g        =0;
qs       =0;
ms       =0;
spinf    =0;
sw       =0;
kpf      =0;
kp       =0;
ygap     =0;
dr       =0;
r_aux    =0;
pinf_4q  =0;
end;

shocks;
var ea;
stderr 0.45;
var eb;
stderr 0.24;
var eg;
stderr 0.52;
var eqs;
stderr 0.45;
var em;
stderr 0.24;
var epinf;
stderr 0.14;
var ew;
stderr 0.24;
end;

@#ifdef SMOOTHER
    varobs 
    % dy 
    dy_4q
    ygap
    dc dinve labobs 
    % pinfobs 
    obs_pinf_4q 
    % dw 
    robs
    ;

    @#ifndef SteepFwd
    calib_smoother(datafile='usmodel_data.xls',diffuse_filter);
    @#endif

    @#ifdef SteepFwd
    calib_smoother(datafile='usmodel_data_SteepFwd.xls',diffuse_filter);
    @#endif

@#endif

@#ifdef IRFs

    stoch_simul(order=1
        ,irf=40
        ,irf_shocks =(em) 
        ,nograph
        ,noprint
        ) 
        ygap pinf_4q obs_r_ann obs_r  
    ;

@#endif