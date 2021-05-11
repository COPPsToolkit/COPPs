%% Sims and WU - Baseline
%
%  Evaluating Central Banks’ Tool Kit: Past, Present, and Future
%  JME 2020
%
%  Implemented for the COPPs toolkit
%  **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
%  "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
% -------------------------------------------------------------------------

var Y C I G L Ld Ih Yw K u w wr mrs Pi Pir pk pw f1 f2 x1 x2 vp vw Q QB RF
RB Rd Rre Rtr M1 M2 mu Lam fw f b re d n lambda phi Omega fcb bcb bG A theta
logY logC logI logLd logPi logRd exr logQ
rtr dr dy pinfobs robs robs_ann qeobs
dy_4q dy_L1 dy_L2 dy_L3 pinf_4q pinf_L1 pinf_L2 pinf_L3
;

varexo eA er et eB eG ef eb et2;

parameters beta kappa spF spB psi epsip epsiw alpha delta0 g bb bcs bcbGs 
by levs sigma Rds RFs RBs Qs QBs M1s M2s pws Ks Ys Is ws delta1 Cs mus mrss 
chi fws res bcbs fcbs fs bGs bs ns ds Delta phis thetas X Omegas lams phip 
phiw psik delta2 eta rhor phipi phiy sr rhof rhob sf sb rhoA sA rhot st 
rhoG sG rhoB sB
dyobsss pinfobsss robss qeobsss
;

Rds = beta^(-1);
RFs = (1+spF)^(1/4)*Rds;
RBs = (1+spB)^(1/4)*Rds;
Qs = (RFs - kappa)^(-1);
QBs = (RBs - kappa)^(-1);
M1s = beta/(Qs*(1-beta*kappa));
M2s = 1 + (M1s-1)*psi;
pws = (epsip-1)/epsip;
Ks = (alpha*pws/(M2s*(1/beta - (1-delta0))))^(1/(1-alpha));
Ys = Ks^(alpha);
Is = delta0*Ks;
ws = (1-alpha)*pws*Ks^(alpha);
delta1 = (alpha*pws*Ks^(alpha-1)/M2s);
Cs = (1-g)*Ys - Is;
mus = (1/Cs)*(1-beta*bb)/(1-bb);
mrss = ((epsiw-1)/epsiw)*ws;
chi = mus*mrss;
fws = psi*Is/(Qs*(1-kappa));
res = bcs*Ys;
bcbs = bcbGs*res/QBs;
fcbs = (res - QBs*bcbs)/Qs;
fs = fws - fcbs;
bGs = by*Ys/QBs;
bs = bGs - bcbs;
ns = (Qs*fs + QBs*bs + res)/levs;
ds = Qs*fs + QBs*bs + res - ns;
Delta = (RBs - Rds)/(RFs - Rds);
phis = (Qs*fs + Delta*QBs*bs)/ns;
thetas = (1 - sigma + phis*beta*(1-sigma)*(RFs - Rds))/((1-sigma)*phis - beta*sigma*phis^2*(RFs - Rds));
X = ns - sigma*( (RFs-Rds)*Qs*fs + (RBs - Rds)*QBs*bs + Rds*ns);
Omegas = 1 - sigma + sigma*phis*thetas;
lams = (thetas/(beta*(RFs - Rds)*(1-sigma + sigma*phis*thetas)) - 1)^(-1);



@#ifdef LOAD_ESTIMATED_PARAMS
load param_sw;
set_param_value('alpha',alpha);
set_param_value('beta',beta);
set_param_value('delta0',delta0);
set_param_value('delta1',delta1);
set_param_value('delta2',delta2);
set_param_value('eta',eta);
set_param_value('kappa',kappa);
set_param_value('spF',spF);
set_param_value('spB',spB);
set_param_value('psi',psi);
set_param_value('epsip',epsip);
set_param_value('epsiw',epsiw);
set_param_value('g',g);
set_param_value('bb',bb);
set_param_value('bcs',bcs);
set_param_value('bcbGs',bcbGs);
set_param_value('by',by);
set_param_value('levs',levs);
set_param_value('sigma',sigma);
set_param_value('Rds',Rds);
set_param_value('RFs',RFs);
set_param_value('RBs',RBs);
set_param_value('Qs',Qs);
set_param_value('QBs',QBs);
set_param_value('M1s',M1s);
set_param_value('M2s',M2s);
set_param_value('pws',pws);
set_param_value('Ks',Ks);
set_param_value('Ys',Ys);
set_param_value('Is',Is);
set_param_value('ws',ws);
set_param_value('Cs',Cs);
set_param_value('mus',mus);
set_param_value('mrss',mrss);
set_param_value('chi',chi);
set_param_value('fws',fws);
set_param_value('res',res);
set_param_value('bcbs',bcbs);
set_param_value('fcbs',fcbs);
set_param_value('fs',fs);
set_param_value('bGs',bGs);
set_param_value('bs',bs);
set_param_value('ns',ns);
set_param_value('ds',ds);
set_param_value('Delta',Delta);
set_param_value('phis',phis);
set_param_value('thetas',thetas);
set_param_value('X',X);
set_param_value('Omegas',Omegas);
set_param_value('lams',lams);
set_param_value('phip',phip);
set_param_value('phiw',phiw);
set_param_value('psik',psik);
set_param_value('eta',eta);
set_param_value('rhor',rhor);
set_param_value('phipi',phipi);
set_param_value('phiy',phiy);
set_param_value('sr',sr);
set_param_value('rhof',rhof);
set_param_value('rhob',rhob);
set_param_value('sf',sf);
set_param_value('sb',sb);
set_param_value('rhoA',rhoA);
set_param_value('sA',sA);
set_param_value('rhot',rhot);
set_param_value('st',st);
set_param_value('rhoG',rhoG);
set_param_value('sG',sG);
set_param_value('rhoB',rhoB);
set_param_value('sB',sB);
set_param_value('dyobsss',0.5023);   % mean - 2.0091
set_param_value('robss',1.3100);     % mean - 5.2399
set_param_value('qeobsss',0);        % ss is zero
@#endif


% COPPs Additions
var pinf r ygap obs_r_ann obs_pinf_4q;
parameters ss_r_ann beta_ss;

demeaned_policy_stst = 100*log(Rds);
ss_r_ann = (robss+demeaned_policy_stst)*4 ;
beta_ss = beta;

pinfobsss = 2 / 4;

var 
rtr_aux 
dbcb
bcb_aux
dqeobs
qeobs_aux
;

bcbs = 0; 

parameters
gamp gamw
;

% SW07 posterior means
gamp = 0.24 * 0;
gamw = 0.58 * 0;
phip = 0.75 + 0.15;
phiw = 0.75 + 0.15;


model;

%%%%%%
% Household
%%%%%%

% (1) SDF
Lam = beta*mu/mu(-1);

% (2) mu
mu = (C - bb*C(-1))^(-1) - beta*bb*(C(+1) - bb*C)^(-1);

% (3) Labor supply
chi*L^(eta) = mrs*mu;

% (4) Bonds
1 = Rd*Lam(+1)*Pi(+1)^(-1);

% (5) Wage-setting
wr = (epsiw/(epsiw-1))*f1/f2;

% (6) f1
% f1 = mrs*w^(epsiw)*Ld + phiw*Lam(+1)*Pi(+1)^(epsiw)*f1(+1);
f1 = mrs*w^(epsiw)*Ld + phiw*Lam(+1)*(Pi(+1)/Pi^gamw)^(epsiw)*f1(+1);

% (7) f2
% f2 = w^(epsiw)*Ld + phiw*Lam(+1)*Pi(+1)^(epsiw-1)*f2(+1);
f2 = w^(epsiw)*Ld + phiw*Lam(+1)*(Pi(+1)/Pi^gamw)^(epsiw-1)*f2(+1);

%%%%%
% Investment firm
%%%%%

% (8) Ihat
Ih = (1 - (psik/2)*(I/I(-1) - 1)^2)*I;

% (9) FOC I
1 = pk*(1 - (psik/2)*(I/I(-1) - 1)^2 - psik*(I/I(-1) - 1)*(I/I(-1))) + Lam(+1)*pk(+1)*psik*(I(+1)/I - 1)*(I(+1)/I)^2;


%%%%%%
% Retail firm
%%%%%%

% (10) Reset inflation
Pir = (epsip/(epsip - 1))*x1/x2;

% (11) x1
% x1 = pw*Y + phip*Lam(+1)*Pi(+1)^(epsip)*x1(+1);
x1 = pw*Y + phip*Lam(+1)*(Pi(+1)/Pi^gamp)^(epsip)*x1(+1);

% (12) x2
% x2 = Y + phip*Lam(+1)*Pi(+1)^(epsip-1)*x2(+1);
x2 = Y + phip*Lam(+1)*(Pi(+1)/Pi^gamp)^(epsip-1)*x2(+1);

%%%%%%%
% Wholesale firm
%%%%%%%%%

% (13) Labor demand
w = (1-alpha)*pw*A*(u*K(-1))^(alpha)*Ld^(-alpha);

% (14) Utilization
pk*M2*(delta1 + delta2*(u - 1)) = alpha*pw*A*(u*K(-1))^(alpha-1)*Ld^(1-alpha);

% (15) Capital
pk*M2 = Lam(+1)*(alpha*pw(+1)*A(+1)*(u(+1)*K)^(alpha-1)*u(+1)*Ld(+1)^(1-alpha) + (1-delta0-delta1*(u(+1) - 1)-(delta2/2)*(u(+1) - 1)^2)*pk(+1)*M2(+1));

% (16) Bonds
Q*M1 = Lam(+1)*Pi(+1)^(-1)*(1+kappa*Q(+1)*M1(+1));

% (17) M1 and M2
(M1-1)/(M2-1) = 1/psi;

% (18) Production
Yw = A*(u*K(-1))^(alpha)*Ld^(1-alpha);

% (19) Capital acumulation
K = Ih + (1-delta0-delta1*(u-1)-(delta2/2)*(u-1)^2)*K(-1);

% (20) Loan in advance
psi*pk*Ih = Q*(fw - kappa*Pi^(-1)*fw(-1));

%%%%%%%%%%%%%%%
% Financial Intermediaries
%%%%%%%%%%%%%%%

% (21) FOC private
Lam(+1)*(RF(+1) - Rd)*Pi(+1)^(-1)*Omega(+1) = theta*lambda/(1+lambda);

% (22) FOC government
Lam(+1)*(RB(+1) - Rd)*Pi(+1)^(-1)*Omega(+1) = Delta*theta*lambda/(1+lambda);

% (23) FOC reserves
Lam(+1)*(Rre - Rd)*Pi(+1)^(-1)*Omega(+1) = 0;

% (24) Omega
Omega = 1- sigma + sigma*phi*theta;

% (25) phi
phi = (Lam(+1)*Pi(+1)^(-1)*Omega(+1)*Rd)/(theta - Lam(+1)*Omega(+1)*Pi(+1)^(-1)*(RF(+1) - Rd));

% (26) Balance sheet
Q*f + QB*b + re = d + n;

% (27) Modified leverage constraint
phi = (Q*f + Delta*QB*b)/n;

% (28) Net worth evolution
n = sigma*Pi^(-1)*((RF - Rd(-1))*Q(-1)*f(-1) + (RB - Rd(-1))*QB(-1)*b(-1) + (Rre(-1) - Rd(-1))*re(-1) + Rd(-1)*n(-1)) + X;

%%%%%%%%%%%
% Aggregate conditions
%%%%%%%%%%

% (34) Price evolution
% 1 = (1-phip)*Pir^(1-epsip) + phip*Pi^(epsip-1);
1 = (1-phip)*Pir^(1-epsip) + phip*Pi(-1)^(gamp*(1-epsip))*Pi^(epsip-1);

% (35) Wage evolution
% w^(1-epsiw) = (1-phiw)*wr^(1-epsiw) + phiw*Pi^(epsiw-1)*w(-1)^(1-epsiw);
w^(1-epsiw) = (1-phiw)*wr^(1-epsiw) + phiw*Pi(-1)^(gamw*(1-epsiw))*Pi^(epsiw-1)*w(-1)^(1-epsiw);

% (36) Aggregate output
Yw = Y*vp;

% (37) Price dispersion
% vp = (1-phip)*Pir^(-epsip) + phip*Pi^(epsip)*vp(-1);
vp = (1-phip)*Pir^(-epsip) + phip*(Pi/(Pi(-1)^gamp))^(epsip)*vp(-1);

% (38) Labor supply / demand
L = Ld*vw;

% (39) Wage dispersion
% vw = (1-phiw)*(wr/w)^(-epsiw) + phiw*(w/w(-1))^(epsiw)*Pi^(epsiw)*vw(-1);
vw = (1-phiw)*(wr/w)^(-epsiw) + phiw*(w/w(-1))^(epsiw)*(Pi/(Pi(-1)^gamw))^(epsiw)*vw(-1);

% (40) Market-clearing private bonds
fw = f + fcb;

% (41) Market-clearing government bonds
bG = b + bcb;

% (42) Aggregate resource constraint
Y = C + I + G;

% (43) Return private bonds
RF = (1+kappa*Q)/Q(-1);

% (44) Return government bonds
RB = (1+kappa*QB)/QB(-1);

% (45) Productivity shock
log(A) = rhoA*log(A(-1)) + sA*eA;

% (46) Credit shock
log(theta) = (1-rhot)*log(thetas) + rhot*log(theta(-1)) + st*et + st*et2;

% (47) Government spending shock
log(G) = (1-rhoG)*log(g*Ys) + rhoG*log(G(-1)) + sG*eG;

% (48) Government debt shock
log(bG) = (1-rhoB)*log(bGs) + rhoB*log(bG(-1)) + sB*eB;

% OTHER STUFF

% (49) excess return
exr = 400*(log(RF(+1)) - log(Rd));

% (50) log output
logY = log(Y);

% (51) log consumption
logC = log(C);

% (52) log investment
logI = log(I);

% (53) Log Labor
logLd = log(Ld);

% (54) Net inflation
logPi = 400*log(Pi);

% (55) Net deposit rate
logRd = 400*log(Rd);

% (56) Log Q
logQ = log(Q);

% Observation equations
dy    = 100*(log(Y)-log(Y(-1))) + dyobsss; 
pinfobs = 100*log(Pi) + pinfobsss;       
robs    = 100*rtr     + robss;           
robs_ann = robs*4;
qeobs   = 100*bcb     + qeobsss;         

dy_L1   = dy(-1);
dy_L2   = dy_L1(-1);
dy_L3   = dy_L2(-1);
dy_4q   = dy + dy_L1 + dy_L2 + dy_L3;

pinf_L1 = pinfobs(-1);
pinf_L2 = pinf_L1(-1);
pinf_L3 = pinf_L2(-1);
pinf_4q = pinfobs + pinf_L1 + pinf_L2 + pinf_L3;

dr = 100*(rtr-rtr_aux);
rtr_aux = rtr(-1);
%%%%%%
%% Central Bank
%%%%%%%%%

rtr = log(Rtr);

% (29) Taylor rule
rtr = (1-rhor)*log(Rds) + rhor*rtr(-1) + (1-rhor)*(phipi*log(Pi) + phiy*(log(Y) - log(Y(-1)))) + sr*er;

% (32) Government holdings: Consider only government bond purchases in the optimal policy choice
bcb = (1-rhob)*bcbs + rhob*bcb(-1) + sb*eb;

% (30)
Rre = Rtr;

% (31) Private holdings
fcb = (1-rhof)*fcbs + rhof*fcb(-1) + sf*ef;

% (33) Balance sheet
Q*fcb + QB*bcb = re;

pinf = log(Pi);      
r = 100*rtr;
ygap = (Y - Ys) / Ys * 100; 
obs_r_ann = robs_ann;
obs_pinf_4q = pinf_4q;
dbcb = (bcb-bcb_aux);
bcb_aux = bcb(-1);
dqeobs = (qeobs-qeobs_aux);
qeobs_aux = qeobs(-1);
end;

initval;
Y = Ys;
Yw = Ys;
L = 1;
Ld = 1;
vp = 1;
vw = 1;
G = g*Ys;
I = Is;
C = Cs;
Ih = Is;
K = Ks;
u = 1;
w = ws;
wr = ws;
mrs = mrss;
Pi = 1;
Pir = 1;
pk = 1;
pw = pws;
f1 = (mrss*ws^(epsiw))/(1-phiw*beta);
f2 = (ws^(epsiw)/(1-phiw*beta));
x1 = (pws*Ys)/(1-phip*beta);
x2 = Ys/(1-phip*beta);
Q = Qs;
QB = QBs;
RF = RFs;
RB = RBs;
Rd = Rds;
Rre = Rds;
Rtr = Rds;
M1 = M1s;
M2 = M2s;
mu = mus;
Lam = beta;
fw = fws;
f = fs;
b = bs;
re = res;
d = ds;
n = ns;
lambda = lams;
phi = phis;
Omega = Omegas;
fcb = fcbs;
bcb = bcbs;
bG = bGs;
A = 1;
theta = thetas;
exr = 400*(log(RFs) - log(Rds));
logY = log(Ys);
logI = log(Is);
logC = log(Cs);
logLd = 0;
logPi = 0;
logRd = 400*log(Rds);
logQ = log(Qs);
end;

% steady;

shocks;
var eA = 1;
var eG = 1;
var er = 1;
var et = 1;
var eB = 1;
var eb = 1;
var ef = 1;
var et2 = 1;
end;

@#ifdef SMOOTHER
    varobs 
    % dy 
    ygap
    obs_pinf_4q 
    robs 
    qeobs
    ;
    calib_smoother(datafile='usmodel_data_fromSW07.xls');
@#endif

