// The basic New Keynesian model
// Gali ch. 3
// Oliver de Groot (July 2020)
//------------------------------------------
// Preamble
//------------------------------------------
// Variables
var     i pi og rn a;
var     r pinf dr r_ann pinf_4q rr rr_ann ui up pinf_ann rn_ann; 
varexo  er eps_a eps_ui eps_up;

// Parameters
parameters beta epsilon theta sigma rho phi alpha phi_pi phi_y eta PSI_yan THETA lambda kappa rho_v rho_a;
parameters sder rho_u;

beta    = 0.99;
sigma   = 1;
phi     = 1;
alpha   = 1/3;
epsilon = 6;
eta     = 4;
theta   = 2/3;
phi_pi  = 1.5;
phi_y   = 0.5/4;
PSI_yan = (1+phi)/(sigma*(1-alpha)+phi+alpha);
THETA   = (1-alpha)/(1-alpha+alpha*epsilon);
lambda  = (1-theta)*(1-beta*theta)*THETA/theta;
kappa   = lambda*(sigma+(phi+alpha)/(1-alpha));
rho     = 1;
rho_v   = 0.5;
rho_a   = 0.9;
sder    = 0.0625;
rho_u   = 0.8;

//------------------------------------------
// Model
//------------------------------------------
model(linear);
    i = 0*rho + phi_pi*pi + phi_y*og + sder*er; 
    
    og = og(+1) - 1/sigma*(i-pi(+1)-rn);
    pi = beta*pi(+1) + kappa*og + ui + up;

    rn = -a;

    // Shocks
    a  = rho_a*a(-1)  + eps_a;  // natural rate
    ui =                eps_ui; // transitory cost-push
    up = rho_u*up(-1) + eps_up; // persistent cost-push

    //OPP
    r       = i;
    pinf    = pi;
    dr      = r-r(-1);
    r_ann   = 4*r;
    pinf_4q = pinf+pinf(-1)+pinf(-2)+pinf(-3);
    rr      = r-pi(1);
    rr_ann  = 4*rr;
    pinf_ann= 4*pinf;
    rn_ann  = 4*rn;
end;

//------------------------------------------
// Steady State
//------------------------------------------
check;

//------------------------------------------
// Shocks
//------------------------------------------
shocks;
    var er    = 0*1;
    var eps_a = 0*1;
    var eps_ui= 1*1;
    var eps_up= 0*1;
end;

//------------------------------------------
// Computation
//------------------------------------------
stoch_simul(irf=200,nograph,noprint);

@#ifdef SMOOTHER
    varobs r pinf og dr;
    calib_smoother(datafile='baseline_costpush_iid.xlsx');
@#endif
