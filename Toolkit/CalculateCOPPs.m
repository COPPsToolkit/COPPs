function [x,shocks_u,shocks_a,exitflag, Lopt] = CalculateCOPPs(m,VarData,~,~,p,ref_x)
%--------------------------------------------------------------------------
%
% NOTES:
%   The model:    H_B x_{t-1} + H_C x_{t} + H_F E_{t}[x_{t+1}] = PSI u_{t}
%   The solution: x_{t} = B X_{t-1} + PHI u_{t}
%
%   Let x = [x_{1} ; ... ; x_{H}] be a (m.nx*H x 1) vector of states
%   and u = [u_{1} ; ... ; u_{H}] be a (m.ns*H x 1) vector of unant shocks
%   and a = [a_{1} ; ... ; a_{H}] be a (m.ns*H x 1) vector of antic shocks
%
%   Then P x = C x_{0} + Q a + R u where, for H = 4:
%
%    P = [  I  0  0  0 ]   C = [ B ]   Q = [ PHI  F*PHI  F^2*PHI  F^3*PHI ]
%        [ -B  I  0  0 ]       [ 0 ]       [  0    PHI    F*PHI   F^2*PHI ]
%        [  0 -B  I  0 ]       [ 0 ]       [  0     0      PHI     F*PHI  ]
%        [  0  0 -B  I ]       [ 0 ]       [  0     0       0       PHI   ]
%
%    R = [ PHI   0    0    0  ]
%        [  0   PHI   0    0  ]
%        [  0    0   PHI   0  ]
%        [  0    0    0   PHI ]
%
%   For a given path of model variables (x), this can be inverted to give
%   the shocks that deliver it, etc.
%
% Authors: O. de Groot, F. Mazelis, R. Motto, A. Ristiniemi, 2019-2021
% Acknowledgements: J. Barrdear
%--------------------------------------------------------------------------

exitflag = 1;

% Inverter = m.aInverter;
P = m.aP;
Q = m.aQ;
%R = m.aR;
%C = m.aC;

% Dimensions
T  = size(VarData,2);
H  = p.H_policy;
D  = p.D_term;
nl = length(m.LossVarIndexes);  % Number of target variables
ns = m.ns;                      % Number of shocks in the model
ni = length(m.InstrShkIndexes); % Number of instrument variables
     
shocks_u = zeros(ns,T);
shocks_a = zeros(ns,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINED OPTIMAL POLICY PROJECTIONS
% Authors: O. de Groot, F. Mazelis, R. Motto, A. Ristiniemi, 2019-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the Effective Lower Bound
if isfield(p.UserProvided,'RateMin')
    if isnumeric(p.PolicyRateStSt)
    rSS  = p.PolicyRateStSt;
    elseif sum(strcmp(p.PolicyRateStSt,m.ParamNames))==1
    rSS  = m.Params(strcmp(p.PolicyRateStSt,m.ParamNames),1);
    else
        error('The steady state of the nominal policy rate params.PolicyRateStSt needs to be provided as a numeric value or alternatively as a character array indicating the corresponding parameter name in the existing mod file. The parameter name you indicated could not be found in the mod file.')
    end
    RateMin = (p.RateMin-rSS)/4;
else
    RateMin = p.RateMin;
end

%% Construct IRFs for instrument[s] (AX1, AX2) and targets (AL)
AX1 = NaN(ni*T,T);
if ni == 2
    AX2 = NaN(ni*T,T);
end
AL = NaN(T,T*nl);

for k = 1:ni
    for i = 1:T
        shocks_a = zeros(m.ns,T);
        shocks_a(m.InstrShkIndexes(k),i) = 1;
        x = P\Q*shocks_a(:);
        x = reshape(x,m.nx,T);
        AX1((k-1)*T+i,:) = x(m.InstrVarIndexes(1),:);       % e.g. Response of instrument 1 (and 2) to instrument 1 shock 
        if ni == 2
            AX2((k-1)*T+i,:) = x(m.InstrVarIndexes(2),:);   % e.g. Response of instrument 1 (and 2) to instrument 2 shock
        end
        xtemp = x(m.LossVarIndexes,:)';
        AL((k-1)*T+i,:)= xtemp(:);
    end
end

if ni==1
    NewShkIndex = 1:H;
elseif ni==2
    NewShkIndex = [1:H,T+1:T+H];
end

AX1     = -AX1';
AX1     = AX1(:,NewShkIndex);
if ni == 2
    AX2 = -AX2';
    AX2 = AX2(:,NewShkIndex);
end
AL = AL(NewShkIndex,:);
LB = ref_x(m.LossVarIndexes,1:T)';
LB = LB(:);

%% Construct loss function weighting matrix
beta_vec = m.DiscountFactor.^(0:T-1);
Omega    = kron(diag(p.LossFunctionWeights),diag(beta_vec));

%% Construct constraints
if strcmp(p.Constraint,'ON')
    if isempty(RateMin)
        rmin = min(ref_x(m.InstrVarIndexes(1),1:T))*ones(1,T);
    else
        rmin = RateMin*ones(1,T);
    end
    RefRtemp = ref_x(m.InstrVarIndexes(1),1:T)-rmin;
    RefRmin  = min(RefRtemp,p.RateDevMin);
    RefRmax  = p.RateDevMax;
    bX       = [RefRmin';RefRmax'];
    AX       = [AX1;-AX1];
end

if ni == 2 && strcmp(p.Constraint,'ON')
    if isempty(p.QeMin)
        qmin = min(ref_x(m.InstrVarIndexes(2),1:T))*ones(1,T);
    else
        qmin = p.QeMin*ones(1,T);
    end
    if isempty(p.QeMax)
        qmax = max(ref_x(m.InstrVarIndexes(2),1:T))*ones(1,T);
    else
        qmax = p.QeMax*ones(1,T);
    end
    RefQ1   =  ref_x(m.InstrVarIndexes(2),1:T)-qmin;
    RefQ2   = -ref_x(m.InstrVarIndexes(2),1:T)+qmax;
    RefQmin = min(p.QeDevMin,RefQ1);
    RefQmax = min(p.QeDevMax,RefQ2);
    bX      = [bX;RefQmin';RefQmax'];
    AX      = [AX;AX2;-AX2];
end

options = optimoptions('quadprog','Display','off');
f.A     = [];
f.b     = [];
f.Aeq   = [];
f.beq   = [];
f.lb    = [];
f.ub    = [];
f.x0    = [];

%% Solve optimal policy problem
if strcmp(p.PolicyType,'COM')
    f.H = AL*Omega*AL';
    f.H = (f.H+f.H')/2; % Ensures Hessian is symmetric
    f.f = LB'*Omega*AL';
    if strcmp(p.Constraint,'ON')
        f.A = AX;
        f.b = bX;
        [shocks,~,exitflag,~,~] = quadprog(f.H,f.f,f.A,f.b,f.Aeq,f.beq,f.lb,f.ub,f.x0,options);
    else
        shocks = f.H\(-f.f');
    end
elseif strcmp(p.PolicyType,'DIS')
    npm    = H/D;
    step   = 0;
    Mdist  = 1;
    xguess = zeros(H*ni,1);
    xnew   = zeros(H*ni,1);
    ALtemp = reshape(AL',T,nl,H*ni);
    
    while Mdist>p.Crit && step<p.MaxIter
        LBguess = LB + (AL')*xguess;
        LBguess = reshape(LBguess,T,nl);
        
        step = step+1;
        for i = 1:npm
            j = (i-1)*D+1;
            LBguesstemp = LBguess(j:end,:);
            ALunant1    = ALtemp(1:end+1-j,:,1:D);
            [d1,d2,d3]  = size(ALunant1);
            ALunant     = reshape(ALunant1,d1*d2,d3,1);
            if ni == 2
                ALunant2 = squeeze(ALtemp(1:end+1-j,:,H+1:H+D));
                [d1,d2,d3]  = size(ALunant2);
                ALunant  = [ALunant,reshape(ALunant2,d1*d2,d3,1)];
            end
            Omegatemp = kron(diag(p.LossFunctionWeights),diag(beta_vec(1:end-j+1)));
            
            if strcmp(p.Constraint,'OFF')
                xnewtemp = -(ALunant'*Omegatemp*ALunant)\(LBguesstemp(:)'*Omegatemp*ALunant)';
            elseif strcmp(p.Constraint,'ON')
                baseR    = ref_x(m.InstrVarIndexes(1),1:T)'-AX1*xguess;
                RefRtemp = baseR'-rmin;
                RefRmin  = min(RefRtemp,(p.RateDevMin'-AX1*xguess)'); 
                RefRmax  = min((p.RateDevMax'+AX1*xguess)',(p.RateDevMax'+AX1*xguess)');              
                if ni == 2
                    baseQ = ref_x(m.InstrVarIndexes(2),1:T)'-AX2*xguess;
                    RefQtempmin =  baseQ'-RefQmin;
                    RefQtempmax = -baseQ'+RefQmax;     
                end             
                f.H = ALunant'*Omegatemp*ALunant;
                f.H = (f.H+f.H')/2; % Ensures Hessian is symmetric
                f.f = LBguesstemp(:)'*Omegatemp*ALunant;
                if ni == 1
                    f.A = [AX1(1:end+1-j,1:D);-AX1(1:end+1-j,1:D)];
                    f.b = [RefRmin(j:end)';RefRmax(j:end)'];
                elseif ni == 2
                    f.A = [AX1(1:end+1-j,[1:D,H+1:H+D]);
                          -AX1(1:end+1-j,[1:D,H+1:H+D]);
                           AX2(1:end+1-j,[1:D,H+1:H+D]);
                          -AX2(1:end+1-j,[1:D,H+1:H+D])];
                    f.b = [RefRmin(j:end)';RefRmax(j:end)';RefQtempmin(j:end)';RefQtempmax(j:end)'];
                end
                f.b(f.b<0)=0;
                [xnewtemp,~,exitflag,~,~] = quadprog(f.H,f.f,f.A,f.b,f.Aeq,f.beq,f.lb,f.ub,f.x0,options);
            end
            xnew(j:j+D-1) = xnewtemp(1:D);
            if ni == 2
                xnew(H+j:H+j+D-1) = xnewtemp(D+1:2*D);
            end
        end
        xupdated = xguess + xnew;
        
        Mdist = max(abs(xnew'));
        fprintf('step %i   %f\n',[step,Mdist])
        xguess = (1-p.Update)*xguess + p.Update*xupdated;
    end
    if step == p.MaxIter
        disp('Discretion did not solve')
        return
    else
        shocks = xguess;
    end
end

%% Calculate loss
Lbase = LB'*Omega*LB;
Lopt  = Lbase + shocks'*AL*Omega*AL'*shocks + 2*LB'*Omega*AL'*shocks;
Lopt  = Lopt/Lbase;
Lbase = 1;

%% Collect shocks and the optimal path
shocks_a = 0*shocks_a;
shocks_a(m.InstrShkIndexes(1),1:H) = shocks(1:H);
if ni==2
    shocks_a(m.InstrShkIndexes(2),1:H) = shocks(H+1:end);
end
x = P\Q*shocks_a(:);
x = reshape(x,m.nx,T);

end