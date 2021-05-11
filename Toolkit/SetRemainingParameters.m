function [params] = SetRemainingParameters(params)

%% Check required parameters

if ~isfield(params,'ModelDirectory')
    error('params.ModelDirectory needs to be provided')
end
if ~isfield(params,'ModelFilename')
    error('params.ModelFilename needs to be provided')
end
if ~isfield(params,'PastPeriods')
    error('params.PastPeriods needs to be provided')
end
if ~isfield(params,'PolicyInstrumentsAndShocks')
    error('params.PolicyInstrumentsAndShocks needs to be provided')
end
if ~isfield(params,'LossFunctionVariables')
    error('params.LossFunctionVariables needs to be provided')
end
if ~isfield(params,'PlannerDiscountFactor')
    error('params.PlannerDiscountFactor needs to be provided')
end

%% Optional parameters

% save user provided parameters
params.UserProvided = rmfield(params,{'ModelDirectory','ModelFilename','PastPeriods','PolicyInstrumentsAndShocks','LossFunctionVariables','PlannerDiscountFactor'});

%% Optimal policy
if ~isfield(params,'PolicyType')
    params.PolicyType = 'COM';     % Choice: COM, DIS      COM: Commitment, DIS: Discretion
end
if ~isfield(params,'LossFunctionWeights')
    params.LossFunctionWeights = ones(1,size(params.LossFunctionVariables,1));
end
if ~isfield(params,'T_loss')
    params.T_loss   = 40;            % Integer value         Number of periods for calculating welfare
end
if ~isfield(params,'H_policy')
    params.H_policy = params.T_loss;        % Integer value         Number of periods of anticipated shocks (must be less than p.T_loss)
end
if ~isfield(params,'T_full')
    params.T_full   = 100; % How far into the future to roll the forecast forward.
end
if ~isfield(params,'D_term')
    params.D_term   = 1; % Number of period in a policymaker's term of office
end
if strcmp(params.PolicyType,'DIS') && mod(params.H_policy,params.D_term) ~= 0
    error('params.D_term needs to be a factor of params.H_policy')
end

%% Inattention parameters (for details, see De Groot, Mazelis - Mitigating the Forward Guidance Puzxzle: Inattention, Credibility, Finite Planning Horizons and Learning)
if ~isfield(params,'TYPE')
    params.TYPE  = 'I';           % Choice: I II III IV   I: Inattention, II: Credibility, III: Finite planning horizon, IV: Learning
end
if ~isfield(params,'alpha')
    params.alpha = 1-0.3;         % TYPE I and II only:   1 = Full attention/credibility, 0 = Complete inattention/incredibility
end
if ~isfield(params,'N_planning')
    params.N_planning = 4;        % Type III only:        Integer values
end
if ~isfield(params,'beta1')
    params.beta1 = 1;             % Type IV  only:        See paper
end
if ~isfield(params,'beta2')
    params.beta2 = 5;             % Type IV  only:        See paper
end

%% Constraints
if ~isfield(params,'Constraint')
    params.Constraint = 'OFF';      % Choice: ON, OFF       ZLB constraint
end
if isfield(params,'RateMin') && ~isfield(params,'PolicyRateStSt')
        error('params.PolicyRateStSt needs to be provided')
elseif ~isfield(params,'RateMin')
    params.RateMin         = 0;        
    params.PolicyRateStSt  = 'NaN';        
end
if ~isfield(params,'RateDevMin')
    params.RateDevMin     = 100*ones(1,params.T_loss);
elseif isfield(params,'RateDevMin') & size(params.RateDevMin,2)==1
    params.RateDevMin     = params.RateDevMin*ones(1,params.T_loss);
elseif isfield(params,'RateDevMin') & size(params.RateDevMin,2)<params.T_loss
    params.RateDevMin(size(params.RateDevMin,2)+1:params.T_loss) = 100;
end
if ~isfield(params,'RateDevMax')
    params.RateDevMax     = 100*ones(1,params.T_loss);
elseif isfield(params,'RateDevMax') & size(params.RateDevMax,2)==1
    params.RateDevMax     = params.RateDevMax*ones(1,params.T_loss);
elseif isfield(params,'RateDevMax') & size(params.RateDevMax,2)<params.T_loss
    params.RateDevMax(size(params.RateDevMax,2)+1:params.T_loss) = 100;
end
if ~isfield(params,'QeMin')
    params.QeMin        = 0;         % Min QE: Leave at 0
end
if ~isfield(params,'QeMax')
    params.QeMax        = 100;       % Max QE: if [], the code will take the max. of the baseline path
end
if ~isfield(params,'QeDevMin')
    params.QeDevMin = 100*ones(1,params.T_loss);       % qe: Max. absolute deviation *below* the baseline path
elseif isfield(params,'QeDevMin') & size(params.QeDevMin,2)==1
    params.QeDevMin     = params.QeDevMin*ones(1,params.T_loss);
elseif isfield(params,'QeDevMin') & size(params.QeDevMin,2)<params.T_loss
    params.QeDevMin(size(params.QeDevMin,2)+1:params.T_loss) = 100;
end
if ~isfield(params,'QeDevMax')
    params.QeDevMax = 100*ones(1,params.T_loss);       % qe: Max. absolute deviation *above* the baseline path
elseif isfield(params,'QeDevMax') & size(params.QeDevMax,2)==1
    params.QeDevMax     = params.QeDevMax*ones(1,params.T_loss);
elseif isfield(params,'QeDevMax') & size(params.QeDevMax,2)<params.T_loss
    params.QeDevMax(size(params.QeDevMax,2)+1:params.T_loss) = 100;
end

%% Plotting
if ~isfield(params,'PlotConstraints')
    params.PlotConstraints        = 'ON';     % Choice: ON, OFF
end

%% Misc
if ~isfield(params,'DynareOptions')
    params.DynareOptions = ' ';
end

%% Algorithm
if ~isfield(params,'MaxIter')
    params.MaxIter      = 1e4;
end
if ~isfield(params,'Update')
    params.Update       = .1 + .4*0;
end
if ~isfield(params,'Crit')
    params.Crit         = 1e-3;
end
if ~isfield(params,'Smoothing')
    params.Smoothing    = 1;
end

