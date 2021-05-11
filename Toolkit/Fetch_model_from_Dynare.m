function [m, Shocks_u, t0_inherited] = Fetch_model_from_Dynare(M_,oo_,PlannerDiscountFactor,PolicyInstrumentsAndShocks,LossVars)
%--------------------------------------------------------------------------
% Fetch_model_from_Dynare    -- Assuming that a model file has already been
%                               processed by Dynare, extract the linear
%                               model from Dynare's own data structures.
%
% INPUTS:
%   -> M_                    : The Dynare model struct
%   -> oo_                   : The Dynare computational output struct
%   -> InstrVarsEqnsAndShockNames : A 3-column cell array with a separate
%      (optional)              row for each policy instrument of the form:
%                             [ PolicyVariable , EquationNo , PolicyShock ]
%                              > PolicyVariable: Name of the variable used
%                                                as an instrument.
%                              > EquationNo: The equation describing the
%                                            current policy rule for the
%                                            instrument.
%                              > PolicyShock: Name of the shock that is
%                                             applied in the policy rule.
%   -> LossVars              : A cell vector of variables that might be
%      (optional)              used in a loss function.  We find the index
%                              for each.
% OUTPUTS:
%   <- m                     : A model file for further processing
%
% NOTES:
%   To understand this function, see the documentation in the Dynare
%   Reference Manual entries on "FILENAME_dynamic.m" and "Typology and
%   ordering of variables".
%
%   Acknowledgements: Modified from code provided by J. Barrdear
%--------------------------------------------------------------------------

% How big is the model?
m.nx = M_.endo_nbr;
m.ns = M_.exo_nbr;

% ---------------------------------------------------------------------
% Extract the (linearised) structural model.  This is obtained by
% calling the Dynare-created function -- FILENAME_dynamic -- where
% FILENAME is the name of the model file passed to Dynare in the first
% place.  This returns the Jacobian of the model (and, if requested,
% the Hessian and third-order derivatives).

% Extract the Jacobian (code copied from model_diagnostics.m in Dynare)
klen = M_.maximum_lag + M_.maximum_lead + 1;
exo_simul = [repmat(oo_.exo_steady_state',klen,1) repmat(oo_.exo_det_steady_state',klen,1)];
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1;
z = repmat(oo_.dr.ys,1,klen);

[~,model_jacobian] = feval([M_.fname '_dynamic'],z(iyr0),exo_simul, ...
    M_.params, oo_.dr.ys, it_);

% Break it out into lag, current, lead and shock loading matrices
m.H_B = zeros(m.nx,m.nx);
m.H_C = zeros(m.nx,m.nx);
m.H_F = zeros(m.nx,m.nx);

if M_.maximum_lag
    LLIrow = 1;
    [~,row,col]=find(M_.lead_lag_incidence(LLIrow,:));
    m.H_B(:,row) = model_jacobian(:,col);
end

LLIrow = 2;
[~,row,col]=find(M_.lead_lag_incidence(LLIrow,:));
m.H_C(:,row) = model_jacobian(:,col);

if M_.maximum_lead
    LLIrow = 3;
    [~,row,col]=find(M_.lead_lag_incidence(LLIrow,:));
    m.H_F(:,row) = model_jacobian(:,col);
end

m.PSI = -model_jacobian(:,nnz(M_.lead_lag_incidence')+1:end);


% ---------------------------------------------------------------------
% Extract the (linearised) model solution.  In theory, since we have
% the structural model, we could just solve it ourselves.  But Dynare
% has already solved it for us, so just pull it out.

% The solution is of the form:  x_{t} = B x_{t-1} + PHI u_{t} where
% x_{t} are deviations from steady state.  This reordering of the state
% variables comes from Ingvar Strid (Riksbank).

%     m.B(:,oo_.dr.state_var) = oo_.dr.ghx(oo_.dr.inv_order_var,:);
%     m.PHI = oo_.dr.ghu(oo_.dr.inv_order_var,:);

m.B   = zeros(m.nx,m.nx);
m.PHI = zeros(m.nx,m.ns);

% The number of 'Effective' or 'real' state variables (economic states)
numberOfDynareStateVars = size(oo_.dr.state_var,2);

stateOrderNew = zeros(numberOfDynareStateVars,1);
iter = 0;
for i=1:m.nx
    stateVarInOrderVar = oo_.dr.order_var(i,1);
    for j=1:numberOfDynareStateVars
        stateVarInStateVar = oo_.dr.state_var(1,j);
        if stateVarInOrderVar == stateVarInStateVar
            iter = iter + 1;
            stateOrderNew(iter,1) = stateVarInOrderVar;
        end
    end
end

m.B(:,stateOrderNew) = oo_.dr.ghx(oo_.dr.inv_order_var,:);
m.PHI = oo_.dr.ghu(oo_.dr.inv_order_var,:);

% The F matrix is used to apply anticipated shocks.  Dynare does not
% produce it, so we create it here.
m.F = -(m.H_C + m.H_F*m.B)\m.H_F;

m.SS = oo_.dr.ys;

%----------------------------------------------------------------------
% Figure out which variables and shocks are policy instruments and
% which are featured in the policymaker's loss function.
m.VarNames   = cellstr(M_.endo_names);
m.ShockNames = cellstr(M_.exo_names);
m.ParamNames = cellstr(M_.param_names);
m.Params     = M_.params;
if nargin >= 3
    if isnumeric(PlannerDiscountFactor)
        m.DiscountFactor = PlannerDiscountFactor;
    elseif sum(strcmp(PlannerDiscountFactor,m.ParamNames))==1
        m.DiscountFactor = m.Params(strcmp(PlannerDiscountFactor,m.ParamNames));
    else
        error('The discount factor params.PlannerDiscountFactor needs to be provided as a numeric value or alternatively as a character array indicating the corresponding parameter name in the existing mod file. The parameter name you indicated could not be found in the mod file.')
    end
end

if nargin >= 4
    m.PolicyInstrumentsAndShocks = PolicyInstrumentsAndShocks;
    
    m.ni = size(PolicyInstrumentsAndShocks,1);
    m.InstrVarIndexes = nan(1,m.ni);
    m.InstrShkIndexes = nan(1,m.ni);
    
    for i = 1:m.ni
        var_name = PolicyInstrumentsAndShocks{i,1};
        shk_name = PolicyInstrumentsAndShocks{i,2};
        
        m.InstrVarIndexes(i) = find(strcmp(var_name,m.VarNames));
        m.InstrShkIndexes(i) = find(strcmp(shk_name,m.ShockNames));
    end
end

if nargin >= 5
    m.nl = size(LossVars,1);
    m.LossVarIndexes = nan(1,m.nl);
    
    for i = 1:m.nl
        var_name = LossVars{i};
        m.LossVarIndexes(i) = find(strcmp(var_name,m.VarNames));
    end
end


%----------------------------------------------------------------------
% If oo_ includes smoothed shocks, bundle them together.
if isfield(oo_,'SmoothedShocks')
    for i = 1:m.ns
        shks = getfield(oo_.SmoothedShocks,char(m.ShockNames(i)));
        if i == 1
            T = length(shks);
            Shocks_u = nan(m.ns,T);
        end
        Shocks_u(i,:) = shks';
    end
    
    for i = 1:m.nx
        x = getfield(oo_.SmoothedVariables,char(m.VarNames(i)));
        if i == 1
            T = length(x);
            SmoothedDeviations = nan(m.ns,T);
        end
        SmoothedDeviations(i,:) = x' - m.SS(i);
    end
    
    t0_inherited = SmoothedDeviations(:,1) - m.PHI*Shocks_u(:,1);
end


end