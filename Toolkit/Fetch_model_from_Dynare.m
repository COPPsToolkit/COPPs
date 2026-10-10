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
%   Requires a version of Dynare that generates the sparse representation
%   of the model (+FILENAME/dynamic_g1.m), i.e. the unstable version.
%
%   Acknowledgements: Modified from code provided by J. Barrdear
%--------------------------------------------------------------------------

% How big is the model?
m.nx = M_.endo_nbr;
m.ns = M_.exo_nbr;

% ---------------------------------------------------------------------
% Extract the (linearised) structural model.  This is obtained by
% calling the Dynare-created function -- +FILENAME/dynamic_g1.m -- where
% FILENAME is the name of the model file passed to Dynare in the first
% place.  This returns the (sparse) Jacobian of the model with respect to
% [y_{t-1}; y_t; y_{t+1}; u_t], endogenous variables in declaration order.

n  = M_.endo_nbr;
g1 = feval([M_.fname '.dynamic_g1'], repmat(oo_.dr.ys,3,1), ...
    [oo_.exo_steady_state; oo_.exo_det_steady_state], M_.params, oo_.dr.ys, ...
    M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);
g1 = full(g1);

% Break it out into lag, current, lead and shock loading matrices
m.H_B = g1(:,1:n);
m.H_C = g1(:,n+1:2*n);
m.H_F = g1(:,2*n+1:3*n);

m.PSI = -g1(:,3*n+(1:m.ns));


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