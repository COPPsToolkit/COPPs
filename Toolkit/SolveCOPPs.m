function [results] = SolveCOPPs(params)
%--------------------------------------------------------------------------
%
% NOTES:
% This file loads the model and baseline forecasts from Dynare and calls
% the function that computes constrained optimal policy paths,
% CalculateCOPPs.
%
% Authors: O. de Groot, F. Mazelis, R. Motto, A. Ristiniemi, 2019-2021
% Acknowledgements: J. Barrdear
%--------------------------------------------------------------------------

COPPs

tic

%% Processing

orig_directory = pwd;
cd(params.ModelDirectory);

if strcmp('.mod',params.ModelFilename(end-4:end)) || strcmp('.MOD',params.ModelFilename(end-4:end))
    filename_stub = params.ModelFilename(1:end-4);
else
    filename_stub = params.ModelFilename;
end

%% Call Dynare to estimate and fetch the model using simple rules for policy
display('Call Dynare to solve the model in its baseline form')
dynare(filename_stub,'-DLOAD_ESTIMATED_PARAMS','-DSMOOTHER',params.DynareOptions);
load([filename_stub,'_results.mat'],'-mat','M_','oo_');

[m_sr, Shocks_u, t0_inherited] = Fetch_model_from_Dynare(M_,oo_,params.PlannerDiscountFactor,params.PolicyInstrumentsAndShocks,params.LossFunctionVariables);
d_sr = Build_d(m_sr,params.PastPeriods,params.T_full,Shocks_u,t0_inherited);
Loss_sr = (m_sr.DiscountFactor.^(0:1:(params.T_full-1)))*(d_sr.Forecast.Deviations(m_sr.LossVarIndexes,:) .^ 2)'*params.LossFunctionWeights';

results.simple_rules.model = m_sr;
results.simple_rules.data  = d_sr;
results.simple_rules.loss  = Loss_sr;

%% COPPs
display(' ')
display('Solving COPPs...')
display(' ')

m = m_sr;
d = d_sr;

[m.aInverter,m.aP,m.aQ,m.aR,m.aC] = BuildInversionMatrices(m,params);

ref_x = d.Forecast.Deviations;
margPath_a = zeros([length(m.InstrVarIndexes) params.T_loss]);
x0 = zeros(m.nx,1);

% Constrained optimisation
[x_a,~,shocks_a,exitflag, Lopt] = CalculateCOPPs(m,margPath_a,{},x0,params,ref_x);

marg_x(:,1:params.T_loss) = x_a;

for h = (params.T_loss+1):params.T_full
    marg_x(:,h) = m.B*marg_x(:,h-1);
end

% Anticipated shocks
shocks_a = [shocks_a , zeros(m.ns,params.T_full-params.T_loss)];

% Projections
Forecasts = d.Forecast;
Forecasts.Deviations = ref_x + marg_x;
Forecasts.Shocks_a = [shocks_a , zeros(m.ns,params.T_full - size(shocks_a,2))];

results.COPPs.data = d_sr;
results.COPPs.data.Forecast = Forecasts;

%% Display output

DisplayOutput(params, exitflag, Lopt)

toc

cd(orig_directory)

end
