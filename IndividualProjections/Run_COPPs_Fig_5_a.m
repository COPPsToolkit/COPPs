%% ########################################################################
%% ################ Constrained Optimal Policy Projections ################
%% ########################################################################
%
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
% This  version: 1.0
%__________________________________________________________________________
%% ########################################################################
%% ##################             Preamble               ##################
%% ########################################################################
%% 
clear;
clc;

% Add toolkit to matlab path 
addpath('..\Toolkit');

%% ########################################################################
%% ##################        Required user input         ##################
%% ########################################################################
%% 

% Model details
params.ModelDirectory         = '..\Models\GaliNK3';
params.ModelFilename          = 'GaliNK3';
params.PlannerDiscountFactor  = 'beta';
params.RateMin                = 0;        % ELB in annual policy rate
params.PolicyRateStSt         = 0;

% Target variables and policy variables
params.LossFunctionVariables       = {'pinf';'og'};
params.PolicyInstrumentsAndShocks  = {'r','er'};

% Timeseries
params.PastPeriods                 = 50;  

%% ########################################################################
%% ##################        Optional user input         ##################
%% ########################################################################
%%
% Forward guidance
params.alpha = 1;

params.T_loss   = 60;
params.H_policy = 1;

kappa   = 0.1275;
epsilon = 6;
alpha_x = kappa/epsilon;

params.LossFunctionWeights         = [1,alpha_x];

params.plotting.NoPlot     = 1;
params.plotting.VarsToPlot = {...
    'pinf_ann',[-3,1]  ,'Inflation'  ;...
    'og'      ,[-1,3]  ,'Output gap' ;...
    'r_ann'   ,[-10,10],'Policy rate';...
    };


%% ########################################################################
%% #########        No user input required below this line        #########
%% ########################################################################
%% Set non-user defined values to default
params = SetRemainingParameters(params);

%% Solve for the policy projections
projections = SolveCOPPs(params);

%% Plot everything
params.plotting.P_past   = 4;
params.plotting.P_future = 13;
PlotCOPPs(projections, params);
