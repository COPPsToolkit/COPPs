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
params.ModelDirectory         = '..\Models\SW07';
params.ModelFilename          = 'Smets_Wouters_2007_GB09';
params.PlannerDiscountFactor  = 'beta_ss';

% Target variables and policy variables
params.LossFunctionVariables       = {'obs_pinf_4q';'ygap';'dr'};
params.PolicyInstrumentsAndShocks  =  {'r' ,'em'};

% Timeseries
params.PastPeriods                 = 219;  % 2009 March

%% ########################################################################
%% ##################        Optional user input         ##################
%% ########################################################################
%%
params.LossFunctionWeights = [1,0.25*1,4];

params.plotting.VarsToPlot = {...
    'obs_pinf_4q' , [ -1 ,  4 ] , 'Inflation (annual, P.P.)'     ; ...
    'ygap'        , [ -9 ,  6 ] , 'Output gap (P.P.)'    ; ...
    'obs_r_ann'   , [ -8 ,  7] , 'Interest rate (annual, P.P.)' ; ...
    'qeobs'       , [ -5 ,  50 ] , 'Asset holdings (% of GDP) '    ; ...
    };

params.T_full              = 78; % How far into the future to roll the forecast forward.
params.plotting.first_date = datetime(1954,7,01); % Set the first date of the series
params.plotting.freq       = 'quarter'; % Frequency of data (month/quarter/year)
params.plotting.NoPlot     = 1;

%% ########################################################################
%% #########        No user input required below this line        #########
%% ########################################################################
%% Set non-user defined values to default
params = SetRemainingParameters(params);

%% Solve for the policy projections
projections = SolveCOPPs(params);

%% Plot everything
PlotCOPPs(projections, params);
