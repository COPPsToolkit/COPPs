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
params.ModelDirectory         = '..\Models\SimsWu';
params.ModelFilename          = 'SWu20';
params.PlannerDiscountFactor  = 'beta_ss';

% Target variables and policy variables
params.LossFunctionVariables       = {'obs_pinf_4q';'ygap';'dr';'bcb';'dbcb'};
params.PolicyInstrumentsAndShocks  =  {'r' ,'er'; 'qeobs','eb'};

% Timeseries
params.PastPeriods         = 219;  % 2009 March

%% ########################################################################
%% ##################        Optional user input         ##################
%% ########################################################################
%%
params.LossFunctionWeights  = [1,0.25,4,0.57,4.58];
params.DynareOptions        = ' -DFlatPC';
params.Constraint           = 'ON';
params.RateMin              = 0.1;        % ELB in annual policy rate
params.PolicyRateStSt       = 'ss_r_ann'; % Parameter name of annual StSt policy rate
params.QeDevMax(1:20)       = [3:3:60];   % qe: Allow max increase of 3PP per quarter
params.RateDevMin           = 0; 
params.RateDevMax           = 0; 

params.plotting.VarsToPlot = {...
    'obs_pinf_4q' , [ -1 ,  4 ]  , 'Inflation (annual, P.P.)'     ; ...
    'ygap'        , [ -9 ,  6 ]  , 'Output gap (P.P.)'            ; ...
    'obs_r_ann'   , [ -8 ,  7]   , 'Interest rate (annual, P.P.)' ; ...
    'qeobs'       , [ -5 ,  50 ] , 'Asset holdings (% of GDP) '   ; ...
    };

paraas.T_full              = 78; % How far into the future to roll the forecast forward.
params.plotting.first_date = datetime(1954,7,01); % Set the first date of the data series
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
