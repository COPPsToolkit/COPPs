
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
addpath('Toolkit');

%% ########################################################################
%% ##################        Required user input         ##################
%% ########################################################################
%% 

% Model details
params.ModelDirectory         = 'Models\SW07';
params.ModelFilename          = 'Smets_Wouters_2007_GB09_mod210206';
params.PlannerDiscountFactor  = 'beta_ss';

% Target variables and policy variables
params.LossFunctionVariables       = {'obs_pinf_4q';'ygap';'dr'};
params.PolicyInstrumentsAndShocks  =  {'r' ,'em'};

% Timeseries
params.PastPeriods          = 219;  % 2009 March
pramas.T_full               = 78; % How far into the future to roll the forecast forward.
params.plotting.first_date  = datetime(1954,7,01); % Set the first date of the series
params.plotting.freq        = 'quarter'; % Frequency of data (month/quarter/year)


%% ########################################################################
%% ##################        Optional user input         ##################
%% ########################################################################
%% 

%% ########################################################################
%% #########        No user input required below this line        #########
%% ########################################################################
%% Set non-user defined values to default
params = SetRemainingParameters(params);

%% Solve for the policy projections
projections = SolveCOPPs(params);

%% Plot everything
PlotCOPPs(projections, params);
